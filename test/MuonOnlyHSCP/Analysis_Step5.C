#include <string>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "tdrstyle.C"
#include "TProfile.h"

#include "Analysis_CommonFunction.h"
#include "Analysis_Global.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"
#include "Analysis_Samples.h"

using namespace std;

/////////////////////////// FUNCTION DECLARATION /////////////////////////////

void CutFlow(string InputPattern, unsigned int CutIndex=0);
void SelectionPlot (string InputPattern, unsigned int CutIndex);
void PredictionAndControlPlot(string InputPattern, string Suffix);
void CosmicBackgroundSystematic(string InputPattern);
void CollisionBackgroundSystematic(string InputPattern);

std::vector<stSignal> signals;

string LegendTitle;

/////////////////////////// CODE PARAMETERS /////////////////////////////

void Analysis_Step5()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);


   GetSignalDefinition(signals);

   string InputDir;				unsigned int CutIndex;
   std::vector<string> Legends;                 std::vector<string> Inputs;

   InputDir = "Results/Eta21/PtMin80/";   
   //CutIndex = 489;
   CutIndex = 150;
   SelectionPlot(InputDir, CutIndex);
   PredictionAndControlPlot(InputDir, "Data/");
   CutFlow(InputDir, CutIndex);
   CosmicBackgroundSystematic(InputDir);
   CollisionBackgroundSystematic(InputDir);
   return;
}

void CutFlow(string InputPattern, unsigned int CutIndex){
  //string Input     = InputPattern + "Histos.root";
   string SavePath  = InputPattern + "/CutFlow/";
   MakeDirectories(SavePath);

   string Input = "/uscms_data/d2/farrell3/WorkArea/SignalSA2012/CMSSW_4_2_8_patch6/src/SUSYBSMAnalysis/HSCP/test/MuonOnly/Results/Eta21/PtMin70/Histos.root";
   TFile* InputFile = new TFile(Input.c_str());
   TFile* InputFileData = new TFile((InputPattern + "Histos_Data.root").c_str());
   TFile* InputFileCosmic = new TFile((InputPattern + "Histos_Cosmic.root").c_str());

   TH1D*  HCuts_Pt       = (TH1D*)GetObjectFromPath(InputFileData, "HCuts_Pt");
   TH1D*  HCuts_TOF      = (TH1D*)GetObjectFromPath(InputFileData, "HCuts_TOF");

   char Buffer[1024]; sprintf(Buffer,"%s/CutFlow_%03i_Pt%03.0f_TOF%04.3f.txt",SavePath.c_str(),CutIndex,HCuts_Pt->GetBinContent(CutIndex+1)/files,HCuts_TOF->GetBinContent(CutIndex+1)/files);
    FILE* pFile = fopen(Buffer,"w");

    stPlots DataPlots, DataPlots_Track, DataPlots_NoTrack, DataPlots_Control;

    stPlots_InitFromFile(InputFile, DataPlots,"Data", InputFileData);
    stPlots_InitFromFile(InputFile, DataPlots_Track,"Data_Track", InputFileData);
    stPlots_InitFromFile(InputFile, DataPlots_NoTrack,"Data_NoTrack", InputFileData);
    stPlots_InitFromFile(InputFile, DataPlots_Control,"Data_Control", InputFileData);

    stPlots_Dump(DataPlots, pFile, CutIndex);
    stPlots_Dump(DataPlots_Track, pFile, CutIndex);
    stPlots_Dump(DataPlots_NoTrack, pFile, CutIndex);
    stPlots_Dump(DataPlots_Control, pFile, CutIndex);

    stPlots_Clear(DataPlots);
    stPlots_Clear(DataPlots_Track);
    stPlots_Clear(DataPlots_NoTrack);
    stPlots_Clear(DataPlots_Control);

    stPlots CosmicPlots, CosmicPlots_Track, CosmicPlots_NoTrack, CosmicPlots_Control;
    stPlots_InitFromFile(InputFile, CosmicPlots,"Cosmic", InputFileCosmic);
    stPlots_InitFromFile(InputFile, CosmicPlots_Track,"Cosmic_Track", InputFileCosmic);
    stPlots_InitFromFile(InputFile, CosmicPlots_NoTrack,"Cosmic_NoTrack", InputFileCosmic);
    stPlots_InitFromFile(InputFile, CosmicPlots_Control,"Cosmic_Control", InputFileCosmic);
    stPlots_Dump(CosmicPlots, pFile, CutIndex);
    stPlots_Dump(CosmicPlots_Track, pFile, CutIndex);
    stPlots_Dump(CosmicPlots_NoTrack, pFile, CutIndex);
    stPlots_Dump(CosmicPlots_Control, pFile, CutIndex);
    stPlots_Clear(CosmicPlots);
    stPlots_Clear(CosmicPlots_Track);
    stPlots_Clear(CosmicPlots_NoTrack);
    stPlots_Clear(CosmicPlots_Control);

    for(unsigned int s=0;s<signals.size();s++){
       if(!signals[s].MakePlot)continue;
       stPlots SignPlots;
       stPlots_InitFromFile(InputFile, SignPlots,signals[s].Name, InputFile);
       stPlots_Dump(SignPlots, pFile, CutIndex);       
       stPlots_Clear(SignPlots);
    }

    fclose(pFile);

}


void SelectionPlot(string InputPattern, unsigned int CutIndex){

  LegendTitle = "SA Only";

  string Input = "/uscms_data/d2/farrell3/WorkArea/SignalSA2012/CMSSW_4_2_8_patch6/src/SUSYBSMAnalysis/HSCP/test/MuonOnly/Results/Eta21/PtMin80/Histos.root";
  //string Input     = "/uscms/home/farrell3/nobackup/WorkArea/EDMMuonOnly/src/SUSYBSMAnalysis/HSCP/test/MuonOnly/" + InputPattern + "Histos.root";
   string SavePath  = InputPattern;
   MakeDirectories(SavePath);

   TFile* InputFile = new TFile(Input.c_str());
   TFile* InputFileData = new TFile((InputPattern + "Histos_Data.root").c_str());
   TFile* InputFileCosmic = new TFile((InputPattern + "Histos_Cosmic.root").c_str()); 

   stPlots DataPlots, DataPlotsNoTrack, DataPlotsTrack, DataPlotsControl, CosmicPlots, CosmicPlotsControl, SignPlots[signals.size()];
   stPlots_InitFromFile(InputFile, DataPlotsTrack,"Data_Track", InputFileData);
   stPlots_InitFromFile(InputFile, DataPlotsNoTrack,"Data_NoTrack", InputFileData);
   stPlots_InitFromFile(InputFile, DataPlotsControl,"Data_Control", InputFileData);
   stPlots_InitFromFile(InputFile, DataPlots,"Data", InputFileData);

   //stPlots_Draw(DataPlotsTrack, SavePath + "/Selection_Data_Track", LegendTitle, CutIndex);
   stPlots_Draw(DataPlotsControl, SavePath + "/Selection_Data_Control", LegendTitle, CutIndex);
   stPlots_Draw(DataPlotsTrack, SavePath + "/Selection_Data", LegendTitle, CutIndex);

   stPlots_InitFromFile(InputFile, CosmicPlots,"Cosmic", InputFileCosmic);
   stPlots_Draw(CosmicPlots, SavePath + "/Selection_Cosmic", LegendTitle, CutIndex);

   for(unsigned int s=0;s<signals.size();s++){
     stPlots_InitFromFile(InputFile, SignPlots[s],signals[s].Name, InputFile);
      if(!signals[s].MakePlot)continue;
      //if (signals[s].Name=="Gluino1000") stPlots_Draw(SignPlots[s], SavePath + "/Selection_" +  signals[s].Name, LegendTitle, CutIndex);
      if (signals[s].Name=="Gluino500") stPlots_Draw(SignPlots[s], SavePath + "/Selection_" +  signals[s].Name, LegendTitle, CutIndex);
      if (signals[s].Name=="GMStau494") stPlots_Draw(SignPlots[s], SavePath + "/Selection_" +  signals[s].Name, LegendTitle, CutIndex);
      //if (signals[s].Name=="GMStau126") stPlots_Draw(SignPlots[s], SavePath + "/Selection_" +  signals[s].Name, LegendTitle, CutIndex);
   }
   //Stplots_Draw(DataPlotsControl, SavePath + "/Selection_DataControl", LegendTitle, CutIndex);

   //stPlots_DrawComparison(SavePath + "/Selection_Comp_Data" , LegendTitle, GluinoCutIndex, &DataPlotsTrack);

   stPlots_DrawComparison(SavePath + "/Selection_Comp_Data" , LegendTitle, CutIndex, &DataPlotsNoTrack, &CosmicPlots);
   stPlots_DrawComparison(SavePath + "/Selection_Comp_Gluino" , LegendTitle, CutIndex, &DataPlots, &CosmicPlots, &SignPlots[2], &SignPlots[7]);
   stPlots_DrawComparison(SavePath + "/Selection_Comp_Signal" , LegendTitle, CutIndex, &SignPlots[2], &SignPlots[18], &SignPlots[23]);

   //stPlots_DrawComparison(SavePath + "/Selection_Check" , LegendTitle, GluinoCutIndex, &DataPlotsNoTrack, &CosmicPlots);
   //stPlots_DrawComparison(SavePath + "/Selection_Comp_Control" , LegendTitle, GluinoCutIndex, &DataPlots, &DataPlotsControl, &SignPlots[0], &SignPlots[7]);
   return;
}



 //////////////////////////////////////////////////     CREATE PLOTS OF CONTROLS AND PREDICTION

void PredictionAndControlPlot(string InputPattern, string Suffix){
   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   LegendTitle = "SA Only";
   string Input     = InputPattern + "Histos_Data.root";
   string SavePath  = InputPattern;
   MakeDirectories(SavePath);

   TFile* InputFile = new TFile(Input.c_str());

   TH1D* CtrlPt_S1_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlPt_S1_TOF");
   TH1D* CtrlPt_S2_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlPt_S2_TOF");
   TH1D* CtrlPt_S3_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlPt_S3_TOF");
   TH1D* CtrlPt_S4_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlPt_S4_TOF");

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlPt_S1_TOF->Integral()>0)CtrlPt_S1_TOF->Scale(1/CtrlPt_S1_TOF->Integral());
   if(CtrlPt_S2_TOF->Integral()>0)CtrlPt_S2_TOF->Scale(1/CtrlPt_S2_TOF->Integral());
   if(CtrlPt_S3_TOF->Integral()>0)CtrlPt_S3_TOF->Scale(1/CtrlPt_S3_TOF->Integral());
   if(CtrlPt_S4_TOF->Integral()>0)CtrlPt_S4_TOF->Scale(1/CtrlPt_S4_TOF->Integral());
   Histos[0] = CtrlPt_S1_TOF;                    legend.push_back(" 70 <p_{T}< 90GeV");
   Histos[1] = CtrlPt_S2_TOF;                    legend.push_back(" 90<p_{T}< 130 GeV");
   Histos[2] = CtrlPt_S3_TOF;                    legend.push_back(" 130<p_{T}< 250 GeV");
   Histos[3] = CtrlPt_S4_TOF;                    legend.push_back(" 250<p_{T}");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", "arbitrary units", 0,2, 0,0); 
   DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlPt_TOFSpectrum");
   delete c1;


   TH1D* CtrlCen_Pt_S1_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_Pt_S1_TOF");
   TH1D* CtrlCen_Pt_S2_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_Pt_S2_TOF");
   TH1D* CtrlCen_Pt_S3_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_Pt_S3_TOF");
   TH1D* CtrlCen_Pt_S4_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_Pt_S4_TOF");

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlCen_Pt_S1_TOF->Integral()>0)CtrlCen_Pt_S1_TOF->Scale(1/CtrlCen_Pt_S1_TOF->Integral());
   if(CtrlCen_Pt_S2_TOF->Integral()>0)CtrlCen_Pt_S2_TOF->Scale(1/CtrlCen_Pt_S2_TOF->Integral());
   if(CtrlCen_Pt_S3_TOF->Integral()>0)CtrlCen_Pt_S3_TOF->Scale(1/CtrlCen_Pt_S3_TOF->Integral());
   if(CtrlCen_Pt_S4_TOF->Integral()>0)CtrlCen_Pt_S4_TOF->Scale(1/CtrlCen_Pt_S4_TOF->Integral());
   Histos[0] = CtrlCen_Pt_S1_TOF;                    legend.push_back(" 70 <p_{T}< 90GeV");
   Histos[1] = CtrlCen_Pt_S2_TOF;                    legend.push_back(" 90<p_{T}< 130 GeV");
   Histos[2] = CtrlCen_Pt_S3_TOF;                    legend.push_back(" 130<p_{T}< 250 GeV");
   Histos[3] = CtrlCen_Pt_S4_TOF;                    legend.push_back(" 250<p_{T}");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", "arbitrary units", 0,2, 0,0);
   DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlCen_Pt_TOFSpectrum");
   delete c1;

   TH1D* CtrlFor_Pt_S1_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_Pt_S1_TOF");
   TH1D* CtrlFor_Pt_S2_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_Pt_S2_TOF");
   TH1D* CtrlFor_Pt_S3_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_Pt_S3_TOF");
   TH1D* CtrlFor_Pt_S4_TOF        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_Pt_S4_TOF");

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlFor_Pt_S1_TOF->Integral()>0)CtrlFor_Pt_S1_TOF->Scale(1/CtrlFor_Pt_S1_TOF->Integral());
   if(CtrlFor_Pt_S2_TOF->Integral()>0)CtrlFor_Pt_S2_TOF->Scale(1/CtrlFor_Pt_S2_TOF->Integral());
   if(CtrlFor_Pt_S3_TOF->Integral()>0)CtrlFor_Pt_S3_TOF->Scale(1/CtrlFor_Pt_S3_TOF->Integral());
   if(CtrlFor_Pt_S4_TOF->Integral()>0)CtrlFor_Pt_S4_TOF->Scale(1/CtrlFor_Pt_S4_TOF->Integral());
   Histos[0] = CtrlFor_Pt_S1_TOF;                    legend.push_back(" 70 <p_{T}< 90GeV");
   Histos[1] = CtrlFor_Pt_S2_TOF;                    legend.push_back(" 90<p_{T}< 130 GeV");
   Histos[2] = CtrlFor_Pt_S3_TOF;                    legend.push_back(" 130<p_{T}< 250 GeV");
   Histos[3] = CtrlFor_Pt_S4_TOF;                    legend.push_back(" 250<p_{T}");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", "arbitrary units", 0,2, 0,0);
   DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlFor_Pt_TOFSpectrum");
   delete c1;


   TH1D* CtrlTOF_S1_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlTOF_S1_Pt"); 
   TH1D* CtrlTOF_S2_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlTOF_S2_Pt");
   TH1D* CtrlTOF_S3_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlTOF_S3_Pt");
   TH1D* CtrlTOF_S4_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlTOF_S4_Pt");

   CtrlTOF_S1_Pt->Rebin(4);  
   CtrlTOF_S2_Pt->Rebin(4);
   CtrlTOF_S3_Pt->Rebin(4);
   CtrlTOF_S4_Pt->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlTOF_S1_Pt->Integral()>0)CtrlTOF_S1_Pt->Scale(1/CtrlTOF_S1_Pt->Integral(0, CtrlTOF_S1_Pt->GetNbinsX()+1)); 
   if(CtrlTOF_S2_Pt->Integral()>0)CtrlTOF_S2_Pt->Scale(1/CtrlTOF_S2_Pt->Integral(0, CtrlTOF_S2_Pt->GetNbinsX()+1));
   if(CtrlTOF_S3_Pt->Integral()>0)CtrlTOF_S3_Pt->Scale(1/CtrlTOF_S3_Pt->Integral(0, CtrlTOF_S3_Pt->GetNbinsX()+1));
   if(CtrlTOF_S4_Pt->Integral()>0)CtrlTOF_S4_Pt->Scale(1/CtrlTOF_S4_Pt->Integral(0, CtrlTOF_S4_Pt->GetNbinsX()+1));
   Histos[0] = CtrlTOF_S2_Pt;                    legend.push_back("1.0<TOF<1.1");
   Histos[1] = CtrlTOF_S3_Pt;                    legend.push_back("0.9<TOF<1.0");
   Histos[2] = CtrlTOF_S4_Pt;                    legend.push_back("TOF<0.9");
   Histos[3] = CtrlTOF_S4_Pt;                    legend.push_back("1.1<TOF");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "P_{t}", "arbitrary units", 80,500, 0,0, false, true, true); 
   //DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlTOF_PtSpectrum");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlTOF_PtSpectrumLog");
   delete c1;

   TH1D* CtrlCen_TOF_S1_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_TOF_S1_Pt");
   TH1D* CtrlCen_TOF_S2_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_TOF_S2_Pt");
   TH1D* CtrlCen_TOF_S3_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_TOF_S3_Pt");
   TH1D* CtrlCen_TOF_S4_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlCen_TOF_S4_Pt");

   CtrlCen_TOF_S1_Pt->Rebin(4);
   CtrlCen_TOF_S2_Pt->Rebin(4);
   CtrlCen_TOF_S3_Pt->Rebin(4);
   CtrlCen_TOF_S4_Pt->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlCen_TOF_S1_Pt->Integral()>0)CtrlCen_TOF_S1_Pt->Scale(1/CtrlCen_TOF_S1_Pt->Integral());
   if(CtrlCen_TOF_S2_Pt->Integral()>0)CtrlCen_TOF_S2_Pt->Scale(1/CtrlCen_TOF_S2_Pt->Integral());
   if(CtrlCen_TOF_S3_Pt->Integral()>0)CtrlCen_TOF_S3_Pt->Scale(1/CtrlCen_TOF_S3_Pt->Integral());
   if(CtrlCen_TOF_S4_Pt->Integral()>0)CtrlCen_TOF_S4_Pt->Scale(1/CtrlCen_TOF_S4_Pt->Integral());
   Histos[3] = CtrlCen_TOF_S1_Pt;                    legend.push_back("1.1<TOF");
   Histos[2] = CtrlCen_TOF_S2_Pt;                    legend.push_back("1.0<TOF<1.1");
   Histos[1] = CtrlCen_TOF_S3_Pt;                    legend.push_back("0.9<TOF<1.0");
   Histos[0] = CtrlCen_TOF_S4_Pt;                    legend.push_back("TOF<0.9");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "P_{t}", "arbitrary units", 80,500, 0,0, false, true, true);
   DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlCen_TOF_PtSpectrum");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlCen_TOF_PtSpectrumLog");
   delete c1;

   TH1D* CtrlFor_TOF_S1_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_TOF_S1_Pt");
   TH1D* CtrlFor_TOF_S2_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_TOF_S2_Pt");
   TH1D* CtrlFor_TOF_S3_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_TOF_S3_Pt");
   TH1D* CtrlFor_TOF_S4_Pt        = (TH1D*)GetObjectFromPath(InputFile, Suffix + "CtrlFor_TOF_S4_Pt");

   CtrlFor_TOF_S1_Pt->Rebin(4);
   CtrlFor_TOF_S2_Pt->Rebin(4);
   CtrlFor_TOF_S3_Pt->Rebin(4);
   CtrlFor_TOF_S4_Pt->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   if(CtrlFor_TOF_S1_Pt->Integral()>0)CtrlFor_TOF_S1_Pt->Scale(1/CtrlFor_TOF_S1_Pt->Integral());
   if(CtrlFor_TOF_S2_Pt->Integral()>0)CtrlFor_TOF_S2_Pt->Scale(1/CtrlFor_TOF_S2_Pt->Integral());
   if(CtrlFor_TOF_S3_Pt->Integral()>0)CtrlFor_TOF_S3_Pt->Scale(1/CtrlFor_TOF_S3_Pt->Integral());
   if(CtrlFor_TOF_S4_Pt->Integral()>0)CtrlFor_TOF_S4_Pt->Scale(1/CtrlFor_TOF_S4_Pt->Integral());
   Histos[3] = CtrlFor_TOF_S1_Pt;                    legend.push_back("1.1<TOF");
   Histos[2] = CtrlFor_TOF_S2_Pt;                    legend.push_back("1.0<TOF<1.1");
   Histos[1] = CtrlFor_TOF_S3_Pt;                    legend.push_back("0.9<TOF<1.0");
   Histos[0] = CtrlFor_TOF_S4_Pt;                    legend.push_back("TOF<0.9");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "P_{t}", "arbitrary units", 80,500, 0,0, false, true, true);
   DrawLegend(Histos,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlFor_TOF_PtSpectrum");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Syst/ControlFor_TOF_PtSpectrumLog");
   delete c1;
}

void CosmicBackgroundSystematic(string InputPattern){

  string SavePath  = InputPattern + "/Syst/";
  MakeDirectories(SavePath);
  TCanvas* c1;
  TH1** Histos = new TH1*[10];
  std::vector<string> legend;

  LegendTitle = "SA Only";

  TFile* InputFileData = new TFile((InputPattern + "Histos_Data.root").c_str());
  TFile* InputFileCosmic = new TFile((InputPattern + "Histos_Cosmic.root").c_str());

  TH1D*  HCuts_Pt       = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_Pt");
  TH1D*  HCuts_TOF      = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_TOF");

  stPlots DataPlots, CosmicPlots;
  stPlots_InitFromFile(InputFileData, DataPlots,"Data_Control", InputFileData);
  stPlots_InitFromFile(InputFileData, CosmicPlots,"Cosmic_Control", InputFileCosmic);

  double PassDz_Cen = CosmicPlots.H_DzCounts_DT->GetBinContent(1);
  double PassDz_For = CosmicPlots.H_DzCounts_CSC->GetBinContent(1);

  TH1F *Pred_TOF100[DzRegions-1];
  TH1F *Pred_TOF110[DzRegions-1];
  TH1F *Pred_TOF120[DzRegions-1];
  TH1F *Pred_TOF130[DzRegions-1];

  TH1F *StatSyst_TOF100 = new TH1F("StatSyst_TOF100", "StatSyst_TOF100", 7, 65, 275);
  TH1F *StatSyst_TOF110 = new TH1F("StatSyst_TOF110", "StatSyst_TOF110", 7, 65, 275);
  TH1F *StatSyst_TOF120 = new TH1F("StatSyst_TOF120", "StatSyst_TOF120", 7, 65, 275);
  TH1F *StatSyst_TOF130 = new TH1F("StatSyst_TOF130", "StatSyst_TOF130", 7, 65, 275);

  TH1F *Stat_TOF100 = new TH1F("Stat_TOF100", "Stat_TOF100", 7, 65, 275);
  TH1F *Stat_TOF110 = new TH1F("Stat_TOF110", "Stat_TOF110", 7, 65, 275);
  TH1F *Stat_TOF120 = new TH1F("Stat_TOF120", "Stat_TOF120", 7, 65, 275);
  TH1F *Stat_TOF130 = new TH1F("Stat_TOF130", "Stat_TOF130", 7, 65, 275);

  TH1F *Syst_TOF100 = new TH1F("Syst_TOF100", "Syst_TOF100", 7, 65, 275);
  TH1F *Syst_TOF110 = new TH1F("Syst_TOF110", "Syst_TOF110", 7, 65, 275);
  TH1F *Syst_TOF120 = new TH1F("Syst_TOF120", "Syst_TOF120", 7, 65, 275);
  TH1F *Syst_TOF130 = new TH1F("Syst_TOF130", "Syst_TOF130", 7, 65, 275);

  for(int Region=1; Region<DzRegions; Region++) {
    string Name="Pred_TOF100_"+RegionNames[Region];
    Pred_TOF100[Region] = new TH1F(Name.c_str(), Name.c_str(), 7, 65, 275);

    Name="Pred_TOF110_"+RegionNames[Region];
    Pred_TOF110[Region] = new TH1F(Name.c_str(), Name.c_str(), 7, 65, 275);

    Name="Pred_TOF120_"+RegionNames[Region];
    Pred_TOF120[Region] = new TH1F(Name.c_str(), Name.c_str(), 7, 65, 275);

    Name="Pred_TOF130_"+RegionNames[Region];
    Pred_TOF130[Region] = new TH1F(Name.c_str(), Name.c_str(), 7, 65, 275);

    double FailDz_Cen = CosmicPlots.H_DzCounts_DT->GetBinContent(Region+1);
    double FailDz_For = CosmicPlots.H_DzCounts_CSC->GetBinContent(Region+1);
    cout << "Fail Cen " << FailDz_Cen << " Pass Cen " << PassDz_Cen << endl;
    double cosmicFraction_Cen = PassDz_Cen/FailDz_Cen;
    double cosmicFractionErr_Cen=sqrt(pow(PassDz_Cen/(FailDz_Cen*FailDz_Cen),2)*FailDz_Cen);

    double cosmicFraction_For = PassDz_For/FailDz_For;
    double cosmicFractionErr_For=sqrt(pow(PassDz_For/(FailDz_For*FailDz_For),2)*FailDz_For);

    for(int Cut=0; Cut<HCuts_Pt->GetNbinsX()+2; Cut++) {
      double D_Control_Cen = DataPlots.H_D_Cen_Syst[Region]->GetBinContent(Cut);
      double NPred_Cosmic_Cen = cosmicFraction_Cen*D_Control_Cen;
      double NPred_Cosmic_Err_Cen = sqrt(cosmicFraction_Cen*cosmicFraction_Cen*D_Control_Cen + cosmicFractionErr_Cen*cosmicFractionErr_Cen*D_Control_Cen*D_Control_Cen);

      double D_Control_For = DataPlots.H_D_For_Syst[Region]->GetBinContent(Cut);
      double NPred_Cosmic_For = cosmicFraction_For*D_Control_For;
      double NPred_Cosmic_Err_For = sqrt(cosmicFraction_For*cosmicFraction_For*D_Control_For + cosmicFractionErr_For*cosmicFractionErr_For*D_Control_For*D_Control_For);

      double NPred = NPred_Cosmic_Cen+NPred_Cosmic_For;
      double NPredErr = sqrt(NPred_Cosmic_Err_Cen*NPred_Cosmic_Err_Cen + NPred_Cosmic_Err_For*NPred_Cosmic_Err_For);

      int Bin=Pred_TOF100[Region]->FindBin(HCuts_Pt->GetBinContent(Cut)/files);

      if(fabs(HCuts_TOF->GetBinContent(Cut)/files-1.00)<0.001) {
	Pred_TOF100[Region]->SetBinContent(Bin, NPred);
	Pred_TOF100[Region]->SetBinError(Bin, NPredErr);
      }

      if(fabs(HCuts_TOF->GetBinContent(Cut)/files-1.1)<0.001) {
        Pred_TOF110[Region]->SetBinContent(Bin, NPred);
        Pred_TOF110[Region]->SetBinError(Bin, NPredErr);
      }

      if(fabs(HCuts_TOF->GetBinContent(Cut)/files-1.2)<0.001) {
        Pred_TOF120[Region]->SetBinContent(Bin, NPred);
        Pred_TOF120[Region]->SetBinError(Bin, NPredErr);
      }

      if(fabs(HCuts_TOF->GetBinContent(Cut)/files-1.3)<0.001) {
        Pred_TOF130[Region]->SetBinContent(Bin, NPred);
        Pred_TOF130[Region]->SetBinError(Bin, NPredErr);
      }
    }
  }

  double UsedRegions=DzRegions-3;
  double endCut=1;

  for(int i=1; i<Pred_TOF100[2]->GetNbinsX()+1; i++) {
    double StatSyst=0, Stat=0, Syst=0, Average=0;
    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      Average+=Pred_TOF100[Region]->GetBinContent(i);
      Stat+=Pred_TOF100[Region]->GetBinError(i)*Pred_TOF100[Region]->GetBinError(i);
    }

    Average=Average/UsedRegions;
    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      StatSyst+=pow(Pred_TOF100[Region]->GetBinContent(i)-Average,2);
    }
    Stat=sqrt(Stat/UsedRegions);
    StatSyst=sqrt(StatSyst/(UsedRegions-1));
    if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

    StatSyst_TOF100->SetBinContent(i, StatSyst/Average);
    Stat_TOF100->SetBinContent(i, Stat/Average);
    Syst_TOF100->SetBinContent(i, Syst/Average);

    StatSyst_TOF100->SetBinError(i, 0.0001);
    Stat_TOF100->SetBinError(i, 0.0001);
    Syst_TOF100->SetBinError(i, 0.0001);
  }

  for(int i=1; i<Pred_TOF110[2]->GetNbinsX()+1; i++) {
    double StatSyst=0, Stat=0, Syst=0, Average=0;
    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      Average+=Pred_TOF110[Region]->GetBinContent(i);
      Stat+=Pred_TOF110[Region]->GetBinError(i)*Pred_TOF110[Region]->GetBinError(i);
    }

    Average=Average/UsedRegions;
    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      StatSyst+=pow(Pred_TOF110[Region]->GetBinContent(i)-Average,2);
    }
    Stat=sqrt(Stat/UsedRegions);
    StatSyst=sqrt(StatSyst/(UsedRegions-1));
    if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

    StatSyst_TOF110->SetBinContent(i, StatSyst/Average);
    Stat_TOF110->SetBinContent(i, Stat/Average);
    Syst_TOF110->SetBinContent(i, Syst/Average);

    StatSyst_TOF110->SetBinError(i, 0.0001);
    Stat_TOF110->SetBinError(i, 0.0001);
    Syst_TOF110->SetBinError(i, 0.0001);
  }

  for(int i=1; i<Pred_TOF120[2]->GetNbinsX()+1; i++) {
    double StatSyst=0, Stat=0, Syst=0, Average=0;

    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      Average+=Pred_TOF120[Region]->GetBinContent(i);
      Stat+=Pred_TOF120[Region]->GetBinError(i)*Pred_TOF120[Region]->GetBinError(i);
    }
    Average=Average/(UsedRegions);

   for(int Region=2; Region<(DzRegions-endCut); Region++) {
      StatSyst+=pow(Pred_TOF120[Region]->GetBinContent(i)-Average,2);
    }
    Stat=sqrt(Stat/UsedRegions);
    StatSyst=sqrt(StatSyst/(UsedRegions-1));
    if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

    StatSyst_TOF120->SetBinContent(i, StatSyst/Average);
    Stat_TOF120->SetBinContent(i, Stat/Average);
    Syst_TOF120->SetBinContent(i, Syst/Average);

    StatSyst_TOF120->SetBinError(i, 0.0001);
    Stat_TOF120->SetBinError(i, 0.0001);
    Syst_TOF120->SetBinError(i, 0.0001);
  }

  for(int i=1; i<Pred_TOF130[2]->GetNbinsX()+1; i++) {
    double StatSyst=0, Stat=0, Syst=0, Average=0;

    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      Average+=Pred_TOF130[Region]->GetBinContent(i);
      Stat+=Pred_TOF130[Region]->GetBinError(i)*Pred_TOF130[Region]->GetBinError(i);
    }
    Average=Average/UsedRegions;

    for(int Region=2; Region<(DzRegions-endCut); Region++) {
      StatSyst+=pow(Pred_TOF130[Region]->GetBinContent(i)-Average,2);
    }

    Stat=sqrt(Stat/UsedRegions);
    StatSyst=sqrt(StatSyst/(UsedRegions-1));
    if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

    StatSyst_TOF130->SetBinContent(i, StatSyst/Average);
    Stat_TOF130->SetBinContent(i, Stat/Average);
    Syst_TOF130->SetBinContent(i, Syst/Average);

    StatSyst_TOF130->SetBinError(i, 0.0001);
    Stat_TOF130->SetBinError(i, 0.0001);
    Syst_TOF130->SetBinError(i, 0.0001);
  }


  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  for(int Region=2; Region<DzRegions; Region++) {
    //Histos[Region-2] = Pred_TOF100[Region];         legend.push_back(LegendNames[Region]);
  }
  Histos[0] = Pred_TOF100[2];         legend.push_back(LegendNames[2]);
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicPrediction_TOF100");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  for(int Region=2; Region<DzRegions; Region++) {
    Histos[Region-2] = Pred_TOF110[Region];         legend.push_back(LegendNames[Region]);
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted", 0,0, 0,5);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicPrediction_TOF110");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  for(int Region=2; Region<DzRegions; Region++) {
    Histos[Region-2] = Pred_TOF120[Region];         legend.push_back(LegendNames[Region]);
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted", 0,0, 0,5);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicPrediction_TOF120");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  for(int Region=2; Region<DzRegions; Region++) {
    Histos[Region-2] = Pred_TOF130[Region];         legend.push_back(LegendNames[Region]);
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted", 0,0, 0,5);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Prediction_TOF130");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = StatSyst_TOF100;         legend.push_back("TOF > 1.0");
  Histos[1] = StatSyst_TOF110;         legend.push_back("TOF > 1.1");
  Histos[2] = StatSyst_TOF120;         legend.push_back("TOF > 1.2");
  //Histos[3] = StatSyst_TOF130;         legend.push_back("TOF > 1.3");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Stat+Syst Rel. Error", 0,0, 0,1.4);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicStatSyst");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Stat_TOF100;         legend.push_back("TOF > 1.0");
  Histos[1] = Stat_TOF110;         legend.push_back("TOF > 1.1");
  Histos[2] = Stat_TOF120;         legend.push_back("TOF > 1.2");
  //Histos[3] = Stat_TOF130;         legend.push_back("TOF > 1.3");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Stat Rel. Error", 0,0, 0,1.4);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicStat");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Syst_TOF100;         legend.push_back("TOF > 1.0");
  Histos[1] = Syst_TOF110;         legend.push_back("TOF > 1.1");
  Histos[2] = Syst_TOF120;         legend.push_back("TOF > 1.2");
  //Histos[3] = Syst_TOF130;         legend.push_back("TOF > 1.3");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Syst Rel. Error", 0,0, 0,1.4);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CosmicSyst");
  delete c1;

}
void CollisionBackgroundSystematic(string InputPattern){

  string SavePath  = InputPattern + "/Syst/";
  MakeDirectories(SavePath);
  TCanvas* c1;
  TH1** Histos = new TH1*[10];
  std::vector<string> legend;

  LegendTitle = "SA Only";

  TH1F* Pred_Def_TOF1p1 = new TH1F("Pred_Def_TOF1p1", "Pred_Def_TOF1p1", 14, 95, 515);
  TH1F* Pred_Low_TOF1p1 = new TH1F("Pred_Low_TOF1p1", "Pred_Low_TOF1p1", 14, 95, 515);
  TH1F* Pred_Rev_TOF1p1 = new TH1F("Pred_Rev_TOF1p1", "Pred_Rev_TOF1p1", 14, 95, 515);

  TH1F *StatSyst_TOF1p1 = new TH1F("StatSyst_TOF1p1", "StatSyst_TOF1p1", 14, 95, 515);
  TH1F *Syst_TOF1p1 = new TH1F("Syst_TOF1p1", "Syst_TOF1p1", 14, 95, 515);
  TH1F *Stat_TOF1p1 = new TH1F("Stat_TOF1p1", "Stat_TOF1p1", 14, 95, 515);

  TH1F* Pred_CSC_Def_TOF1p1 = new TH1F("Pred_CSC_Def_TOF1p1", "Pred_CSC_Def_TOF1p1", 14, 95, 515);
  TH1F* Pred_CSC_Low_TOF1p1 = new TH1F("Pred_CSC_Low_TOF1p1", "Pred_CSC_Low_TOF1p1", 14, 95, 515);
  TH1F* Pred_CSC_Rev_TOF1p1 = new TH1F("Pred_CSC_Rev_TOF1p1", "Pred_CSC_Rev_TOF1p1", 14, 95, 515);

  TH1F* Pred_DT_Def_TOF1p1 = new TH1F("Pred_DT_Def_TOF1p1", "Pred_DT_Def_TOF1p1", 14, 95, 515);
  TH1F* Pred_DT_Low_TOF1p1 = new TH1F("Pred_DT_Low_TOF1p1", "Pred_DT_Low_TOF1p1", 14, 95, 515);
  TH1F* Pred_DT_Rev_TOF1p1 = new TH1F("Pred_DT_Rev_TOF1p1", "Pred_DT_Rev_TOF1p1", 14, 95, 515);

  TH1F* Pred_Def_TOF1p15 = new TH1F("Pred_Def_TOF1p15", "Pred_Def_TOF1p15", 14, 95, 515);
  TH1F* Pred_Low_TOF1p15 = new TH1F("Pred_Low_TOF1p15", "Pred_Low_TOF1p15", 14, 95, 515);
  TH1F* Pred_Rev_TOF1p15 = new TH1F("Pred_Rev_TOF1p15", "Pred_Rev_TOF1p15", 14, 95, 515);

  TH1F *StatSyst_TOF1p15 = new TH1F("StatSyst_TOF1p15", "StatSyst_TOF1p15", 14, 95, 515);
  TH1F *Syst_TOF1p15 = new TH1F("Syst_TOF1p15", "Syst_TOF1p15", 14, 95, 515);
  TH1F *Stat_TOF1p15 = new TH1F("Stat_TOF1p15", "Stat_TOF1p15", 14, 95, 515);

  TH1F* Pred_Def_TOF1p2 = new TH1F("Pred_Def_TOF1p2", "Pred_Def_TOF1p2", 14, 95, 515);
  TH1F* Pred_Low_TOF1p2 = new TH1F("Pred_Low_TOF1p2", "Pred_Low_TOF1p2", 14, 95, 515);
  TH1F* Pred_Rev_TOF1p2 = new TH1F("Pred_Rev_TOF1p2", "Pred_Rev_TOF1p2", 14, 95, 515);

  TH1F *StatSyst_TOF1p2 = new TH1F("StatSyst_TOF1p2", "StatSyst_TOF1p2", 14, 95, 515);
  TH1F *Syst_TOF1p2 = new TH1F("Syst_TOF1p2", "Syst_TOF1p2", 14, 95, 515);
  TH1F *Stat_TOF1p2 = new TH1F("Stat_TOF1p2", "Stat_TOF1p2", 14, 95, 515);


  TFile* InputFile = new TFile((InputPattern + "Histos_Data.root").c_str());
  TFile* InputFileCosmic = new TFile((InputPattern + "Histos_Cosmic.root").c_str());

  TH1D*  HCuts_Pt       = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_Pt");
  TH1D*  HCuts_TOF      = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_TOF");

   TH1D*  H_A_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_Cen");
   TH1D*  H_B_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_Cen");
   TH1D*  H_C_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_Cen");

   TH1D*  H_A_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_For");
   TH1D*  H_B_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_For");
   TH1D*  H_C_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_For");

   TH1D*  H_A_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_Cen");
   TH1D*  H_B_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_Cen");
   TH1D*  H_C_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_Cen");

   TH1D*  H_A_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_For");
   TH1D*  H_B_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_For");
   TH1D*  H_C_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_For");

   TH1D*  H_A_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_Cen_Low");
   TH1D*  H_B_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_Cen_Low");
   TH1D*  H_C_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_Cen_Low");
   TH1D*  H_D_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_Cen_Low");

   TH1D*  H_A_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_For_Low");
   TH1D*  H_B_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_For_Low");
   TH1D*  H_C_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_For_Low");
   TH1D*  H_D_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_For_Low");

   TH1D*  H_A_Control_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_Cen_Low");
   TH1D*  H_B_Control_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_Cen_Low");
   TH1D*  H_C_Control_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_Cen_Low");
   TH1D*  H_D_Control_Cen_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_D_Cen_Low");

   TH1D*  H_A_Control_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_For_Low");
   TH1D*  H_B_Control_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_For_Low");
   TH1D*  H_C_Control_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_For_Low");
   TH1D*  H_D_Control_For_Low  = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_D_For_Low");

   TH1D*  H_Cosmic_FailDz      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/FailDz");
   TH1D*  H_Cosmic_Preselected      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/Preselected");

   double cosmicsPassCuts = H_Cosmic_Preselected->GetBinContent(1);
   double cosmicsFailCuts = H_Cosmic_FailDz->GetBinContent(1);

   double cosmicFraction=(cosmicsPassCuts/cosmicsFailCuts);
   double cosmicFractionErr=sqrt(pow(1./cosmicsFailCuts,2)*cosmicsPassCuts + pow(cosmicsPassCuts/(cosmicsFailCuts*cosmicsFailCuts),2)*cosmicsFailCuts);

  for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){

      double A_Cen=0, B_Cen=0, C_Cen=0;
      double A_For=0, B_For=0, C_For=0;

      A_Cen=H_A_Cen->GetBinContent(CutIndex+1);
      B_Cen=H_B_Cen->GetBinContent(CutIndex+1);
      C_Cen=H_C_Cen->GetBinContent(CutIndex+1);

      A_For=H_A_For->GetBinContent(CutIndex+1);
      B_For=H_B_For->GetBinContent(CutIndex+1);
      C_For=H_C_For->GetBinContent(CutIndex+1);

      double A_Cen_Control=0, B_Cen_Control=0, C_Cen_Control=0;
      double A_For_Control=0, B_For_Control=0, C_For_Control=0;

      A_Cen_Control=H_A_Control_Cen->GetBinContent(CutIndex+1);
      B_Cen_Control=H_B_Control_Cen->GetBinContent(CutIndex+1);
      C_Cen_Control=H_C_Control_Cen->GetBinContent(CutIndex+1);

      A_For_Control=H_A_Control_For->GetBinContent(CutIndex+1);
      B_For_Control=H_B_Control_For->GetBinContent(CutIndex+1);
      C_For_Control=H_C_Control_For->GetBinContent(CutIndex+1);

      double A_Cen_Err = A_Cen + cosmicFractionErr*cosmicFractionErr*A_Cen_Control*A_Cen_Control + cosmicFraction*cosmicFraction*A_Cen_Control;
      //double B_Cen_Err = B_Cen + cosmicFractionErr*cosmicFractionErr*B_Cen_Control*B_Cen_Control + cosmicFraction*cosmicFraction*B_Cen_Control;
      double C_Cen_Err = C_Cen + cosmicFractionErr*cosmicFractionErr*C_Cen_Control*C_Cen_Control + cosmicFraction*cosmicFraction*C_Cen_Control;

      A_Cen=A_Cen-cosmicFraction*A_Cen_Control;
      B_Cen=B_Cen-cosmicFraction*B_Cen_Control;
      C_Cen=C_Cen-cosmicFraction*C_Cen_Control;

      double A_For_Err = A_For + cosmicFractionErr*cosmicFractionErr*A_For_Control*A_For_Control + cosmicFraction*cosmicFraction*A_For_Control;
      //double B_For_Err = B_For + cosmicFractionErr*cosmicFractionErr*B_For_Control*B_For_Control + cosmicFraction*cosmicFraction*B_For_Control;
      double C_For_Err = C_For + cosmicFractionErr*cosmicFractionErr*C_For_Control*C_For_Control + cosmicFraction*cosmicFraction*C_For_Control;

      A_For=A_For-cosmicFraction*A_For_Control;
      B_For=B_For-cosmicFraction*B_For_Control;
      C_For=C_For-cosmicFraction*C_For_Control;

      double NPred_Cen=0;
      double Perr_Cen=0;

      double NPred_For=0;
      double Perr_For=0;

      if(A_Cen>0) {
        NPred_Cen    = ((C_Cen*B_Cen)/A_Cen);
        if(NPred_Cen<0) NPred_Cen=0;
        //if(A_Cen>0 && B_Cen>0 && C_Cen>0) Perr_Cen = sqrt( (pow(B_Cen/A_Cen,2)*C_Cen_Err) + (pow(C_Cen/A_Cen,2)*B_Cen_Err) + (pow((B_Cen*C_Cen/(A_Cen*A_Cen)),2)*A_Cen_Err) );
        if(A_Cen>0 && B_Cen>0 && C_Cen>0) Perr_Cen = sqrt( (pow(B_Cen/A_Cen,2)*C_Cen_Err) + (pow((B_Cen*C_Cen/(A_Cen*A_Cen)),2)*A_Cen_Err) );
        else Perr_Cen=0;
      }

      if(A_For>0) {
        NPred_For    = ((C_For*B_For)/A_For);
        if(NPred_For<0) NPred_For=0;
        //if(A_For>0 && B_For>0 && C_For>0) Perr_For = sqrt( (pow(B_For/A_For,2)*C_For_Err) + (pow(C_For/A_For,2)*B_For_Err) + (pow((B_For*C_For/(A_For*A_For)),2)*A_For_Err));
        if(A_For>0 && B_For>0 && C_For>0) Perr_For = sqrt( (pow(B_For/A_For,2)*C_For_Err) + (pow((B_For*C_For/(A_For*A_For)),2)*A_For_Err));
        else Perr_For=0;
      }

      double NPred_Coll = NPred_Cen+NPred_For;
      double Perr_Coll = sqrt(Perr_Cen*Perr_Cen + Perr_For*Perr_For);

      double A_Cen_Low=0, B_Cen_Low=0, C_Cen_Low=0, D_Cen_Low=0;
      double A_For_Low=0, B_For_Low=0, C_For_Low=0, D_For_Low=0;

      A_Cen_Low=H_A_Cen_Low->GetBinContent(CutIndex+1);
      B_Cen_Low=H_B_Cen_Low->GetBinContent(CutIndex+1);
      C_Cen_Low=H_C_Cen_Low->GetBinContent(CutIndex+1);
      D_Cen_Low=H_D_Cen_Low->GetBinContent(CutIndex+1);

      A_For_Low=H_A_For_Low->GetBinContent(CutIndex+1);
      B_For_Low=H_B_For_Low->GetBinContent(CutIndex+1);
      C_For_Low=H_C_For_Low->GetBinContent(CutIndex+1);
      D_For_Low=H_D_For_Low->GetBinContent(CutIndex+1);

      double A_Cen_Low_Control=0, B_Cen_Low_Control=0, C_Cen_Low_Control=0, D_Cen_Low_Control=0;
      double A_For_Low_Control=0, B_For_Low_Control=0, C_For_Low_Control=0, D_For_Low_Control=0;

      A_Cen_Low_Control=H_A_Control_Cen_Low->GetBinContent(CutIndex+1);
      B_Cen_Low_Control=H_B_Control_Cen_Low->GetBinContent(CutIndex+1);
      C_Cen_Low_Control=H_C_Control_Cen_Low->GetBinContent(CutIndex+1);
      D_Cen_Low_Control=H_D_Control_Cen_Low->GetBinContent(CutIndex+1);

      A_For_Low_Control=H_A_Control_For_Low->GetBinContent(CutIndex+1);
      B_For_Low_Control=H_B_Control_For_Low->GetBinContent(CutIndex+1);
      C_For_Low_Control=H_C_Control_For_Low->GetBinContent(CutIndex+1);
      D_For_Low_Control=H_D_Control_For_Low->GetBinContent(CutIndex+1);

      double A_Cen_Low_Err = A_Cen_Low + cosmicFractionErr*cosmicFractionErr*A_Cen_Low_Control*A_Cen_Low_Control + cosmicFraction*cosmicFraction*A_Cen_Low_Control;
      double B_Cen_Low_Err = B_Cen_Low + cosmicFractionErr*cosmicFractionErr*B_Cen_Low_Control*B_Cen_Low_Control + cosmicFraction*cosmicFraction*B_Cen_Low_Control;
      double C_Cen_Low_Err = C_Cen_Low + cosmicFractionErr*cosmicFractionErr*C_Cen_Low_Control*C_Cen_Low_Control + cosmicFraction*cosmicFraction*C_Cen_Low_Control;
      double D_Cen_Low_Err = D_Cen_Low + cosmicFractionErr*cosmicFractionErr*D_Cen_Low_Control*D_Cen_Low_Control + cosmicFraction*cosmicFraction*D_Cen_Low_Control;

      A_Cen_Low=A_Cen_Low-cosmicFraction*A_Cen_Low_Control;
      B_Cen_Low=B_Cen_Low-cosmicFraction*B_Cen_Low_Control;
      C_Cen_Low=C_Cen_Low-cosmicFraction*C_Cen_Low_Control;
      D_Cen_Low=D_Cen_Low-cosmicFraction*D_Cen_Low_Control;

      double A_For_Low_Err = A_For_Low + cosmicFractionErr*cosmicFractionErr*A_For_Low_Control*A_For_Low_Control + cosmicFraction*cosmicFraction*A_For_Low_Control;
      double B_For_Low_Err = B_For_Low + cosmicFractionErr*cosmicFractionErr*B_For_Low_Control*B_For_Low_Control + cosmicFraction*cosmicFraction*B_For_Low_Control;
      double C_For_Low_Err = C_For_Low + cosmicFractionErr*cosmicFractionErr*C_For_Low_Control*C_For_Low_Control + cosmicFraction*cosmicFraction*C_For_Low_Control;
      double D_For_Low_Err = D_For_Low + cosmicFractionErr*cosmicFractionErr*D_For_Low_Control*D_For_Low_Control + cosmicFraction*cosmicFraction*D_For_Low_Control;

      A_For_Low=A_For_Low-cosmicFraction*A_For_Low_Control;
      B_For_Low=B_For_Low-cosmicFraction*B_For_Low_Control;
      C_For_Low=C_For_Low-cosmicFraction*C_For_Low_Control;
      D_For_Low=D_For_Low-cosmicFraction*D_For_Low_Control;

      double NPred_Cen_Low=0;
      double Perr_Cen_Low=0;

      double NPred_For_Low=0;
      double Perr_For_Low=0;

      if(B_Cen_Low>0) {
        NPred_Cen_Low    = ((D_Cen_Low*B_Cen)/B_Cen_Low);
        if(NPred_Cen_Low<0) NPred_Cen_Low=0;
        //if(B_Cen_Low>0 && B_Cen>0 && D_Cen_Low>0) Perr_Cen_Low = sqrt( (pow(B_Cen/B_Cen_Low,2)*D_Cen_Low_Err) + (pow(D_Cen_Low/B_Cen_Low,2)*B_Cen_Err) + (pow((B_Cen*(D_Cen_Low)/(B_Cen_Low*B_Cen_Low)),2)*B_Cen_Low_Err) );
        if(B_Cen_Low>0 && B_Cen>0 && D_Cen_Low>0) Perr_Cen_Low = sqrt( (pow(B_Cen/B_Cen_Low,2)*D_Cen_Low_Err) + (pow((B_Cen*(D_Cen_Low)/(B_Cen_Low*B_Cen_Low)),2)*B_Cen_Low_Err));
        else Perr_Cen_Low=0;
      }

      if(B_For_Low>0) {
        NPred_For_Low    = ((D_For_Low*B_For)/B_For_Low);
        if(NPred_For_Low<0) NPred_For_Low=0;
        //if(B_For_Low>0 && B_For>0 && D_For_Low>0) Perr_For_Low = sqrt( (pow(B_For/B_For_Low,2)*D_For_Low_Err) + (pow(D_For_Low/B_For_Low,2)*B_For_Err) + (pow((B_For*(D_For_Low)/(B_For_Low*B_For_Low)),2)*B_For_Low_Err) );
        if(B_For_Low>0 && B_For>0 && D_For_Low>0) Perr_For_Low = sqrt( (pow(B_For/B_For_Low,2)*D_For_Low_Err) + (pow((B_For*(D_For_Low)/(B_For_Low*B_For_Low)),2)*B_For_Low_Err));
        else Perr_For_Low=0;
      }

      double NPred_Coll_Low = NPred_Cen_Low+NPred_For_Low;
      double Perr_Coll_Low = sqrt(Perr_Cen_Low*Perr_Cen_Low + Perr_For_Low*Perr_For_Low);

      double NPred_Cen_Rev=0;
      double Perr_Cen_Rev=0;

      double NPred_For_Rev=0;
      double Perr_For_Rev=0;

      if(A_Cen_Low>0) {
        NPred_Cen_Rev    = ((C_Cen_Low*B_Cen)/A_Cen_Low);
        if(NPred_Cen_Rev<0) NPred_Cen_Rev=0;
        //if(A_Cen_Low>0 && B_Cen>0 && C_Cen_Low>0) Perr_Cen_Rev = sqrt( (pow(B_Cen/A_Cen_Low,2)*C_Cen_Low_Err) + (pow(C_Cen_Low/A_Cen_Low,2)*B_Cen_Err) + (pow((B_Cen*(C_Cen_Low)/(A_Cen_Low*A_Cen_Low)),2)*A_Cen_Low_Err));
        if(A_Cen_Low>0 && B_Cen>0 && C_Cen_Low>0) Perr_Cen_Rev = sqrt( (pow(B_Cen/A_Cen_Low,2)*C_Cen_Low_Err) + (pow((B_Cen*(C_Cen_Low)/(A_Cen_Low*A_Cen_Low)),2)*A_Cen_Low_Err));
        else Perr_Cen_Rev=0;
      }

      if(A_For_Low>0) {
        NPred_For_Rev    = ((C_For_Low*B_For)/A_For_Low);
        if(NPred_For_Rev<0) NPred_For_Rev=0;
        //if(A_For_Low>0 && B_For>0 && C_For_Low>0) Perr_For_Rev = sqrt( (pow(B_For/A_For_Low,2)*C_For_Low_Err) + (pow(C_For_Low/A_For_Low,2)*B_For_Err) + (pow((B_For*(C_For_Low)/(A_For_Low*A_For_Low)),2)*A_For_Low_Err) );
        if(A_For_Low>0 && B_For>0 && C_For_Low>0) Perr_For_Rev = sqrt( (pow(B_For/A_For_Low,2)*C_For_Low_Err) + (pow((B_For*(C_For_Low)/(A_For_Low*A_For_Low)),2)*A_For_Low_Err));
        else Perr_For_Rev=0;
      }

      double NPred_Coll_Rev = NPred_Cen_Rev+NPred_For_Rev;
      double Perr_Coll_Rev = sqrt(Perr_Cen_Rev*Perr_Cen_Rev + Perr_For_Rev*Perr_For_Rev);

      int Bin=Pred_Def_TOF1p1->FindBin(HCuts_Pt->GetBinContent(CutIndex)/files);

      if(fabs(HCuts_TOF->GetBinContent(CutIndex)/files-1.1)<0.0001) {
	cout << "A Cen " << A_Cen << " for pt cut of " << HCuts_Pt->GetBinContent(CutIndex)/files << endl;
	Pred_Def_TOF1p1->SetBinContent(Bin, NPred_Coll);
	Pred_Def_TOF1p1->SetBinError(Bin, Perr_Coll);

	Pred_Low_TOF1p1->SetBinContent(Bin, NPred_Coll_Low);
        Pred_Low_TOF1p1->SetBinError(Bin, Perr_Coll_Low);

	Pred_Rev_TOF1p1->SetBinContent(Bin, NPred_Coll_Rev);
        Pred_Rev_TOF1p1->SetBinError(Bin, Perr_Coll_Rev);

        Pred_CSC_Def_TOF1p1->SetBinContent(Bin, NPred_Cen);
        Pred_CSC_Def_TOF1p1->SetBinError(Bin, Perr_Cen);

        Pred_CSC_Low_TOF1p1->SetBinContent(Bin, NPred_Cen_Low);
        Pred_CSC_Low_TOF1p1->SetBinError(Bin, Perr_Cen_Low);

        Pred_CSC_Rev_TOF1p1->SetBinContent(Bin, NPred_Cen_Rev);
        Pred_CSC_Rev_TOF1p1->SetBinError(Bin, Perr_Cen_Rev);

        Pred_DT_Def_TOF1p1->SetBinContent(Bin, NPred_For);
        Pred_DT_Def_TOF1p1->SetBinError(Bin, Perr_For);

        Pred_DT_Low_TOF1p1->SetBinContent(Bin, NPred_For_Low);
        Pred_DT_Low_TOF1p1->SetBinError(Bin, Perr_For_Low);

        Pred_DT_Rev_TOF1p1->SetBinContent(Bin, NPred_For_Rev);
        Pred_DT_Rev_TOF1p1->SetBinError(Bin, Perr_For_Rev);
      }

      if(fabs(HCuts_TOF->GetBinContent(CutIndex)/files-1.15)<0.0001) {
        Pred_Def_TOF1p15->SetBinContent(Bin, NPred_Coll);
	Pred_Def_TOF1p15->SetBinError(Bin, Perr_Coll);

        Pred_Low_TOF1p15->SetBinContent(Bin, NPred_Coll_Low);
        Pred_Low_TOF1p15->SetBinError(Bin, Perr_Coll_Low);

        Pred_Rev_TOF1p15->SetBinContent(Bin, NPred_Coll_Rev);
        Pred_Rev_TOF1p15->SetBinError(Bin, Perr_Coll_Rev);
      }

      if(fabs(HCuts_TOF->GetBinContent(CutIndex)/files-1.2)<0.0001) {
        Pred_Def_TOF1p2->SetBinContent(Bin, NPred_Coll);
	Pred_Def_TOF1p2->SetBinError(Bin, Perr_Coll);

        Pred_Low_TOF1p2->SetBinContent(Bin, NPred_Coll_Low);
        Pred_Low_TOF1p2->SetBinError(Bin, Perr_Coll_Low);

        Pred_Rev_TOF1p2->SetBinContent(Bin, NPred_Coll_Rev);
        Pred_Rev_TOF1p2->SetBinError(Bin, Perr_Coll_Rev);
      }
  }

  for(int i=0; i<Pred_Def_TOF1p1->GetNbinsX()+2; i++) {

  double StatSyst=0, Stat=0, Syst=0, Average=0;

  Average+=Pred_Def_TOF1p1->GetBinContent(i);
  Average+=Pred_Low_TOF1p1->GetBinContent(i);
  Average+=Pred_Rev_TOF1p1->GetBinContent(i);

  Stat+=Pred_Def_TOF1p1->GetBinError(i)*Pred_Def_TOF1p1->GetBinError(i);
  Stat+=Pred_Low_TOF1p1->GetBinError(i)*Pred_Low_TOF1p1->GetBinError(i);
  Stat+=Pred_Rev_TOF1p1->GetBinError(i)*Pred_Rev_TOF1p1->GetBinError(i);

  Average=Average/3.;

  StatSyst+=pow(Pred_Def_TOF1p1->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Low_TOF1p1->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Rev_TOF1p1->GetBinContent(i)-Average,2);  

  Stat=sqrt(Stat/3.);
  StatSyst=sqrt(StatSyst/2.);
  if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

  StatSyst_TOF1p1->SetBinContent(i, StatSyst/Average);
  Stat_TOF1p1->SetBinContent(i, Stat/Average);
  Syst_TOF1p1->SetBinContent(i, Syst/Average);

  StatSyst_TOF1p1->SetBinError(i, 0.0001);
  Stat_TOF1p1->SetBinError(i, 0.0001);
  Syst_TOF1p1->SetBinError(i, 0.0001);
  }

  for(int i=0; i<Pred_Def_TOF1p15->GetNbinsX()+2; i++) {

  double StatSyst=0, Stat=0, Syst=0, Average=0;

  Average+=Pred_Def_TOF1p15->GetBinContent(i);
  Average+=Pred_Low_TOF1p15->GetBinContent(i);
  Average+=Pred_Rev_TOF1p15->GetBinContent(i);

  Stat+=Pred_Def_TOF1p15->GetBinError(i)*Pred_Def_TOF1p15->GetBinError(i);
  Stat+=Pred_Low_TOF1p15->GetBinError(i)*Pred_Low_TOF1p15->GetBinError(i);
  Stat+=Pred_Rev_TOF1p15->GetBinError(i)*Pred_Rev_TOF1p15->GetBinError(i);

  Average=Average/3.;

  StatSyst+=pow(Pred_Def_TOF1p15->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Low_TOF1p15->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Rev_TOF1p15->GetBinContent(i)-Average,2);  

  Stat=sqrt(Stat/3.);
  StatSyst=sqrt(StatSyst/2.);
  if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

  StatSyst_TOF1p15->SetBinContent(i, StatSyst/Average);
  Stat_TOF1p15->SetBinContent(i, Stat/Average);
  Syst_TOF1p15->SetBinContent(i, Syst/Average);

  StatSyst_TOF1p15->SetBinError(i, 0.0001);
  Stat_TOF1p15->SetBinError(i, 0.0001);
  Syst_TOF1p15->SetBinError(i, 0.0001);
  }

  for(int i=0; i<Pred_Def_TOF1p2->GetNbinsX()+2; i++) {

  double StatSyst=0, Stat=0, Syst=0, Average=0;

  Average+=Pred_Def_TOF1p2->GetBinContent(i);
  Average+=Pred_Low_TOF1p2->GetBinContent(i);
  Average+=Pred_Rev_TOF1p2->GetBinContent(i);

  Stat+=Pred_Def_TOF1p2->GetBinError(i)*Pred_Def_TOF1p2->GetBinError(i);
  Stat+=Pred_Low_TOF1p2->GetBinError(i)*Pred_Low_TOF1p2->GetBinError(i);
  Stat+=Pred_Rev_TOF1p2->GetBinError(i)*Pred_Rev_TOF1p2->GetBinError(i);

  Average=Average/3.;

  StatSyst+=pow(Pred_Def_TOF1p2->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Low_TOF1p2->GetBinContent(i)-Average,2);
  StatSyst+=pow(Pred_Rev_TOF1p2->GetBinContent(i)-Average,2);  

  Stat=sqrt(Stat/3.);
  StatSyst=sqrt(StatSyst/2.);
  if(StatSyst*StatSyst>Stat*Stat) Syst=sqrt(StatSyst*StatSyst - Stat*Stat);

  StatSyst_TOF1p2->SetBinContent(i, StatSyst/Average);
  Stat_TOF1p2->SetBinContent(i, Stat/Average);
  Syst_TOF1p2->SetBinContent(i, Syst/Average);

  StatSyst_TOF1p2->SetBinError(i, 0.0001);
  Stat_TOF1p2->SetBinError(i, 0.0001);
  Syst_TOF1p2->SetBinError(i, 0.0001);
  }

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Pred_Def_TOF1p1;         legend.push_back("TOF > 1.0 && TOF < 1.1");
  Histos[1] = Pred_Rev_TOF1p1;         legend.push_back("TOF > 0.9 && TOF < 1.0");
  Histos[2] = Pred_Low_TOF1p1;         legend.push_back("TOF < 0.9");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted Tracks", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(true);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Coll_PredTracks1p1");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Pred_CSC_Def_TOF1p1;         legend.push_back("TOF > 1.0 && TOF < 1.1");
  Histos[1] = Pred_CSC_Rev_TOF1p1;         legend.push_back("TOF > 0.9 && TOF < 1.0");
  Histos[2] = Pred_CSC_Low_TOF1p1;         legend.push_back("TOF < 0.9");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted Tracks", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(true);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Coll_Cen_PredTracks1p1");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Pred_DT_Def_TOF1p1;         legend.push_back("TOF > 1.0 && TOF < 1.1");
  Histos[1] = Pred_DT_Rev_TOF1p1;         legend.push_back("TOF > 0.9 && TOF < 1.0");
  Histos[2] = Pred_DT_Low_TOF1p1;         legend.push_back("TOF < 0.9");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted Tracks", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(true);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Coll_For_PredTracks1p1");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Pred_Def_TOF1p15;         legend.push_back("TOF > 1.0 && TOF < 1.15");
  Histos[1] = Pred_Rev_TOF1p15;         legend.push_back("TOF > 0.9 && TOF < 1.0");
  Histos[2] = Pred_Low_TOF1p15;         legend.push_back("TOF < 0.9");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted Tracks", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(true);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Coll_PredTracks1p15");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Pred_Def_TOF1p2;         legend.push_back("TOF > 1.0 && TOF < 1.2");
  Histos[1] = Pred_Rev_TOF1p2;         legend.push_back("TOF > 0.9 && TOF < 1.0");
  Histos[2] = Pred_Low_TOF1p2;         legend.push_back("TOF < 0.9");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Predicted Tracks", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.8, 0.9, 0.4, 0.05);
  c1->SetLogy(true);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"Coll_PredTracks1p2");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Stat_TOF1p1;         legend.push_back("TOF>1.1");
  Histos[1] = Stat_TOF1p15;        legend.push_back("TOF>1.15");
  Histos[2] = Stat_TOF1p2;         legend.push_back("TOF>1.2");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Stat Rel. Error", 0,0, 0,0.06);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CollStat");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = StatSyst_TOF1p1;         legend.push_back("TOF>1.1");
  Histos[1] = StatSyst_TOF1p15;        legend.push_back("TOF>1.15");
  Histos[2] = StatSyst_TOF1p2;         legend.push_back("TOF>1.2");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Stat+Syst Rel. Error", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CollStatSyst");
  delete c1;

  c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
  Histos[0] = Syst_TOF1p1;         legend.push_back("TOF>1.1");
  Histos[1] = Syst_TOF1p15;        legend.push_back("TOF>1.15");
  Histos[2] = Syst_TOF1p2;         legend.push_back("TOF>1.2");
  DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Pt Cut", "Syst Rel. Error", 0,0, 0,0);
  DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
  c1->SetLogy(false);
  DrawPreliminary(IntegratedLuminosity);
  SaveCanvas(c1,SavePath,"CollSyst");
  delete c1;

}
