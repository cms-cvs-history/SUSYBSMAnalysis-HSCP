#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCutG.h" 
#include "TPaveText.h"
#include "tdrstyle.C"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_Samples.h"
//#include "CL95.h"
#include "roostats_cl95.C"
#include "nSigma.C"

using namespace std;



struct stAllInfo{
   double Mass;
   double XSec_Th;
   double XSec_Err;
   double XSec_Exp;
   double XSec_ExpUp;
   double XSec_ExpDown;
   double XSec_Exp2Up;
   double XSec_Exp2Down;
   double XSec_Obs;
   double Eff;
  double Eff_SystTOF;
  double Eff_SystPt;
   double Significance;
   double Index;
   double WP_Pt;
   double WP_TOF;
   float  NData;
  float  NData_Track;
  float  NData_NoTrack;
  float  NData_Cen;
  float  NData_For;
   float  NPred;
  float  NPred_Coll;
  float  NPred_Cosmic;
  float  NPred_Cen;
  float  NPred_For;
   float  NPredErr;
   float  NSign;

   stAllInfo(string path=""){
     Mass=-1; XSec_Th=-1; XSec_Err=-1; XSec_Exp=-1; XSec_ExpUp=-1;XSec_ExpDown=-1;XSec_Exp2Up=-1;XSec_Exp2Down=-1; XSec_Obs=-1; Eff=-1; 
      if(path=="")return;
      FILE* pFile = fopen(path.c_str(),"r");
      if(!pFile){printf("Can't open %s\n",path.c_str()); return;}
      fscanf(pFile,"Mass            : %lf\n",&Mass);
      fscanf(pFile,"Index           : %lf\n",&Index);
      fscanf(pFile,"WP_Pt           : %lf\n",&WP_Pt);
      fscanf(pFile,"WP_TOF          : %lf\n",&WP_TOF);
      fscanf(pFile,"Eff             : %lf\n",&Eff);
      fscanf(pFile,"Eff_SystTOF     : %lf\n",&Eff_SystTOF);
      fscanf(pFile,"Eff_SystPt      : %lf\n",&Eff_SystPt);
      fscanf(pFile,"Signif          : %lf\n",&Significance);
      fscanf(pFile,"XSec_Th         : %lf\n",&XSec_Th);
      fscanf(pFile,"XSec_Exp        : %lf\n",&XSec_Exp);
      fscanf(pFile,"XSec_ExpUp      : %lf\n",&XSec_ExpUp);
      fscanf(pFile,"XSec_ExpDown    : %lf\n",&XSec_ExpDown);
      fscanf(pFile,"XSec_Exp2Up     : %lf\n",&XSec_Exp2Up);
      fscanf(pFile,"XSec_Exp2Down   : %lf\n",&XSec_Exp2Down);
      fscanf(pFile,"XSec_Obs        : %lf\n",&XSec_Obs);
      fscanf(pFile,"NData           : %E\n" ,&NData);
      fscanf(pFile,"NData_Cen       : %E\n" ,&NData_Cen);
      fscanf(pFile,"NData_For       : %E\n" ,&NData_For);
      fscanf(pFile,"NData_Track     : %E\n" ,&NData_Track);
      fscanf(pFile,"NData_NoTrack   : %E\n" ,&NData_NoTrack);
      fscanf(pFile,"NPred           : %E\n" ,&NPred);
      fscanf(pFile,"NPred_Coll      : %E\n" ,&NPred_Coll);
      fscanf(pFile,"NPred_Cosmic    : %E\n" ,&NPred_Cosmic);
      fscanf(pFile,"NPred_Cen       : %E\n" ,&NPred_Cen);
      fscanf(pFile,"NPred_For       : %E\n" ,&NPred_For);
      fscanf(pFile,"NPredErr        : %E\n" ,&NPredErr);
      fscanf(pFile,"NSign           : %E\n" ,&NSign);
      fclose(pFile);
   }
};


struct stGraph{
   TGraph* GluinoF0;
   TGraph* GluinoF1;
   TGraph* GluinoF5;
   TGraph* GluinoTh;
   TCutG*  GluinoThErr;
};

double PlotMinScale = 0.0005;
double PlotMaxScale = 3;

TGraph* MakePlot(FILE* pFile, string InputPattern, string ModelName, int XSectionType=2, string Mass0="", string Mass1="", string Mass2="", string Mass3="", string Mass4="", string Mass5="", string Mass6="", string Mass7="", string Mass8="", string Mass9="",string Mass10="", string Mass11="", string Mass12="", string Mass13="");


void Exclusion(string pattern, string modelName, string signal, double Ratio_0C=-1, double Ratio_1C=-1, double Ratio_2C=-1);
void FindCutPoint(string pattern, string modelName, string signal, double Ratio_0C=-1, double Ratio_1C=-1, double Ratio_2C=-1);
int      JobIdToIndex(string JobId);

double FindIntersection(TGraph* obs, TGraph* th, double Min, double Max, double Step, double ThUncertainty=0, bool debug=false);
int ReadXSection(string InputFile, double* Mass, double* XSec, double* Low, double* High,  double* ErrLow, double* ErrHigh);
TCutG* GetErrorBand(string name, int N, double* Mass, double* Low, double* High, double MinLow=PlotMinScale, double MaxHigh=PlotMaxScale);
void CheckSignalUncertainty(FILE* pFile, string InputPattern);
void DrawModelLimitWithBand(string InputPattern, string inputmodel);
std::vector<string> GetModels(string inputmodel);
string GetModelName(string inputmodel);

char Buffer[2048];

int    CurrentSampleIndex;
string OutputPath;

std::vector<stSignal> signals;
//std::vector<double> signalsMeanHSCPPerEvent;

double RescaleFactor=1;
double RescaleError=0.1;
int Mode=0;
void Analysis_Step6(string MODE_="COMPILE", string InputPattern="", string modelName="", string signal="", double Ratio_0C=-1, double Ratio_1C=-1, double Ratio_2C=-1){
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.12);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
   gStyle->SetNdivisions(550,"Y");

   if(MODE_=="COMPILE")return;

   if(MODE_=="ANALYSE"){
      Exclusion(InputPattern, modelName, signal, Ratio_0C, Ratio_1C, Ratio_2C);
      return;
   }

   if(MODE_=="FINDCUT"){
     FindCutPoint(InputPattern, modelName, signal, Ratio_0C, Ratio_1C, Ratio_2C);
     return;
   }

   string Pattern  = "Results/Eta21/PtMin80/";

   string outpath = string("Results/EXCLUSION/");
   MakeDirectories(outpath);

   std::vector<string> ModelNames;
   ModelNames.push_back("All");

   std::vector<string> Models;
   Models.push_back("Gluinof1");
   Models.push_back("Gluinof5");
   Models.push_back("Gluinof10");

   TCanvas* c1;

   FILE* pFile = fopen((string("Analysis_Step6_Result.txt")).c_str(),"w");

   fprintf(pFile, "\\documentclass{article}\n");
   fprintf(pFile, "\\begin{document}\n\n");
   fprintf(pFile, "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");

   TGraph* Obs_GluinoF1  = MakePlot(pFile,Pattern,"Gluino (f=10\\%)", 2, "Gluino300_f1" , "Gluino400_f1" , "Gluino500_f1" , "Gluino600_f1" , "Gluino700_f1", "Gluino800_f1", "Gluino900_f1", "Gluino1000_f1", "Gluino1100_f1" , "Gluino1200_f1");
   TGraph* Obs_GluinoF5  = MakePlot(pFile,Pattern,"Gluino (f=50\\%)", 2, "Gluino300_f5" , "Gluino400_f5" , "Gluino500_f5" , "Gluino600_f5" , "Gluino700_f5", "Gluino800_f5", "Gluino900_f5", "Gluino1000_f5", "Gluino1100_f5" , "Gluino1200_f5");
   TGraph* Obs_GluinoF10 = MakePlot(pFile,Pattern,"Gluino (f=100\\%)",2, "Gluino300_f10" , "Gluino400_f10" , "Gluino500_f10" , "Gluino600_f10" , "Gluino700_f10", "Gluino800_f10", "Gluino900_f10", "Gluino1000_f10", "Gluino1100_f10" , "Gluino1200_f10");

   TGraph* Obs_Stop = MakePlot(pFile,Pattern,"Stop",2, "Stop300" , "Stop400" , "Stop500" , "Stop600" , "Stop700", "Stop800", "Stop900", "Stop1000", "Stop1100" , "Stop1200");

   fprintf(pFile,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(pFile, "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");

   fprintf(pFile,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(pFile,"\\end{document}\n\n");


   CheckSignalUncertainty(pFile,Pattern);

   double ThGluinoMass [100]; double ThGluinoXSec [100];  double ThGluinoLow  [100]; double ThGluinoHigh [100]; double ThGluinoErrLow  [100];  double ThGluinoErrHigh [100];
   int ThGluinoN = ReadXSection("gluino_XSec.txt", ThGluinoMass,ThGluinoXSec,ThGluinoLow,ThGluinoHigh, ThGluinoErrLow, ThGluinoErrHigh);
   TGraph* GluinoXSec    = new TGraph(ThGluinoN,ThGluinoMass,ThGluinoXSec);
   GluinoXSec->SetTitle("");
   GluinoXSec->GetYaxis()->SetTitleOffset(1.70);
   TCutG* GluinoXSecErr = GetErrorBand("gluinoErr",ThGluinoN,ThGluinoMass,ThGluinoLow,ThGluinoHigh);

   fprintf(pFile,"-----------------------\nNO TH UNCERTAINTY ACCOUNTED FOR   \n-------------------------\n");

   fprintf(pFile,"-----------------------\n0%% TK ONLY       \n-------------------------\n");
   fprintf(pFile,"MASS EXCLUDED UP TO %8.3fGeV for GluinoF1 \n", FindIntersection(Obs_GluinoF1,  GluinoXSec, 300, 1200, 1, 0.00));
   fprintf(pFile,"MASS EXCLUDED UP TO %8.3fGeV for GluinoF5 \n", FindIntersection(Obs_GluinoF5,  GluinoXSec, 300, 1200, 1, 0.00));
   fprintf(pFile,"MASS EXCLUDED UP TO %8.3fGeV for GluinoF10 \n", FindIntersection(Obs_GluinoF10,  GluinoXSec, 300, 1200, 1, 0.00));

   fclose(pFile);

   GluinoXSec      ->SetLineColor(4);  GluinoXSec      ->SetMarkerColor(4);   GluinoXSec      ->SetLineWidth(1);   GluinoXSec      ->SetLineStyle(3);  GluinoXSec      ->SetMarkerStyle(1);
   Obs_GluinoF1 ->SetLineColor(4);  Obs_GluinoF1 ->SetMarkerColor(4);   Obs_GluinoF1 ->SetLineWidth(2);   Obs_GluinoF1 ->SetLineStyle(1);  Obs_GluinoF1 ->SetMarkerStyle(22);
   Obs_GluinoF5 ->SetLineColor(4);  Obs_GluinoF5 ->SetMarkerColor(4);   Obs_GluinoF5 ->SetLineWidth(2);   Obs_GluinoF5 ->SetLineStyle(1);  Obs_GluinoF5 ->SetMarkerStyle(23);
   Obs_GluinoF10->SetLineColor(4);  Obs_GluinoF10->SetMarkerColor(4);   Obs_GluinoF10->SetLineWidth(2);   Obs_GluinoF10->SetLineStyle(1);  Obs_GluinoF10->SetMarkerStyle(23);

   TLegend* LEGTh = new TLegend(0.15,0.7,0.42,0.9);
   LEGTh->SetHeader("Theoretical Prediction");
   LEGTh->SetFillColor(0);
   LEGTh->SetBorderSize(0);
   TGraph* GlThLeg = (TGraph*) GluinoXSec->Clone("GluinoThLeg");
   GlThLeg->SetFillColor(GluinoXSecErr->GetFillColor());
   LEGTh->AddEntry(GlThLeg, "gluino (NLO+NLL)" ,"LF");

   c1 = new TCanvas("c1", "c1",600,600);
   TMultiGraph* MG = new TMultiGraph();
   MG->Add(GluinoXSec      ,"L");
   MG->Add(Obs_GluinoF1      ,"LP");
   MG->Add(Obs_GluinoF5      ,"LP");
   MG->Add(Obs_GluinoF10     ,"LP");
   MG->Draw("A");
   GluinoXSecErr->Draw("f");

   MG->Draw("same");
   MG->SetTitle("");
   MG->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
   MG->GetYaxis()->SetTitle("#sigma (pb)");
   MG->GetYaxis()->SetTitleOffset(1.70);
   MG->GetYaxis()->SetRangeUser(PlotMinScale,PlotMaxScale);
   
   DrawPreliminary(IntegratedLuminosity);
   
   TLegend* LEG = new TLegend(0.42,0.58,0.77,0.9);
   LEG->SetHeader("Tk - Only");
   LEG->SetFillColor(0); 
   LEG->SetBorderSize(0);
   LEG->AddEntry(Obs_GluinoF5 , "gluino;100% #tilde{g}g"    ,"LP");
   LEG->AddEntry(Obs_GluinoF5 , "gluino; 50% #tilde{g}g"    ,"LP");
   LEG->AddEntry(Obs_GluinoF1 , "gluino; 10% #tilde{g}g"    ,"LP");
   LEG->Draw();

   LEGTh->Draw();

//   c1->SetGridx(true);
//   c1->SetGridy(true);
   SaveCanvas(c1, outpath, string("Exclusion"));
   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("ExclusionLog"));
   delete c1;


   DrawModelLimitWithBand(Pattern, "Gluinof10");

   return; 
}


void CheckSignalUncertainty(FILE* pFile, string InputPattern){

  std::vector<string> Models;
  Models.push_back("Gluino300_f1");
  Models.push_back("Gluino400_f1");
  Models.push_back("Gluino500_f1");
  Models.push_back("Gluino600_f1");
  Models.push_back("Gluino700_f1");
  Models.push_back("Gluino800_f1");
  Models.push_back("Gluino900_f1");
  Models.push_back("Gluino1000_f1");
  Models.push_back("Gluino1100_f1");
  Models.push_back("Gluino1200_f1");
   Models.push_back("Gluino300_f5");
   Models.push_back("Gluino400_f5");
   Models.push_back("Gluino500_f5");
   Models.push_back("Gluino600_f5");
   Models.push_back("Gluino700_f5");
   Models.push_back("Gluino800_f5");
   Models.push_back("Gluino900_f5");
   Models.push_back("Gluino1000_f5");
   Models.push_back("Gluino1100_f5");
   Models.push_back("Gluino1200_f5");
   Models.push_back("Gluino300_f10");
   Models.push_back("Gluino400_f10");
   Models.push_back("Gluino500_f10");
   Models.push_back("Gluino600_f10");
   Models.push_back("Gluino700_f10");
   Models.push_back("Gluino800_f10");
   Models.push_back("Gluino900_f10");
   Models.push_back("Gluino1000_f10");
   Models.push_back("Gluino1100_f10");
   Models.push_back("Gluino1200_f10");

     fprintf(pFile, "%20s   Eff   --> PtScale |  TOFScale || TotalUncertainty\n","Model");

   for(unsigned int s=0;s<Models.size();s++){
        stAllInfo tmp(InputPattern+"/EXCLUSION" + "/"+Models[s]+".txt");
        double Pt = tmp.Eff - tmp.Eff_SystPt;
        double TOF = tmp.Eff - tmp.Eff_SystTOF;

	double Pttemp=max(Pt, 0.0);
        double TOFtemp=max(TOF, 0.0);

	
	fprintf(pFile, "%20s   %7.3f --> %7.3f  |  %7.3f || %7.3f\n",+Models[s].c_str(), tmp.Eff, Pt/tmp.Eff, TOF/tmp.Eff, sqrt(Pttemp*Pttemp + TOFtemp*TOFtemp)/tmp.Eff);

   }
}

TGraph* MakePlot(FILE* pFile, string InputPattern, string ModelName, int XSectionType, string Mass0, string Mass1, string Mass2, string Mass3, string Mass4, string Mass5, string Mass6, string Mass7, string Mass8, string Mass9,string Mass10, string Mass11, string Mass12, string Mass13){
   unsigned int N=0;
   stAllInfo Infos[14];

   if(Mass0!=""){Infos[0] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass0+".txt"); N=1;}
   if(Mass1!=""){Infos[1] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass1+".txt"); N=2;}
   if(Mass2!=""){Infos[2] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass2+".txt"); N=3;}
   if(Mass3!=""){Infos[3] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass3+".txt"); N=4;}
   if(Mass4!=""){Infos[4] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass4+".txt"); N=5;}
   if(Mass5!=""){Infos[5] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass5+".txt"); N=6;}
   if(Mass6!=""){Infos[6] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass6+".txt"); N=7;}
   if(Mass7!=""){Infos[7] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass7+".txt"); N=8;}
   if(Mass8!=""){Infos[8] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass8+".txt"); N=9;}
   if(Mass9!=""){Infos[9] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass9+".txt"); N=10;}
   if(Mass10!=""){Infos[10] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass10+".txt"); N=11;}
   if(Mass11!=""){Infos[11] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass11+".txt"); N=12;}
   if(Mass12!=""){Infos[12] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass12+".txt"); N=13;}
   if(Mass13!=""){Infos[13] = stAllInfo(InputPattern+"/EXCLUSION/"+Mass13+".txt"); N=14;}

   double Mass   [14];for(unsigned int i=0;i<14;i++){Mass   [i]=Infos[i].Mass;    }
   double XSecTh [14];for(unsigned int i=0;i<14;i++){XSecTh [i]=Infos[i].XSec_Th; }
   double XSecObs[14];for(unsigned int i=0;i<14;i++){XSecObs[i]=Infos[i].XSec_Obs;}
   double XSecExp[14];for(unsigned int i=0;i<14;i++){XSecExp[i]=Infos[i].XSec_Exp;}

   if(XSectionType>0)for(unsigned int i=0;i<N;i++)printf("%-18s %5.0f --> Pt>%+6.1f & TOF>%+4.3f --> NData=%2.0f  NPred=%6.1E+-%6.1E  NSign=%6.1E (Eff=%3.2f) Local Significance %3.2f\n",ModelName.c_str(),Infos[i].Mass,Infos[i].WP_Pt,Infos[i].WP_TOF, Infos[i].NData, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NSign, Infos[i].Eff, Infos[i].Significance);

   if(XSectionType>0){
   for(unsigned int i=0;i<N;i++){
     fprintf(pFile   ,"%-20s & %4.0f & %6.0f & %5.3f & %6.3f $\\pm$ %6.3f & %2.0f & %4.3f & %6.1E & %6.1E & %6.1E & %3.2f \\\\\n", ModelName.c_str(), Infos[i].Mass,  Infos[i].WP_Pt,Infos[i].WP_TOF, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NData, Infos[i].Eff, Infos[i].XSec_Th,Infos[i].XSec_Exp, Infos[i].XSec_Obs, Infos[i].Significance);
   }
   }
   
   TGraph* graph = NULL;
   if(XSectionType==0)graph = new TGraph(N,Mass,XSecTh);
   if(XSectionType==1)graph = new TGraph(N,Mass,XSecExp);
   if(XSectionType==2)graph = new TGraph(N,Mass,XSecObs);
   graph->SetTitle("");
   graph->GetYaxis()->SetTitle("CrossSection ( pb )");
   graph->GetYaxis()->SetTitleOffset(1.70);
   return graph;
}

void Exclusion(string pattern, string modelName, string signal, double Ratio_0C, double Ratio_1C, double Ratio_2C){
   GetSignalDefinition(signals);
   CurrentSampleIndex        = JobIdToIndex(signal); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return;} 

   cout << "Finding limit for " << signal << endl;

   stAllInfo toReturn(pattern+"/EXCLUSION" + "/Gluino1000_f10.txt");

   int CutIndex=toReturn.Index;
   toReturn.Mass      = signals[JobIdToIndex(signal)].Mass;
   toReturn.XSec_Th   = signals[JobIdToIndex(signal)].XSec;
   toReturn.XSec_Err  = signals[JobIdToIndex(signal)].XSec * 0.15;


   double RatioValue[] = {Ratio_0C, Ratio_1C, Ratio_2C};

   string outpath = pattern + "/EXCLUSION/";
   MakeDirectories(outpath);

   //FILE* pFile = fopen((outpath+"/"+modelName+".info").c_str(),"w");
   //if(!pFile)printf("Can't open file : %s\n",(outpath+"/"+modelName+".info").c_str());

   string InputPathSign     = pattern + "Histos.root";
   TFile* InputFileSign     = new TFile(InputPathSign.c_str());

   TH1D* TotalE          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/TotalE");

   TH1D* HSCPE[4];
   HSCPE[0]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/HSCPE");
   HSCPE[1]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC0/HSCPE");
   HSCPE[2]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC1/HSCPE");
   HSCPE[3]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC2/HSCPE");

   TH1D* HSCPE_SystTOF[4];
   HSCPE_SystTOF[0]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/HSCPE_SystTOF");
   HSCPE_SystTOF[1]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC0/HSCPE_SystTOF");
   HSCPE_SystTOF[2]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC1/HSCPE_SystTOF");
   HSCPE_SystTOF[3]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC2/HSCPE_SystTOF");

   TH1D* HSCPE_SystPt[4];
   HSCPE_SystPt[0]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/HSCPE_SystPt");
   HSCPE_SystPt[1]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC0/HSCPE_SystPt");
   HSCPE_SystPt[2]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC1/HSCPE_SystPt");
   HSCPE_SystPt[3]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC2/HSCPE_SystPt");

      double Eff = 0;
      double Eff_SystTOF = 0;
      double Eff_SystPt = 0;

      for(unsigned int i=0;i<3;i++){
	CurrentSampleIndex = JobIdToIndex(signal); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return; }
	double INTERN_Eff = HSCPE[i+1]->GetBinContent(CutIndex+1)/TotalE->Integral();
	Eff+= INTERN_Eff   * RatioValue[i];

        double INTERN_Eff_SystTOF = HSCPE_SystTOF[i+1]->GetBinContent(CutIndex+1)/TotalE->Integral();
        Eff_SystTOF+= INTERN_Eff_SystTOF   * RatioValue[i];

        double INTERN_Eff_SystPt = HSCPE_SystPt[i+1]->GetBinContent(CutIndex+1)/TotalE->Integral();
        Eff_SystPt+= INTERN_Eff_SystPt   * RatioValue[i];
      }

     toReturn.Eff = Eff;
     toReturn.Eff_SystTOF = Eff_SystTOF;
     toReturn.Eff_SystPt = Eff_SystPt;
     toReturn.NSign = Eff*(signals[CurrentSampleIndex].XSec*IntegratedLuminosity);

     //fclose(pFile);   


   LimitResult CLMResults;
   //double signalUncertainty=0.10;

   //CLMResults =  roostats_limit(IntegratedLuminosity, IntegratedLuminosity*0.045, Eff, Eff*signalUncertainty,toReturn.NPred, toReturn.NPredErr, toReturn.NData, false, 1, "cls", "", 12345);

   toReturn.XSec_Exp  = CLMResults.GetExpectedLimit();
   toReturn.XSec_ExpUp    = CLMResults.GetOneSigmaHighRange();
   toReturn.XSec_ExpDown  = CLMResults.GetOneSigmaLowRange();
   toReturn.XSec_Exp2Up   = CLMResults.GetTwoSigmaHighRange();
   toReturn.XSec_Exp2Down = CLMResults.GetTwoSigmaLowRange();
   toReturn.XSec_Obs  = CLMResults.GetObservedLimit();
   toReturn.Significance = nSigma(toReturn.NPred, toReturn.NData, toReturn.NPredErr/toReturn.NPred);

     FILE* pFile2 = fopen((outpath+"/"+modelName+".txt").c_str(),"w");
     if(!pFile2)printf("Can't open file : %s\n",(outpath+"/"+modelName+".txt").c_str());
     fprintf(pFile2,"Mass         : %f\n",signals[JobIdToIndex(signal)].Mass);
     fprintf(pFile2,"Index        : %f\n",toReturn.Index);
     fprintf(pFile2,"WP_Pt        : %f\n",toReturn.WP_Pt);
     fprintf(pFile2,"WP_TOF       : %f\n",toReturn.WP_TOF);
     fprintf(pFile2,"Eff          : %f\n",toReturn.Eff);
     fprintf(pFile2,"Eff_SystTOF  : %f\n",toReturn.Eff_SystTOF);
     fprintf(pFile2,"Eff_SystPt   : %f\n",toReturn.Eff_SystPt);
     fprintf(pFile2,"Signif       : %f\n",toReturn.Significance);
     fprintf(pFile2,"XSec_Th      : %f\n",toReturn.XSec_Th);
     fprintf(pFile2,"XSec_Exp     : %f\n",toReturn.XSec_Exp);
     fprintf(pFile2,"XSec_ExpUp   : %f\n",toReturn.XSec_ExpUp);
     fprintf(pFile2,"XSec_ExpDown : %f\n",toReturn.XSec_ExpDown);
     fprintf(pFile2,"XSec_Exp2Up  : %f\n",toReturn.XSec_Exp2Up);
     fprintf(pFile2,"XSec_Exp2Down: %f\n",toReturn.XSec_Exp2Down);
     fprintf(pFile2,"XSec_Obs     : %f\n",toReturn.XSec_Obs);     
     fprintf(pFile2,"NData        : %+6.2E\n",toReturn.NData);
     fprintf(pFile2,"NData_Cen    : %+6.2E\n",toReturn.NData_Cen);
     fprintf(pFile2,"NData_For    : %+6.2E\n",toReturn.NData_For);
     fprintf(pFile2,"NData_Track  : %+6.2E\n",toReturn.NData_Track);
     fprintf(pFile2,"NData_NoTrack: %+6.2E\n",toReturn.NData_NoTrack);
     fprintf(pFile2,"NPred        : %+6.2E\n",toReturn.NPred);
     fprintf(pFile2,"NPred_Coll   : %+6.2E\n",toReturn.NPred_Coll);
     fprintf(pFile2,"NPred_Cosmic : %+6.2E\n",toReturn.NPred_Cosmic);
     fprintf(pFile2,"NPred_Cen    : %+6.2E\n",toReturn.NPred_Cen);
     fprintf(pFile2,"NPred_For    : %+6.2E\n",toReturn.NPred_For);
     fprintf(pFile2,"NPredErr     : %+6.2E\n",toReturn.NPredErr);
     fprintf(pFile2,"NSign        : %+6.2E\n",toReturn.NSign);

     fclose(pFile2);
     return;
}



void FindCutPoint(string pattern, string modelName, string signal, double Ratio_0C, double Ratio_1C, double Ratio_2C){

   GetSignalDefinition(signals);
   CurrentSampleIndex        = JobIdToIndex(signal); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return;  } 

   stAllInfo toReturn;
   toReturn.Mass      = signals[JobIdToIndex(signal)].Mass;
   toReturn.Index     = 0;
   toReturn.WP_Pt     = 0;
   toReturn.WP_TOF    = 0;
   toReturn.XSec_Th   = signals[JobIdToIndex(signal)].XSec;
   toReturn.XSec_Err  = signals[JobIdToIndex(signal)].XSec * 0.15;
   toReturn.XSec_Exp  = 1E50;
   toReturn.XSec_ExpUp    = 1E50;
   toReturn.XSec_ExpDown  = 1E50;
   toReturn.XSec_Exp2Up   = 1E50;
   toReturn.XSec_Exp2Down = 1E50;
   toReturn.XSec_Obs  = 1E50;
   toReturn.Eff       = 0;
   toReturn.NData     = 0;
   toReturn.NData_Track  = 0;
   toReturn.NData_NoTrack= 0;
   toReturn.NData_Cen     = 0;
   toReturn.NData_For     = 0;
   toReturn.NPred     = 1E50;
   toReturn.NPred_Coll  = 1E50;
   toReturn.NPred_Cosmic= 1E50;
   toReturn.NPred_Cen  = 1E50;
   toReturn.NPred_For  = 1E50;
   toReturn.NPredErr  = 1E50;
   toReturn.NSign     = 1E50;



   double RatioValue[] = {Ratio_0C, Ratio_1C, Ratio_2C};

   double MaxSOverB=-1; 
   int MaxSOverBIndex=-1;

   string outpath = pattern + "/EXCLUSION/";
   MakeDirectories(outpath);

   FILE* pFile = fopen((outpath+"/"+modelName+".info").c_str(),"w");
   if(!pFile)printf("Can't open file : %s\n",(outpath+"/"+modelName+".info").c_str());

   string InputPath = pattern + "Histos_Data.root";

   TFile* InputFile     = new TFile(InputPath.c_str());

   TH1D*  H_B     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B");
   TH1D*  H_C     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C");
   TH1D*  H_D     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D");

   TH1D*  H_A_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_Cen");
   TH1D*  H_B_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_Cen");
   TH1D*  H_C_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_Cen");
   TH1D*  H_D_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_Cen");

   TH1D*  H_A_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_For");
   TH1D*  H_B_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_For");
   TH1D*  H_C_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_For");
   TH1D*  H_D_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_For");

   TH1D*  H_B_Control     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B");

   TH1D*  H_A_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_Cen");
   TH1D*  H_B_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_Cen");
   TH1D*  H_C_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_Cen");

   TH1D*  H_A_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_A_For");
   TH1D*  H_B_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_B_For");
   TH1D*  H_C_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_C_For");

   TH1D*  H_D_Control_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_D_Cen");
   TH1D*  H_D_Control_For     = (TH1D*)GetObjectFromPath(InputFile, "Data_Control/H_D_For");
   TH1D*  H_D_Track     = (TH1D*)GetObjectFromPath(InputFile, "Data_Track/H_D");
   TH1D*  H_D_NoTrack     = (TH1D*)GetObjectFromPath(InputFile, "Data_NoTrack/H_D");

   string InputPathCosmic     = pattern + "Histos_Cosmic.root";
   TFile* InputFileCosmic     = new TFile(InputPathCosmic.c_str());

   TH1D*  HCuts_Pt      = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_Pt");
   TH1D*  HCuts_TOF     = (TH1D*)GetObjectFromPath(InputFileCosmic, "HCuts_TOF");

   TH1D*  H_Cosmic_FailDz_Cen      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/FailDz_DT");
   TH1D*  H_Cosmic_Preselected_Cen      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/Preselected_DT");
   TH1D*  H_Cosmic_FailDz_For      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/FailDz_CSC");
   TH1D*  H_Cosmic_Preselected_For      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/Preselected_CSC");

   TH1D*  H_B_Cen_Cosmic     = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/H_B_Cen");
   TH1D*  H_D_Cen_Cosmic      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/H_D_Cen" );
   TH1D*  H_B_For_Cosmic      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/H_B_For");
   TH1D*  H_D_For_Cosmic      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/H_D_For");
   TH1D*  H_D_Cosmic      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic/H_D" );

   TH1D*  H_B_Cosmic_Control     = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic_Control/H_B");
   TH1D*  H_D_Cosmic_Control     = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic_Control/H_D");
   TH1D*  H_B_Cen_Cosmic_Control     = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic_Control/H_B_Cen");
   TH1D*  H_B_For_Cosmic_Control      = (TH1D*)GetObjectFromPath(InputFileCosmic, "Cosmic_Control/H_B_For");


   string InputPathSign     = pattern + "Histos.root";
   TFile* InputFileSign     = new TFile(InputPathSign.c_str());

   TH1D* TotalE          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/TotalE");
   TH1D* HSCPE[4];
   HSCPE[0]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "/HSCPE");
   HSCPE[1]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC0/HSCPE");
   HSCPE[2]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC1/HSCPE");
   HSCPE[3]          = (TH1D*)GetObjectFromPath(InputFileSign, signals[CurrentSampleIndex].Name + "_NC2/HSCPE");

   ///##############################################################################"

   double cosmicsPassCuts_Cen  = H_Cosmic_Preselected_Cen->GetBinContent(1);
   double cosmicsFailCuts_Cen = H_Cosmic_FailDz_Cen->GetBinContent(1);

   double cosmicsPassCuts_For  = H_Cosmic_Preselected_For->GetBinContent(1);
   double cosmicsFailCuts_For = H_Cosmic_FailDz_For->GetBinContent(1);

   double cosmicFraction_Cen=(cosmicsPassCuts_Cen/cosmicsFailCuts_Cen);
   double cosmicFractionErr_Cen=sqrt(pow(1./cosmicsFailCuts_Cen,2)*cosmicsPassCuts_Cen +pow(cosmicsPassCuts_Cen/(cosmicsFailCuts_Cen*cosmicsFailCuts_Cen),2)*cosmicsFailCuts_Cen);
   double cosmicFraction_For=(cosmicsPassCuts_For/cosmicsFailCuts_For);
   double cosmicFractionErr_For=sqrt(pow(1./cosmicsFailCuts_For,2)*cosmicsPassCuts_For +pow(cosmicsPassCuts_For/(cosmicsFailCuts_For*cosmicsFailCuts_For),2)*cosmicsFailCuts_For);

   //Going to first loop and find the cut with the min S over sqrt(B) because this is quick and normally gives a cut with a reach near the minimum
   stAllInfo CutInfo[HCuts_Pt->GetNbinsX()];
   for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++) CutInfo[CutIndex]=toReturn;
   for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){

      if(H_C->GetBinContent(CutIndex+1)<15 || H_B->GetBinContent(CutIndex+1)<15)continue;

      double B_Cen_Cosmic=H_B_Cen_Cosmic->GetBinContent(CutIndex+1);
      double D_Cen_Cosmic=H_D_Cen_Cosmic->GetBinContent(CutIndex+1);
      double B_For_Cosmic=H_B_For_Cosmic->GetBinContent(CutIndex+1);
      double D_For_Cosmic=H_D_For_Cosmic->GetBinContent(CutIndex+1);
      double D_Cosmic=H_D_Cosmic->GetBinContent(CutIndex+1);

      double B_Cen_Cosmic_Control=H_B_Cen_Cosmic_Control->GetBinContent(CutIndex+1);
      double B_For_Cosmic_Control=H_B_For_Cosmic_Control->GetBinContent(CutIndex+1);
      double B_Cosmic_Control=H_B_Cosmic_Control->GetBinContent(CutIndex+1);
      double D_Cosmic_Control=H_D_Cosmic_Control->GetBinContent(CutIndex+1);

      double NData_Track = H_D_Track->GetBinContent(CutIndex+1);
      double NData_NoTrack = H_D_NoTrack->GetBinContent(CutIndex+1);
      double NData = H_D->GetBinContent(CutIndex+1);
      double NData_Cen = H_D_Cen->GetBinContent(CutIndex+1);
      double NData_For = H_D_For->GetBinContent(CutIndex+1);

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

      double B_Control=H_B_Control->GetBinContent(CutIndex+1);

      double D_Control_Cen=H_D_Control_Cen->GetBinContent(CutIndex+1);
      double D_Control_For=H_D_Control_For->GetBinContent(CutIndex+1);

      double A_Cen_Err = A_Cen + cosmicFractionErr_Cen*cosmicFractionErr_Cen*A_Cen_Control*A_Cen_Control + cosmicFraction_Cen*cosmicFraction_Cen*A_Cen_Control;
      double B_Cen_Err = B_Cen + cosmicFractionErr_Cen*cosmicFractionErr_Cen*B_Cen_Control*B_Cen_Control + cosmicFraction_Cen*cosmicFraction_Cen*B_Cen_Control;
      double C_Cen_Err = C_Cen + cosmicFractionErr_Cen*cosmicFractionErr_Cen*C_Cen_Control*C_Cen_Control + cosmicFraction_Cen*cosmicFraction_Cen*C_Cen_Control;

      A_Cen=A_Cen-cosmicFraction_Cen*A_Cen_Control;
      B_Cen=B_Cen-cosmicFraction_Cen*B_Cen_Control;
      C_Cen=C_Cen-cosmicFraction_Cen*C_Cen_Control;

      double A_For_Err = A_For + cosmicFractionErr_For*cosmicFractionErr_For*A_For_Control*A_For_Control + cosmicFraction_For*cosmicFraction_For*A_For_Control;
      double B_For_Err = B_For + cosmicFractionErr_For*cosmicFractionErr_For*B_For_Control*B_For_Control + cosmicFraction_For*cosmicFraction_For*B_For_Control;
      double C_For_Err = C_For + cosmicFractionErr_For*cosmicFractionErr_For*C_For_Control*C_For_Control + cosmicFraction_For*cosmicFraction_For*C_For_Control;

      A_For=A_For-cosmicFraction_For*A_For_Control;
      B_For=B_For-cosmicFraction_For*B_For_Control;
      C_For=C_For-cosmicFraction_For*C_For_Control;

      double NPred_Cen=0;
      double Perr_Cen=0;

      double NPred_For=0;
      double Perr_For=0;

      if(A_Cen>0) {
        NPred_Cen    = ((C_Cen*B_Cen)/A_Cen);
        if(NPred_Cen<0) NPred_Cen=0;
        if(A_Cen>0 && B_Cen>0 && C_Cen>0) Perr_Cen = sqrt( (pow(B_Cen/A_Cen,2)*C_Cen_Err) + (pow(C_Cen/A_Cen,2)*B_Cen_Err) + (pow((B_Cen*(C_Cen)/(A_Cen*A_Cen)),2)*A_Cen_Err) );
        else Perr_Cen=0;
      }

      if(A_For>0) {
        NPred_For    = ((C_For*B_For)/A_For);
        if(NPred_For<0) NPred_For=0;
        if(A_For>0 && B_For>0 && C_For>0) Perr_For = sqrt( (pow(B_For/A_For,2)*C_For_Err) + (pow(C_For/A_For,2)*B_For_Err) + (pow((B_For*(C_For)/(A_For*A_For)),2)*A_For_Err) );
        else Perr_For=0;
      }

      double NPred_Coll = NPred_Cen+NPred_For;

      double NPred_Cosmic = 0;
      //if(B_Cen_Cosmic>0 && B_For_Cosmic>0) NPred_Cosmic = B_Cen_Control*cosmicFraction_Cen*D_Cen_Cosmic/B_Cen_Cosmic + B_For_Control*cosmicFraction_For*D_For_Cosmic/B_For_Cosmic;
      if((B_Cosmic_Control+D_Cosmic_Control)>0) NPred_Cosmic = B_Control*D_Cosmic/(B_Cosmic_Control+D_Cosmic_Control);
      double NPred = NPred_Coll+NPred_Cosmic;
      cout << "One cosmic in control forward region would be " << D_Cosmic/(B_Cosmic_Control+D_Cosmic_Control) << endl;

      double Coll_Syst=0.08;
      double Perr_Coll = sqrt(Perr_Cen*Perr_Cen + Perr_For*Perr_For + NPred_Coll*NPred_Coll*Coll_Syst*Coll_Syst);
      double Cosmic_Syst=0.6;

      double Perr_Cosmic = D_Cosmic/(B_Cosmic_Control+D_Cosmic_Control);
      //double Perr_Cosmic = sqrt(cosmicFraction_Cen*cosmicFraction_Cen*D_Control_Cen + cosmicFractionErr_Cen*cosmicFractionErr_Cen*D_Control_Cen*D_Control_Cen + 
      //				cosmicFraction_For*cosmicFraction_For*D_Control_For + cosmicFractionErr_For*cosmicFractionErr_For*D_Control_For*D_Control_For +
      //				D_Control_Cen*cosmicFraction_Cen*Cosmic_Syst*D_Control_Cen*cosmicFraction_Cen*Cosmic_Syst);

	double NPredErr=sqrt(Perr_Coll*Perr_Coll+Perr_Cosmic*Perr_Cosmic);


      if(NPred<=0){continue;} //Is <=0 only when prediction failed or is not meaningful (i.e. WP=(0,0,0) )

      double Eff       = 0;

      bool Mix012C = (RatioValue[0]>0 || RatioValue[1]>0 || RatioValue[2]>0);

      if(Mix012C) {
	for(unsigned int i=0;i<3;i++){
	  CurrentSampleIndex        = JobIdToIndex(signal); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return; }
	  double INTERN_Eff       = HSCPE[i+1]->GetBinContent(CutIndex+1)/TotalE->Integral();
	  Eff                      += INTERN_Eff   * RatioValue[i];
	}
      }
      else {
	Eff=HSCPE[0]->GetBinContent(CutIndex+1)/TotalE->Integral();
      }

      if(Eff==0)continue;

      fprintf(pFile  ,"%10s: Testing CutIndex=%4i (Pt>%6.2f TOF>%6.3f) Ndata=%+6.2E Ndata Cen=%+6.2E Ndata For=%+6.2E Ndata No Track=%+6.2E NPred=%6.3E+-%6.3E NPred Coll=%6.3E NPred Cosmic=%6.3E NPred Cen=%6.3E NPred For=%6.3E SignalEff=%6.3f\n\n",signal.c_str(),CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1)/files, HCuts_TOF->GetBinContent(CutIndex+1)/files,NData,NData_Cen, NData_For,NData_NoTrack, NPred, NPredErr,NPred_Coll, NPred_Cosmic, NPred_Cen, NPred_For, Eff);fflush(stdout);
      fprintf(stdout ,"%10s: Testing CutIndex=%4i (Pt>%6.2f TOF>%6.3f) Ndata=%+6.2E Ndata Cen=%+6.2E Ndata For=%+6.2E Ndata No Track=%+6.2E NPred=%6.3E+-%6.3E NPred Coll=%6.3E NPred Cosmic=%6.3E NPred Cen=%6.3E NPred For=%6.3E SignalEff=%6.3f\n\n",signal.c_str(),CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1)/files, HCuts_TOF->GetBinContent(CutIndex+1)/files,NData,NData_Cen, NData_For,NData_NoTrack,NPred, NPredErr,NPred_Coll, NPred_Cosmic, NPred_Cen, NPred_For, Eff);fflush(stdout);
      //Need to divide cut values by files because of the merging step combininng files


      if(Eff/sqrt(max(0.1, NPred))>MaxSOverB) {MaxSOverB=Eff/sqrt(max(0.1, NPred)); MaxSOverBIndex=CutIndex;}

     toReturn.Index     = CutIndex;
     toReturn.WP_Pt     = HCuts_Pt ->GetBinContent(CutIndex+1)/files;
     toReturn.WP_TOF    = HCuts_TOF->GetBinContent(CutIndex+1)/files;
     toReturn.XSec_Th   = signals[JobIdToIndex(signal)].XSec;
     toReturn.XSec_Err  = signals[JobIdToIndex(signal)].XSec * 0.15;
     toReturn.Eff       = Eff;
     toReturn.NData     = NData;
     toReturn.NData_Track = NData_Track;
     toReturn.NData_NoTrack = NData_NoTrack;
     toReturn.NData_Cen     = NData_Cen;
     toReturn.NData_For     = NData_For;
     toReturn.NPred     = NPred;
     toReturn.NPred_Coll = NPred_Coll;
     toReturn.NPred_Cosmic = NPred_Cosmic;
     toReturn.NPred_Cen = NPred_Cen;
     toReturn.NPred_For = NPred_For;
     toReturn.NPredErr  = NPredErr;
     toReturn.NSign     = Eff*(signals[CurrentSampleIndex].XSec*IntegratedLuminosity);

     CutInfo[CutIndex]=toReturn;
   }


   fclose(pFile);   
   //Find reach for point with best S Over sqrt(B) first.
   double NPredSB=CutInfo[MaxSOverBIndex].NPred;
   double NPredErrSB=CutInfo[MaxSOverBIndex].NPredErr;
   double EffSB=CutInfo[MaxSOverBIndex].Eff;
   double FiveSigma=1E50;
   for (int n_obs=5; n_obs<1000; n_obs++) {
     if(nSigma(NPredSB, n_obs, NPredErrSB/NPredSB)>=5) {
       FiveSigma=n_obs;
       break;
     }
   }

   double MinReach=(FiveSigma-NPredSB)/(EffSB*IntegratedLuminosity);
   cout << endl << "Min Reach " << MinReach << endl;
   toReturn=CutInfo[MaxSOverBIndex]; // In case this point does give the best reach avoids rounding errors


   for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){
     double NPred=CutInfo[CutIndex].NPred;
     double NPredErr=CutInfo[CutIndex].NPredErr;
     double Eff=CutInfo[CutIndex].Eff;
     if(Eff==0) continue;  //Eliminate points where prediction could not be made
     FiveSigma=1E50;
     for (int n_obs=5; n_obs<1000; n_obs++) {
       if(n_obs<(NPred+3*sqrt(NPred))) continue;    //5 sigma implies more than 5 times sqrt(B) excess so can cut these points, put it at 3 to be safe
       double thisReach=(n_obs-NPred)/(Eff*IntegratedLuminosity);
       if(thisReach>=MinReach) break;    // This selection point will not give the optimum reach so move on
       if(nSigma(NPred, n_obs, NPredErr/NPred)>=5) {
	 FiveSigma=n_obs;
	 break;
       }
     }
     double Reach=(FiveSigma-NPred)/(Eff*IntegratedLuminosity);

     if(Reach>MinReach) continue;
     MinReach=Reach;
     toReturn=CutInfo[CutIndex];
   }

   toReturn.Significance = nSigma(toReturn.NPred, toReturn.NData, toReturn.NPredErr/toReturn.NPred);

     FILE* pFile2 = fopen((outpath+"/"+modelName+".txt").c_str(),"w");
     if(!pFile2)printf("Can't open file : %s\n",(outpath+"/"+modelName+".txt").c_str());
     fprintf(pFile2,"Mass            : %f\n",signals[JobIdToIndex(signal)].Mass);
     fprintf(pFile2,"Index           : %f\n",toReturn.Index);
     fprintf(pFile2,"WP_Pt           : %f\n",toReturn.WP_Pt);
     fprintf(pFile2,"WP_TOF          : %f\n",toReturn.WP_TOF);
     fprintf(pFile2,"Eff             : %f\n",toReturn.Eff);
     fprintf(pFile2,"Signif          : %f\n",toReturn.Significance);
     fprintf(pFile2,"XSec_Th         : %f\n",toReturn.XSec_Th);
     fprintf(pFile2,"XSec_Exp        : %f\n",toReturn.XSec_Exp);
     fprintf(pFile2,"XSec_ExpUp      : %f\n",toReturn.XSec_ExpUp);
     fprintf(pFile2,"XSec_ExpDown    : %f\n",toReturn.XSec_ExpDown);
     fprintf(pFile2,"XSec_Exp2Up     : %f\n",toReturn.XSec_Exp2Up);
     fprintf(pFile2,"XSec_Exp2Down   : %f\n",toReturn.XSec_Exp2Down);
     fprintf(pFile2,"XSec_Obs        : %f\n",toReturn.XSec_Obs);     
     fprintf(pFile2,"NData           : %+6.2E\n",toReturn.NData);
     fprintf(pFile2,"NData_Cen       : %+6.2E\n",toReturn.NData_Cen);
     fprintf(pFile2,"NData_For       : %+6.2E\n",toReturn.NData_For);
     fprintf(pFile2,"NData_Track     : %+6.2E\n",toReturn.NData_Track);
     fprintf(pFile2,"NData_NoTrack   : %+6.2E\n",toReturn.NData_NoTrack);
     fprintf(pFile2,"NPred           : %+6.2E\n",toReturn.NPred);
     fprintf(pFile2,"NPred_Coll      : %+6.2E\n",toReturn.NPred_Coll);
     fprintf(pFile2,"NPred_Cosmic    : %+6.2E\n",toReturn.NPred_Cosmic);
     fprintf(pFile2,"NPred_Cen       : %+6.2E\n",toReturn.NPred_Cen);
     fprintf(pFile2,"NPred_For       : %+6.2E\n",toReturn.NPred_For);
     fprintf(pFile2,"NPredErr        : %+6.2E\n",toReturn.NPredErr);
     fprintf(pFile2,"NSign           : %+6.2E\n",toReturn.NSign);

     fclose(pFile2);
     return;
}





int JobIdToIndex(string JobId){
   for(unsigned int s=0;s<signals.size();s++){
      if(signals[s].Name==JobId)return s;
   }return -1;
}

double FindIntersection(TGraph* obs, TGraph* th, double Min, double Max, double Step, double ThUncertainty, bool debug){

   double Intersection = -1;

   double ThShift = 1.0-ThUncertainty;
   double PreviousX = Min;
   double PreviousV = obs->Eval(PreviousX, 0, "") - (ThShift * th->Eval(PreviousX, 0, "")) ;
   if(PreviousV>0)return -1;
   for(double x=Min+=Step;x<Max;x+=Step){                 
      double V = obs->Eval(x, 0, "") - (ThShift * th->Eval(x, 0, "") );
      if(debug){
         printf("%7.2f --> Obs=%6.2E  Th=%6.2E",x,obs->Eval(x, 0, ""),ThShift * th->Eval(x, 0, ""));
         if(V>=0)printf("   X\n");
         else printf("\n");
      }
      if(V<0){
         PreviousX = x;
         PreviousV = V;
      }else{
         Intersection = PreviousX;
      }
   }
   return Intersection;
}



int ReadXSection(string InputFile, double* Mass, double* XSec, double* Low, double* High, double* ErrLow, double* ErrHigh)
{
   FILE* pFile = fopen(InputFile.c_str(),"r");
   if(!pFile){ 
      printf("Not Found: %s\n",InputFile.c_str());
      return -1;
   }

   float tmpM, tmpX, tmpL, tmpH;
   
   int NPoints = 0;
   while ( ! feof (pFile) ){
     fscanf(pFile,"%f %E %E %E\n",&tmpM,&tmpX,&tmpH,&tmpL);
     Mass   [NPoints] = tmpM;
     XSec   [NPoints] = tmpX;
     Low    [NPoints] = tmpL;
     High   [NPoints] = tmpH;
     ErrLow [NPoints] = tmpX-tmpL;
     ErrHigh[NPoints] = tmpH-tmpX;
     NPoints++;

     //printf("%fGeV --> Error = %f\n", tmpM, 0.5*(tmpH-tmpL)/tmpX);
   }

   fclose(pFile);

   return NPoints;
}


TCutG* GetErrorBand(string name, int N, double* Mass, double* Low, double* High, double MinLow, double MaxHigh)
{
   TCutG* cutg = new TCutG(name.c_str(),2*N);
   cutg->SetFillColor(kGreen-7);
   for(int i=0;i<N;i++){
      double Min = std::max(Low[i],MinLow);
      cutg->SetPoint( i,Mass[i], Min);
   }
   for(int i=0;i<N;i++){
      double Max = std::min(High[N-1-i],MaxHigh);
      cutg->SetPoint(N+i,Mass[N-1-i], Max);
   }
   return cutg;
}

void DrawModelLimitWithBand(string InputPattern, string inputmodel)
{
   std::vector<string> Models;
   string modelname;
   if(inputmodel == "Gluinof1"){
      Models.push_back("Gluino300_f1");
      Models.push_back("Gluino400_f1");
      Models.push_back("Gluino500_f1");
      Models.push_back("Gluino600_f1");
      Models.push_back("Gluino700_f1");
      Models.push_back("Gluino800_f1");
      Models.push_back("Gluino900_f1");
      Models.push_back("Gluino1000_f1");
      Models.push_back("Gluino1100_f1");
      Models.push_back("Gluino1200_f1");
      modelname="gluino; 10% #tilde{g}g (NLO+NLL)";
   }
   else if(inputmodel == "Gluinof5"){
      Models.push_back("Gluino300_f5");
      Models.push_back("Gluino400_f5");
      Models.push_back("Gluino500_f5");
      Models.push_back("Gluino600_f5");
      Models.push_back("Gluino700_f5");
      Models.push_back("Gluino800_f5");
      Models.push_back("Gluino900_f5");
      Models.push_back("Gluino1000_f5");
      Models.push_back("Gluino1100_f5");
      Models.push_back("Gluino1200_f5");
      modelname="gluino; 50% #tilde{g}g (NLO+NLL)";
   }

   else if(inputmodel == "Gluinof10"){
     Models.push_back("Gluino300_f10");
     Models.push_back("Gluino400_f10");
     Models.push_back("Gluino500_f10");
     Models.push_back("Gluino600_f10");
     Models.push_back("Gluino700_f10");
     Models.push_back("Gluino800_f10");
     Models.push_back("Gluino900_f10");
     Models.push_back("Gluino1000_f10");
     Models.push_back("Gluino1100_f10");
     Models.push_back("Gluino1200_f10");
     modelname="gluino; 100% #tilde{g}g (NLO+NLL)";
   }

   else{cout<<"no model specified"<<endl;}

   unsigned int N = Models.size();
   stAllInfo Infos;double Mass[N], XSecTh[N], XSecExp[N],XSecObs[N], XSecExpUp[N],XSecExpDown[N],XSecExp2Up[N],XSecExp2Down[N];
   for(unsigned int i=0;i<N;i++){
      Infos = stAllInfo(InputPattern+"EXCLUSION/" + Models[i] +".txt");
      Mass[i]=Infos.Mass;
      XSecTh [i]=Infos.XSec_Th;
      XSecObs[i]=Infos.XSec_Obs;
      XSecExp[i]=Infos.XSec_Exp;
      XSecExpUp[i]=Infos.XSec_ExpUp;
      XSecExpDown[i]=Infos.XSec_ExpDown;
      XSecExp2Up[i]=Infos.XSec_Exp2Up;
      XSecExp2Down[i]=Infos.XSec_Exp2Down;
   }

   TGraph* graphtheory = new TGraph(N,Mass,XSecTh);
   TGraph* graphobs = new TGraph(N,Mass,XSecObs);
   TGraph* graphexp = new TGraph(N,Mass,XSecExp);
   TCutG*  ExpErr = GetErrorBand("ExpErr",N,Mass,XSecExpDown,XSecExpUp);
   TCutG*  Exp2SigmaErr = GetErrorBand("Exp2SigmaErr",N,Mass,XSecExp2Down,XSecExp2Up);

   graphtheory->SetLineStyle(3);
   graphtheory->SetFillColor(kBlue);
   graphexp->SetLineStyle(4); 
   graphexp->SetLineColor(kRed);
   graphexp->SetMarkerStyle(); 
   graphexp->SetMarkerSize(0.); 
   Exp2SigmaErr->SetFillColor(kYellow);
   Exp2SigmaErr->SetLineColor(kWhite);
   ExpErr->SetFillColor(kGreen);
   ExpErr->SetLineColor(kWhite);
   graphobs->SetLineColor(kBlack);
   graphobs->SetLineWidth(2);
   graphobs->SetMarkerColor(kBlack);
   graphobs->SetMarkerStyle(23);

   TCanvas* c1 = new TCanvas("c1", "c1",600,600);
   TMultiGraph* MG = new TMultiGraph();

   MG->Add(graphexp      ,"LP");
   MG->Add(graphobs      ,"LP");
   MG->Add(graphtheory      ,"L");
   MG->Draw("A");
   Exp2SigmaErr->Draw("f");
   ExpErr  ->Draw("f");
   MG->Draw("same");
   MG->SetTitle("");
   MG->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
   MG->GetYaxis()->SetTitle("#sigma (pb)");
   MG->GetYaxis()->SetTitleOffset(1.70);
   MG->GetYaxis()->SetRangeUser(PlotMinScale,PlotMaxScale);
   DrawPreliminary(IntegratedLuminosity);
   
   TLegend* LEG = new TLegend(0.40,0.65,0.8,0.90);
   string headerstr;
   headerstr = "95% C.L. Limits (SA Only)";
   LEG->SetHeader(headerstr.c_str());
   LEG->SetFillColor(0); 
   LEG->SetBorderSize(0);
   LEG->AddEntry(graphtheory,  modelname.c_str() ,"L");
   LEG->AddEntry(graphexp, "Expected"       ,"L");
   LEG->AddEntry(ExpErr, "Expected #pm 1#sigma","F");
   LEG->AddEntry(Exp2SigmaErr, "Expected #pm 2#sigma "       ,"F");
   LEG->AddEntry(graphobs, "Observed"       ,"LP");
   LEG->Draw();

   c1->SetLogy(true);


   SaveCanvas(c1,"Results/EXCLUSION/", string(inputmodel + "ExclusionLog"));
   delete c1;


}
std::vector<string> GetModels(string inputmodel)
{
   std::vector<string> Models;
   string modelname;
   if(inputmodel == "Gluinof1"){
      Models.push_back("Gluino300_f1");
      Models.push_back("Gluino400_f1");
      Models.push_back("Gluino500_f1");
      Models.push_back("Gluino600_f1");
      Models.push_back("Gluino700_f1");
      Models.push_back("Gluino800_f1");
      Models.push_back("Gluino900_f1");
      Models.push_back("Gluino1000_f1");
      Models.push_back("Gluino1100_f1");
      Models.push_back("Gluino1200_f1");
      modelname="gluino; 10% #tilde{g}g (NLO+NLL)";
   }
   else if(inputmodel == "Gluinof5"){
      Models.push_back("Gluino300_f5");
      Models.push_back("Gluino400_f5");
      Models.push_back("Gluino500_f5");
      Models.push_back("Gluino600_f5");
      Models.push_back("Gluino700_f5");
      Models.push_back("Gluino800_f5");
      Models.push_back("Gluino900_f5");
      Models.push_back("Gluino1000_f5");
      Models.push_back("Gluino1100_f5");
      Models.push_back("Gluino1200_f5");
     modelname="gluino; 50% #tilde{g}g (NLO+NLL)";
   }
   else if(inputmodel == "Gluinof10"){
     Models.push_back("Gluino300_f10");
     Models.push_back("Gluino400_f10");
     Models.push_back("Gluino500_f10");
     Models.push_back("Gluino600_f10");
     Models.push_back("Gluino700_f10");
     Models.push_back("Gluino800_f10");
     Models.push_back("Gluino900_f10");
     Models.push_back("Gluino1000_f10");
     Models.push_back("Gluino1100_f10");
     Models.push_back("Gluino1200_f10");
     modelname="gluino; 10% #tilde{g}g (NLO+NLL)";
   }

   else{cout<<"no model specified"<<endl;}
   return Models;

}
string GetModelName(string inputmodel)
{
   string modelname;
   if(inputmodel == "Gluinof1"){
      modelname="gluino; 10% #tilde{g}g";
   }
   else if(inputmodel == "Gluinof5"){
     modelname="gluino; 50% #tilde{g}g";
   }
   else if(inputmodel == "Gluinof10"){
     modelname="gluino; 100% #tilde{g}g";
   }
   else{cout<<"no model specified"<<endl;}
   return modelname;


}
