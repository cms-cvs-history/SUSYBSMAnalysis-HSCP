
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
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "tdrstyle.C"
#include "Analysis_PlotFunction.h"
#include "Analysis_Samples.h"
#include "Analysis_Global.h"

std::vector<stSignal> signals;

void Make2DPlot_Core(string ResultPattern);
void MakeCompPlot(string Histo, string DirName, string InputPattern1, string InputPattern2="", string InputPattern3="", string InputPattern4="");
//void CheckPredictionRescale(string Histo, string DirName, string InputPattern1, string InputPattern2="", string InputPattern3="", string InputPattern4="");
void CheckPredictionRescale(string Histo, string DirName, string InputPattern1);
void MakeHitSplit_Plot(string InputPattern);

int JobIdToIndex(string JobId);
TF1* GetMassLine(double M, bool MC);

void Make2DPlot(){
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.13);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
/*
     MakeHitSplit_Plot("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt+00/WPI+00/");
     MakeHitSplit_Plot("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt+00/WPI+00/");

     Make2DPlot_Core("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt+00/WPI+00/");
     Make2DPlot_Core("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt+00/WPI+00/");
     Make2DPlot_Core("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-20/WPI-30/");
     Make2DPlot_Core("SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt-35/WPI-35/");

     MakeCompPlot("", "TkOnlyClusterCleaning", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt-50/WPI+00/", "SplitMode2/MinHit01/Sele_dedxASmi/Mass_dedxCNPHarm2/Type0/WPPt-50/WPI+00/");
     MakeCompPlot("", "TkMuonClusterCleaning", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-50/WPI-50/", "SplitMode2/MinHit01/Sele_dedxASmi/Mass_dedxCNPHarm2/Type1/WPPt-50/WPI-50/");
*/

     CheckPredictionRescale("TkOnly_WP20_20", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt-20/WPI-20/");
     CheckPredictionRescale("TkOnly_WP15_15", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type0/WPPt-15/WPI-15/");


     CheckPredictionRescale("TkMuon_WP20_05", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-20/WPI-05/");
     CheckPredictionRescale("TkMuon_WP15_05", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-15/WPI-05/");
     CheckPredictionRescale("TkMuon_WP05_05", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-05/WPI-05/");
     CheckPredictionRescale("TkMuon_WP10_10", "Prediction", "SplitMode2/MinHit01/Sele_dedxSTASmi/Mass_dedxSTCNPHarm2/Type1/WPPt-10/WPI-10/");

}




void Make2DPlot_Core(string InputPattern){
   TCanvas* c1;
   TLegend* leg;
 
   GetSignalDefinition(signals);

   string Input = "Results/ANALYSE/" + InputPattern + "DumpHistos.root";
   string outpath = string("Results/PLOT/") + InputPattern;
   MakeDirectories(outpath);


   TFile* InputFile = new TFile(Input.c_str());
   TH1D* Data_Mass    = (TH1D*)GetObjectFromPath(InputFile, "Mass_Data");
   TH2D* Data_PIs     = (TH2D*)GetObjectFromPath(InputFile, "Data_AS_PIs");
   TH2D* Data_PIm     = (TH2D*)GetObjectFromPath(InputFile, "Data_AS_PIm");
   TH1D* Stop130_Mass = (TH1D*)GetObjectFromPath(InputFile, "Mass_Stop130");
   TH2D* Stop130_PIs  = (TH2D*)GetObjectFromPath(InputFile, "Stop130_AS_PIs");
   TH2D* Stop130_PIm  = (TH2D*)GetObjectFromPath(InputFile, "Stop130_AS_PIm");
   TH1D* Stop200_Mass = (TH1D*)GetObjectFromPath(InputFile, "Mass_Stop200");
   TH2D* Stop200_PIs  = (TH2D*)GetObjectFromPath(InputFile, "Stop200_AS_PIs");
   TH2D* Stop200_PIm  = (TH2D*)GetObjectFromPath(InputFile, "Stop200_AS_PIm");
   TH1D* Stop300_Mass = (TH1D*)GetObjectFromPath(InputFile, "Mass_Stop300");
   TH2D* Stop300_PIs  = (TH2D*)GetObjectFromPath(InputFile, "Stop300_AS_PIs");
   TH2D* Stop300_PIm  = (TH2D*)GetObjectFromPath(InputFile, "Stop300_AS_PIm");
   TH1D* Stop500_Mass = (TH1D*)GetObjectFromPath(InputFile, "Mass_Stop500");
   TH2D* Stop500_PIs  = (TH2D*)GetObjectFromPath(InputFile, "Stop500_AS_PIs");
   TH2D* Stop500_PIm  = (TH2D*)GetObjectFromPath(InputFile, "Stop500_AS_PIm");
   TH1D* Stop800_Mass = (TH1D*)GetObjectFromPath(InputFile, "Mass_Stop800");
   TH2D* Stop800_PIs  = (TH2D*)GetObjectFromPath(InputFile, "Stop800_AS_PIs");
   TH2D* Stop800_PIm  = (TH2D*)GetObjectFromPath(InputFile, "Stop800_AS_PIm");


//   Data_Mass    = (TH1D*) Data_Mass->Rebin(4);
   Stop130_Mass = (TH1D*) Stop130_Mass->Rebin(4);
   Stop200_Mass = (TH1D*) Stop200_Mass->Rebin(4);
   Stop300_Mass = (TH1D*) Stop300_Mass->Rebin(4);
   Stop500_Mass = (TH1D*) Stop500_Mass->Rebin(4);
   Stop800_Mass = (TH1D*) Stop800_Mass->Rebin(4);


   double Min = 1E-9;
   double Max = 1;

   c1 = new TCanvas("c1","c1", 800, 600);
//   c1->SetGridx(true);
   c1->SetGridy(true);
   Stop130_Mass->SetAxisRange(Min,Max,"Y");
   Stop130_Mass->SetTitle("");
   Stop130_Mass->SetStats(kFALSE);
   Stop130_Mass->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   Stop130_Mass->GetYaxis()->SetTitle("Entries in 8.4nb^{-1}");
   Stop130_Mass->SetLineWidth(2);
   Stop130_Mass->SetLineColor(Color[0]);
   Stop130_Mass->SetMarkerColor(Color[0]);
   Stop130_Mass->SetMarkerStyle(Marker[0]);
//   Stop130_Mass->SetMarkerSize(0.3);
   Stop130_Mass->Draw("E1");
   Stop200_Mass->Draw("E1 same");
   Stop200_Mass->SetLineWidth(2);
   Stop200_Mass->SetLineColor(Color[1]);
   Stop200_Mass->SetMarkerColor(Color[1]);
   Stop200_Mass->SetMarkerStyle(Marker[1]);
//   Stop200_Mass->SetMarkerSize(0.3);
   Stop300_Mass->Draw("E1 same");
   Stop300_Mass->SetLineWidth(2);
   Stop300_Mass->SetLineColor(Color[2]);
   Stop300_Mass->SetMarkerColor(Color[2]);
   Stop300_Mass->SetMarkerStyle(Marker[2]);
//   Stop300_Mass->SetMarkerSize(0.3);
   Stop500_Mass->Draw("E1 same");
   Stop500_Mass->SetLineWidth(2);
   Stop500_Mass->SetLineColor(Color[3]);
   Stop500_Mass->SetMarkerColor(Color[3]);
   Stop500_Mass->SetMarkerStyle(Marker[3]);
//   Stop500_Mass->SetMarkerSize(0.3);
   Stop800_Mass->Draw("E1 same");
   Stop800_Mass->SetLineWidth(2);
   Stop800_Mass->SetLineColor(Color[4]);
   Stop800_Mass->SetMarkerColor(Color[4]);
   Stop800_Mass->SetMarkerStyle(Marker[4]);
//   Stop800_Mass->SetMarkerSize(0.3);
   c1->SetLogy(true);


   TLine* lineStop130 = new TLine(130, Min, 130, Max);
   lineStop130->SetLineWidth(2);
   lineStop130->SetLineColor(Color[0]);
   lineStop130->SetLineStyle(2);
   lineStop130->Draw("same");

   TLine* lineStop200 = new TLine(200, Min, 200, Max);
   lineStop200->SetLineWidth(2);
   lineStop200->SetLineColor(Color[1]);
   lineStop200->SetLineStyle(2);
   lineStop200->Draw("same");

   TLine* lineStop300 = new TLine(300, Min, 300, Max);
   lineStop300->SetLineWidth(2);
   lineStop300->SetLineColor(Color[2]);
   lineStop300->SetLineStyle(2);
   lineStop300->Draw("same");

   TLine* lineStop500 = new TLine(500, Min, 500, Max);
   lineStop500->SetLineWidth(2);
   lineStop500->SetLineColor(Color[3]);
   lineStop500->SetLineStyle(2);
   lineStop500->Draw("same");

   TLine* lineStop800 = new TLine(800, Min, 800, Max);
   lineStop800->SetLineWidth(2);
   lineStop800->SetLineColor(Color[4]);
   lineStop800->SetLineStyle(2);
   lineStop800->Draw("same");

   leg = new TLegend(0.80,0.93,0.80 - 0.20,0.93 - 6*0.03);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(Stop130_Mass, "Stop130"   ,"P");
   leg->AddEntry(Stop200_Mass, "Stop200"   ,"P");
   leg->AddEntry(Stop300_Mass, "Stop300"   ,"P");
   leg->AddEntry(Stop500_Mass, "Stop500"   ,"P");
   leg->AddEntry(Stop800_Mass, "Stop800"   ,"P");
   leg->Draw();
   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "Stop_Mass");
   delete c1;

   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetLogz(true);
   c1->SetGridx(true);
   c1->SetGridy(true);
   Data_PIs->SetTitle("");
   Data_PIs->SetStats(kFALSE);
   Data_PIs->GetXaxis()->SetTitle("P (GeV/c)");
   Data_PIs->GetYaxis()->SetTitle("dE/dx discriminator");
   Data_PIs->SetMarkerSize (0.2);
   Data_PIs->SetMarkerColor(Color[4]);
   Data_PIs->SetFillColor(Color[4]);
   Data_PIs->Draw("COLZ");
   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "Data_PIs");
   delete c1;

   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetLogz(true);
   c1->SetGridx(true);
   c1->SetGridy(true);
   Data_PIm->SetTitle("");
   Data_PIm->SetStats(kFALSE);
   Data_PIm->GetXaxis()->SetTitle("P (GeV/c)");
   Data_PIm->GetYaxis()->SetTitle("dE/dx estimator (MeV/cm)");
   Data_PIm->SetAxisRange(0,25,"Y");
   Data_PIm->SetMarkerSize (0.2);
   Data_PIm->SetMarkerColor(Color[4]);
   Data_PIm->SetFillColor(Color[4]);
   Data_PIm->Draw("COLZ");
   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "Data_PIm");
   delete c1;



   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetGridx(true);
   c1->SetGridy(true);
   Stop800_PIs->SetTitle("");
   Stop800_PIs->SetStats(kFALSE);
   Stop800_PIs->GetXaxis()->SetTitle("P (GeV/c)");
   Stop800_PIs->GetYaxis()->SetTitle("dE/dx discriminator");
   Stop800_PIs->Scale(1/Stop800_PIs->Integral());
   Stop800_PIs->SetMarkerSize (0.2);
   Stop800_PIs->SetMarkerColor(Color[4]);
   Stop800_PIs->SetFillColor(Color[4]);
   Stop800_PIs->Draw("");
   Stop500_PIs->Scale(1/Stop500_PIs->Integral());
   Stop500_PIs->SetMarkerSize (0.2);
   Stop500_PIs->SetMarkerColor(Color[3]);
   Stop500_PIs->SetFillColor(Color[3]);
   Stop500_PIs->Draw("same");
   Stop300_PIs->Scale(1/Stop300_PIs->Integral());
   Stop300_PIs->SetMarkerSize (0.2);
   Stop300_PIs->SetMarkerColor(Color[2]);
   Stop300_PIs->SetFillColor(Color[2]);
   Stop300_PIs->Draw("same");
   Stop200_PIs->Scale(1/Stop200_PIs->Integral());
   Stop200_PIs->SetMarkerSize (0.2);
   Stop200_PIs->SetMarkerColor(Color[1]);
   Stop200_PIs->SetFillColor(Color[1]);
   Stop200_PIs->Draw("same");
   Stop130_PIs->Scale(1/Stop130_PIs->Integral());
   Stop130_PIs->SetMarkerSize (0.2);
   Stop130_PIs->SetMarkerColor(Color[0]);
   Stop130_PIs->SetFillColor(Color[0]);
   Stop130_PIs->Draw("same");

   leg = new TLegend(0.80,0.93,0.80 - 0.20,0.93 - 6*0.03);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(Stop130_PIs, "Stop130"   ,"F");
   leg->AddEntry(Stop200_PIs, "Stop200"   ,"F");
   leg->AddEntry(Stop300_PIs, "Stop300"   ,"F");
   leg->AddEntry(Stop500_PIs, "Stop500"   ,"F");
   leg->AddEntry(Stop800_PIs, "Stop800"   ,"F");
   leg->Draw();
   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "Stop_PIs");
   delete c1;



   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetGridx(true);
   c1->SetGridy(true);
   Stop800_PIm->SetTitle("");
   Stop800_PIm->SetStats(kFALSE);
   Stop800_PIm->GetXaxis()->SetTitle("P (GeV/c)");
   Stop800_PIm->GetYaxis()->SetTitle("dE/dx estimator (MeV/cm)");
   Stop800_PIm->SetAxisRange(0,25,"Y");
   Stop800_PIm->Scale(1/Stop800_PIm->Integral());
   Stop800_PIm->SetMarkerSize (0.2);
   Stop800_PIm->SetMarkerColor(Color[4]);
   Stop800_PIm->SetFillColor(Color[4]);
   Stop800_PIm->Draw("");
   Stop500_PIm->Scale(1/Stop500_PIm->Integral());
   Stop500_PIm->SetMarkerSize (0.2);
   Stop500_PIm->SetMarkerColor(Color[3]);
   Stop500_PIm->SetFillColor(Color[3]);
   Stop500_PIm->Draw("same");
   Stop300_PIm->Scale(1/Stop300_PIm->Integral());
   Stop300_PIm->SetMarkerSize (0.2);
   Stop300_PIm->SetMarkerColor(Color[2]);
   Stop300_PIm->SetFillColor(Color[2]);
   Stop300_PIm->Draw("same");
   Stop200_PIm->Scale(1/Stop200_PIm->Integral());
   Stop200_PIm->SetMarkerSize (0.2);
   Stop200_PIm->SetMarkerColor(Color[1]);
   Stop200_PIm->SetFillColor(Color[1]);
   Stop200_PIm->Draw("same");
   Stop130_PIm->Scale(1/Stop130_PIm->Integral());
   Stop130_PIm->SetMarkerSize (0.2);
   Stop130_PIm->SetMarkerColor(Color[0]);
   Stop130_PIm->SetFillColor(Color[0]);
   Stop130_PIm->Draw("same");

   TF1* Stop800Line = GetMassLine(800, true);
   Stop800Line->SetLineColor(kMagenta-7);
   Stop800Line->SetLineWidth(2);
   Stop800Line->Draw("same");
   TF1* Stop500Line = GetMassLine(500, true);
   Stop500Line->SetLineColor(kGreen-7);
   Stop500Line->SetLineWidth(2);
   Stop500Line->Draw("same");
   TF1* Stop300Line = GetMassLine(300, true);
   Stop300Line->SetLineColor(kGray+3);
   Stop300Line->SetLineWidth(2);
   Stop300Line->Draw("same");
   TF1* Stop200Line = GetMassLine(200, true);
   Stop200Line->SetLineColor(kBlue-7);
   Stop200Line->SetLineWidth(2);
   Stop200Line->Draw("same");
   TF1* Stop130Line = GetMassLine(130, true);
   Stop130Line->SetLineColor(kRed-7);
   Stop130Line->SetLineWidth(2);
   Stop130Line->Draw("same");

   leg = new TLegend(0.80,0.93,0.80 - 0.20,0.93 - 6*0.03);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(Stop130_PIm, "Stop130"   ,"F");
   leg->AddEntry(Stop200_PIm, "Stop200"   ,"F");
   leg->AddEntry(Stop300_PIm, "Stop300"   ,"F");
   leg->AddEntry(Stop500_PIm, "Stop500"   ,"F");
   leg->AddEntry(Stop800_PIm, "Stop800"   ,"F");
   leg->Draw();
   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "Stop_PIm");
   delete c1;
}

void MakeCompPlot(string Histo, string DirName, string InputPattern1, string InputPattern2, string InputPattern3, string InputPattern4){
   string outpath = string("Results/PLOT/") + DirName;
   MakeDirectories(outpath);

   TH1D* Histo1;
   TH1D* Histo2;
   TFile* InputFile1;
   TFile* InputFile2;
   string Input;

   TH1** Histos = new TH1*[10];
   std::vector<string> legend;
   TCanvas* c1;


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Data_BS_Is"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(4);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Data_BS_Is"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx discriminator", "#Tracks", 0,0, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);   
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Is_Data");
   delete c1;

   printf("INTEGRAL OF DATA from 0.5 to 1 --> %6.2E\n",Histo1->Integral(Histo1->GetXaxis()->FindBin(0.5), Histo1->GetXaxis()->FindBin(1.0)));
   printf("INTEGRAL OF DATA from 0.5 to 1 --> %6.2E\n",Histo2->Integral(Histo2->GetXaxis()->FindBin(0.5), Histo2->GetXaxis()->FindBin(1.0)));

   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Gluino200_BS_Is"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(4);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Gluino200_BS_Is"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx discriminator", "#Tracks", 0,0, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Is_Gluino200");
   delete c1;

   printf("INTEGRAL OF S from 0.5 to 1 --> %6.2E\n",Histo1->Integral(Histo1->GetXaxis()->FindBin(0.5), Histo1->GetXaxis()->FindBin(1.0)));
   printf("INTEGRAL OF S from 0.5 to 1 --> %6.2E\n",Histo2->Integral(Histo2->GetXaxis()->FindBin(0.5), Histo2->GetXaxis()->FindBin(1.0)));


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Gluino900_BS_Is"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(4);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Gluino900_BS_Is"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(4);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx discriminator", "#Tracks", 0,0, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Is_Gluino900");
   delete c1;





   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Data_BS_Im"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(1);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Data_BS_Im"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(1);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx estimator (MeV/cm)", "#Tracks", 0,25, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);   
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Im_Data");
   delete c1;


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Gluino200_BS_Im"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(1);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Gluino200_BS_Im"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(1);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx estimator (MeV/cm)", "#Tracks", 0,25, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Im_Gluino200");
   delete c1;


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   Histo1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Gluino900_BS_Im"))->Clone("Hist1");
   Histo1 = (TH1D*)Histo1->Rebin(1);

   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   Histo2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Gluino900_BS_Im"))->Clone("Hist2");
   Histo2 = (TH1D*)Histo2->Rebin(1);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)Histo1;                      legend.push_back("With ClusterCleaning");
   Histos[1] = (TH1*)Histo2;                      legend.push_back("Without Cluster Cleaning");
   Histo1->SetMaximum(Histo1->GetMaximum()*1.1);
   DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1",  "dE/dx estimator (Mev/cm)", "#Tracks", 0,25, 0,0);
   Histo2->SetLineStyle(2);
   DrawLegend((TObject**)Histos,legend,"","LP",0.69, 0.93, 0.30, 0.05);
   c1->SetLogy(true);
   c1->Modified();
   DrawPreliminary(-1);
   SaveCanvas(c1,outpath,"Im_Gluino900");
   delete c1;
}


void CheckPredictionRescale(string Histo, string DirName, string InputPattern1){
   string outpath = string("Results/PLOT/") + DirName;
   MakeDirectories(outpath);

   TH1D* Histo1;
   TFile* InputFile1;
   string Input;

   TH1** Histos = new TH1*[10];
   std::vector<string> legend;
   TCanvas* c1;


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   TH1D* Pred1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_Pred"))->Clone("Pred1");
   Pred1 = (TH1D*)Pred1->Rebin(4);
   TH1D* Data1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_Data"))->Clone("Data1");
   Data1 = (TH1D*)Data1->Rebin(4);
   TH1D* MCTr1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_MCTr"))->Clone("MCTr1");
   MCTr1 = (TH1D*)MCTr1->Rebin(4);
   MCTr1->Scale(Data1->Integral()/MCTr1->Integral());


   double RescaleFactor = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin( 75.0))/Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin( 75.0));
//   double RescaleFactor = 1.5;
   TH1D* Resc1 = (TH1D*)(Pred1->Clone("Resc1"));
   TH1D* Resc2 = (TH1D*)(Pred1->Clone("Resc2"));
   Resc1->Scale(RescaleFactor);
   printf("SCALE FACTOR = %f\n",RescaleFactor);
   Resc2->Scale(1.5);


   double D,P,R;
   printf("%s\n",InputPattern1.c_str());
   D = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin(999.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin(999.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin( 0.0), Resc1->GetXaxis()->FindBin(999.0));
//   printf("INTEGRAL in [ 0,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin( 75.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin( 75.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin( 0.0), Resc1->GetXaxis()->FindBin( 75.0));
   printf("INTEGRAL in [ 0, 75] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data1->Integral(Data1->GetXaxis()->FindBin(75.0), Data1->GetXaxis()->FindBin(999.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin(75.0), Pred1->GetXaxis()->FindBin(999.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin(75.0), Resc1->GetXaxis()->FindBin(999.0));
   printf("INTEGRAL in [75,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);


   D = Data1->Integral(Data1->GetXaxis()->FindBin(75.0), Data1->GetXaxis()->FindBin(100.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin(75.0), Pred1->GetXaxis()->FindBin(100.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin(75.0), Resc1->GetXaxis()->FindBin(100.0));
   printf("INTEGRAL in [75,100] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);


   D = Data1->Integral(Data1->GetXaxis()->FindBin(100.0), Data1->GetXaxis()->FindBin(999.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin(100.0), Pred1->GetXaxis()->FindBin(999.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin(100.0), Resc1->GetXaxis()->FindBin(999.0));
   printf("INTEGRAL in [100,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);


   double Max = std::max(Data1->GetMaximum(), Pred1->GetMaximum());
   Max        = std::max(MCTr1->GetMaximum(), Max);
   Max       *= 1.5;
   double Min = std::min(0.01,Pred1->GetMaximum());
   Min       *= 0.05;

   TLegend* leg;
   c1 = new TCanvas("c1","c1,",600,600);

   MCTr1->GetXaxis()->SetNdivisions(505);
   MCTr1->SetTitle("");
   MCTr1->SetStats(kFALSE);
   MCTr1->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   MCTr1->GetYaxis()->SetTitle("#Tracks");
   MCTr1->GetYaxis()->SetTitleOffset(1.50);
   MCTr1->SetLineColor(39);
   MCTr1->SetFillColor(64);
   MCTr1->SetMarkerSize(1);
   MCTr1->SetMarkerStyle(1);
   MCTr1->SetMarkerColor(39);
   MCTr1->SetMaximum(Max);
   MCTr1->SetMinimum(Min);
   MCTr1->Draw("HIST");
   TH1D* MCTr1Err = (TH1D*)MCTr1->Clone("MCTr1_Mass_Err");
   MCTr1Err->SetLineColor(46);
   MCTr1Err->Draw("E1 same");

   Pred1->SetMarkerStyle(21);
   Pred1->SetMarkerColor(8);
   Pred1->SetMarkerSize(0.5);
   Pred1->SetLineColor(8);
   Pred1->SetFillColor(0);
   Pred1->Draw("E1 same");

   Resc1->SetMarkerStyle(22);
   Resc1->SetMarkerColor(2);
   Resc1->SetMarkerSize(0.5);
   Resc1->SetLineColor(2);
   Resc1->SetFillColor(0);
   Resc1->Draw("E1 same");

   Data1->SetMarkerStyle(20);
   Data1->SetMarkerColor(1);
   Data1->SetMarkerSize(0.5);
   Data1->SetLineColor(1);
   Data1->SetFillColor(0);
   Data1->Draw("E1 same");

   leg = new TLegend(0.79,0.93,0.49,0.73);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(MCTr1, "MC"          ,"F");
   leg->AddEntry(Pred1, "Prediction Absolute"  ,"P");
   leg->AddEntry(Resc1, "Prediction Rescaled"  ,"P");
   leg->AddEntry(Data1, "Data"        ,"P");
   leg->Draw();

   DrawPreliminary(-1);
   c1->SetLogy(true);
   SaveCanvas(c1, outpath, Histo);

//   c1->SetLogy(true);
//   SaveCanvas(c1, outpath, Histo+"B");
   delete c1;

  TH1D* Ratio1       = (TH1D*)Resc1->Clone("Ratio1");
  TH1D* Ratio2       = (TH1D*)Resc2->Clone("Ratio2");
  TH1D* DataWithStat = (TH1D*)Data1->Clone("DataWithStat");
  for(unsigned int i=0;i<Ratio1->GetNbinsX();i++){
     if(Data1->GetBinContent(i)<2){
        DataWithStat->SetBinContent(i,0);
        DataWithStat->SetBinError(i,0);
     }
  }
  Ratio1->Divide(DataWithStat);
  Ratio2->Divide(DataWithStat);


  /*
  TH1D* Ratio1 = (TH1D*)Data1->Clone("Ratio1");
  for(unsigned int i=0;i<Ratio1->GetNbinsX();i++){
     if(Resc1->GetBinContent(i)>0 && Data1->GetBinContent(i)>1){
        Ratio1->SetBinContent(i,Data1->GetBinContent(i)/Resc1->GetBinContent(i));
        Ratio1->SetBinError(i,Ratio1->GetBinContent(i)*0.5);
     }else{
        Ratio1->SetBinContent(i,0);
        Ratio1->SetBinError(i,0);
     }
  }
  */

   c1 = new TCanvas("c1","c1,",600,600);
   Ratio1->SetAxisRange(0,500,"X");
   Ratio1->SetAxisRange(0,2,"Y");
   Ratio1->SetTitle("");
   Ratio1->SetStats(kFALSE);
   Ratio1->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   Ratio1->GetYaxis()->SetTitle("#Predicted / #Observed Tracks");
   Ratio1->GetYaxis()->SetTitleOffset(1.50);
   Ratio1->SetMarkerStyle(21);
   Ratio1->SetMarkerColor(4);
   Ratio1->SetMarkerSize(1);
   Ratio1->SetLineColor(4);
   Ratio1->SetFillColor(0);
   Ratio1->Draw("E1");

  
   TBox* b = new TBox(0,0.5,520,1.5);
   b->SetFillStyle(3003);
   b->SetFillColor(8);
   b->Draw("same");

   TLine* l = new TLine(0,1.0,520,1.0);
   l->Draw("same");

   Ratio1->Draw("same E1");

   Ratio2->SetAxisRange(0,500,"X");
   Ratio2->SetTitle("");
   Ratio2->SetStats(kFALSE);
   Ratio2->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   Ratio2->GetYaxis()->SetTitle("#Predicted / #Observed Tracks");
   Ratio2->GetYaxis()->SetTitleOffset(1.50);
   Ratio2->SetMarkerStyle(22);
   Ratio2->SetMarkerColor(2);
   Ratio2->SetMarkerSize(1);
   Ratio2->SetLineColor(2);
   Ratio2->SetFillColor(0);
   Ratio2->Draw("E1 same");

   leg = new TLegend(0.79,0.93,0.54,0.75);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("Rescale Factor");   
   char buf[255];
   sprintf(buf,"%4.2f",RescaleFactor);
   leg->AddEntry(Ratio1, buf     ,"P");
   leg->AddEntry(Ratio2, "1.50"  ,"P");
   leg->Draw();


   DrawPreliminary(-1);  
   SaveCanvas(c1, outpath, Histo + "Rescale1");
   delete c1;
}


/*
void CheckPredictionRescale(string Histo, string DirName, string InputPattern1, string InputPattern2, string InputPattern3, string InputPattern4){
   string outpath = string("Results/PLOT/") + DirName;
   MakeDirectories(outpath);

   TH1D* Histo1;
   TH1D* Histo2;
   TFile* InputFile1;
   TFile* InputFile2;
   string Input;

   TH1** Histos = new TH1*[10];
   std::vector<string> legend;
   TCanvas* c1;


   Input = string("Results/ANALYSE/") + InputPattern1 + "DumpHistos.root";
   InputFile1 = new TFile(Input.c_str());
   TH1D* Pred1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_Pred"))->Clone("Pred1");
   Pred1 = (TH1D*)Pred1->Rebin(4);
   TH1D* Data1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_Data"))->Clone("Data1");
   Data1 = (TH1D*)Data1->Rebin(4);
   TH1D* MCTr1 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile1, "Mass_MCTr"))->Clone("MCTr1");
   MCTr1 = (TH1D*)MCTr1->Rebin(4);
   MCTr1->Scale(Data1->Integral()/MCTr1->Integral());


   Input = string("Results/ANALYSE/") + InputPattern2 + "DumpHistos.root";
   InputFile2 = new TFile(Input.c_str());
   TH1D* Pred2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Mass_Pred"))->Clone("Pred2");
   Pred2 = (TH1D*)Pred2->Rebin(4);
   TH1D* Data2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Mass_Data"))->Clone("Data2");
   Data2 = (TH1D*)Data2->Rebin(4);
   TH1D* MCTr2 = (TH1D*)((TH1D*)GetObjectFromPath(InputFile2, "Mass_MCTr"))->Clone("MCTr2");
   MCTr2 = (TH1D*)MCTr2->Rebin(4);
   MCTr2->Scale(Data2->Integral()/MCTr2->Integral());


//   double RescaleFactor = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin( 75.0))/Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin( 75.0));
   double RescaleFactor = 1.5;
   TH1D* Resc1 = (TH1D*)(Pred1->Clone("Resc1"));
   Resc1->Scale(RescaleFactor);
   TH1D* Resc2 = (TH1D*)(Pred2->Clone("Resc2"));
   Resc2->Scale(RescaleFactor);

   printf("SCALE FACTOR = %f\n",RescaleFactor);

   double D,P,R;
   printf("%s\n",InputPattern1.c_str());
   D = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin(999.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin(999.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin( 0.0), Resc1->GetXaxis()->FindBin(999.0));
//   printf("INTEGRAL in [ 0,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data1->Integral(Data1->GetXaxis()->FindBin( 0.0), Data1->GetXaxis()->FindBin( 75.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin( 0.0), Pred1->GetXaxis()->FindBin( 75.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin( 0.0), Resc1->GetXaxis()->FindBin( 75.0));
   printf("INTEGRAL in [ 0, 75] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data1->Integral(Data1->GetXaxis()->FindBin(75.0), Data1->GetXaxis()->FindBin(999.0));
   P = Pred1->Integral(Pred1->GetXaxis()->FindBin(75.0), Pred1->GetXaxis()->FindBin(999.0));
   R = Resc1->Integral(Resc1->GetXaxis()->FindBin(75.0), Resc1->GetXaxis()->FindBin(999.0));
   printf("INTEGRAL in [75,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   printf("%s\n",InputPattern2.c_str());
   D = Data2->Integral(Data2->GetXaxis()->FindBin( 0.0), Data2->GetXaxis()->FindBin(999.0));
   P = Pred2->Integral(Pred2->GetXaxis()->FindBin( 0.0), Pred2->GetXaxis()->FindBin(999.0));
   R = Resc2->Integral(Resc2->GetXaxis()->FindBin( 0.0), Resc2->GetXaxis()->FindBin(999.0));
//   printf("INTEGRAL in [ 0,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data2->Integral(Data2->GetXaxis()->FindBin( 0.0), Data2->GetXaxis()->FindBin( 75.0));
   P = Pred2->Integral(Pred2->GetXaxis()->FindBin( 0.0), Pred2->GetXaxis()->FindBin( 75.0));
   R = Resc2->Integral(Resc2->GetXaxis()->FindBin( 0.0), Resc2->GetXaxis()->FindBin( 75.0));
   printf("INTEGRAL in [ 0, 75] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   D = Data2->Integral(Data2->GetXaxis()->FindBin(75.0), Data2->GetXaxis()->FindBin(999.0));
   P = Pred2->Integral(Pred2->GetXaxis()->FindBin(75.0), Pred2->GetXaxis()->FindBin(999.0));
   R = Resc2->Integral(Resc2->GetXaxis()->FindBin(75.0), Resc2->GetXaxis()->FindBin(999.0));
   printf("INTEGRAL in [75,999] --> P = %9.3f R = %9.3f D = %9.3f  D/P = %8.3f  D/R = %8.3f\n", P, R, D, D/P, D/R);

   double Max = std::max(Data1->GetMaximum(), Pred1->GetMaximum());
   Max        = std::max(MCTr1->GetMaximum(), Max);
   Max        = std::max(Data2->GetMaximum(), Max);
   Max        = std::max(MCTr2->GetMaximum(), Max);
   Max       *= 1.5;
   double Min = std::min(0.01,Pred1->GetMaximum());
   Min        = std::min( Min,Pred2->GetMaximum());
   Min       *= 0.05;

   TLegend* leg;
   c1 = new TCanvas("c1","c1,",1200,600);
   c1->Divide(2,1);
   c1->cd(1);

   MCTr1->GetXaxis()->SetNdivisions(505);
   MCTr1->SetTitle("");
   MCTr1->SetStats(kFALSE);
   MCTr1->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   MCTr1->GetYaxis()->SetTitle("#Tracks");
   MCTr1->GetYaxis()->SetTitleOffset(1.50);
   MCTr1->SetLineColor(39);
   MCTr1->SetFillColor(64);
   MCTr1->SetMarkerSize(1);
   MCTr1->SetMarkerStyle(1);
   MCTr1->SetMarkerColor(39);
   MCTr1->SetMaximum(Max);
   MCTr1->SetMinimum(Min);
   MCTr1->Draw("HIST");
   TH1D* MCTr1Err = (TH1D*)MCTr1->Clone("MCTr1_Mass_Err");
   MCTr1Err->SetLineColor(46);
   MCTr1Err->Draw("E1 same");

   Pred1->SetMarkerStyle(21);
   Pred1->SetMarkerColor(8);
   Pred1->SetMarkerSize(0.5);
   Pred1->SetLineColor(8);
   Pred1->SetFillColor(0);
   Pred1->Draw("E1 same");

   Resc1->SetMarkerStyle(22);
   Resc1->SetMarkerColor(2);
   Resc1->SetMarkerSize(0.5);
   Resc1->SetLineColor(2);
   Resc1->SetFillColor(0);
   Resc1->Draw("E1 same");

   Data1->SetMarkerStyle(20);
   Data1->SetMarkerColor(1);
   Data1->SetMarkerSize(0.5);
   Data1->SetLineColor(1);
   Data1->SetFillColor(0);
   Data1->Draw("E1 same");

   leg = new TLegend(0.79,0.93,0.49,0.73);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(MCTr1, "MC"          ,"F");
   leg->AddEntry(Pred1, "Prediction Absolute"  ,"P");
   leg->AddEntry(Resc1, "Prediction Rescaled"  ,"P");
   leg->AddEntry(Data1, "Data"        ,"P");
   leg->Draw();

   (c1->cd(1))->SetLogy(true);
   //SaveCanvas(c1, outpath, Histo+"B");
   //delete c1;

   //c1 = new TCanvas("c1","c1,",600,600);
   c1->cd(2);

   MCTr2->GetXaxis()->SetNdivisions(505);
   MCTr2->SetTitle("");
   MCTr2->SetStats(kFALSE);
   MCTr2->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   MCTr2->GetYaxis()->SetTitle("#Tracks");
   MCTr2->GetYaxis()->SetTitleOffset(1.50);
   MCTr2->SetLineColor(39);
   MCTr2->SetFillColor(64);
   MCTr2->SetMarkerSize(1);
   MCTr2->SetMarkerStyle(1);
   MCTr2->SetMarkerColor(39);
   MCTr2->SetMaximum(Max);
   MCTr2->SetMinimum(Min);
   MCTr2->Draw("HIST");
   TH1D* MCTr2Err = (TH1D*)MCTr2->Clone("MCTr2_Mass_Err");
   MCTr2Err->SetLineColor(46);
   MCTr2Err->Draw("E1 same");

   Pred2->SetMarkerStyle(21);
   Pred2->SetMarkerColor(8);
   Pred2->SetMarkerSize(0.5);
   Pred2->SetLineColor(8);
   Pred2->SetFillColor(0);
   Pred2->Draw("E1 same");

   Resc2->SetMarkerStyle(22);
   Resc2->SetMarkerColor(2);
   Resc2->SetMarkerSize(0.5);
   Resc2->SetLineColor(2);
   Resc2->SetFillColor(0);
   Resc2->Draw("E1 same");

   Data2->SetMarkerStyle(20);
   Data2->SetMarkerColor(1);
   Data2->SetMarkerSize(0.5);
   Data2->SetLineColor(1);
   Data2->SetFillColor(0);
   Data2->Draw("E1 same");

   leg = new TLegend(0.79,0.93,0.49,0.73);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(MCTr2, "MC"          ,"F");
   leg->AddEntry(Pred2, "Prediction Absolute"  ,"P");
   leg->AddEntry(Resc2, "Prediction Rescaled"  ,"P");
   leg->AddEntry(Data2, "Data"        ,"P");
   leg->Draw();

   (c1->cd(2))->SetLogy(true);
   SaveCanvas(c1, outpath, Histo);

//   c1->SetLogy(true);
//   SaveCanvas(c1, outpath, Histo+"B");
   delete c1;






  TH1D* Ratio1 = (TH1D*)Data1->Clone("Ratio1");
  TH1D* Ratio2 = (TH1D*)Data1->Clone("Ratio2");
  for(unsigned int i=0;i<Ratio1->GetNbinsX();i++){
     if(Resc1->GetBinContent(i)>0 && Data1->GetBinContent(i)>1){
        Ratio1->SetBinContent(i,Data1->GetBinContent(i)/Resc1->GetBinContent(i));
        Ratio1->SetBinError(i,Ratio1->GetBinContent(i)*0.5);
     }else{
        Ratio1->SetBinContent(i,0);
        Ratio1->SetBinError(i,0);
     }

     if(Resc2->GetBinContent(i)>0 && Data2->GetBinContent(i)>1){
        Ratio2->SetBinContent(i,Data2->GetBinContent(i)/Resc2->GetBinContent(i));
        Ratio2->SetBinError(i,Ratio2->GetBinContent(i)*0.5);
     }else{
        Ratio2->SetBinContent(i,0);
        Ratio2->SetBinError(i,0);
     }
  }

   c1 = new TCanvas("c1","c1,",600,600);
   Ratio1->SetAxisRange(0,500,"X");
   Ratio1->SetTitle("");
   Ratio1->SetStats(kFALSE);
   Ratio1->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   Ratio1->GetYaxis()->SetTitle("#Observed / #Predicted Tracks");
   Ratio1->GetYaxis()->SetTitleOffset(1.50);
   Ratio1->SetMarkerStyle(22);
   Ratio1->SetMarkerColor(2);
   Ratio1->SetMarkerSize(1);
   Ratio1->SetLineColor(2);
   Ratio1->SetFillColor(0);
   Ratio1->Draw("E1");
   SaveCanvas(c1, outpath, Histo + "Rescale1");
   delete c1;



   c1 = new TCanvas("c1","c1,",600,600);
   Ratio2->SetAxisRange(0,500,"X");
   Ratio2->SetTitle("");
   Ratio2->SetStats(kFALSE);
   Ratio2->GetXaxis()->SetTitle("Reconstructed Mass (GeV/c^{2})");
   Ratio2->GetYaxis()->SetTitle("#Observed / #Predicted Tracks");
   Ratio2->GetYaxis()->SetTitleOffset(1.50);
   Ratio2->SetMarkerStyle(22);
   Ratio2->SetMarkerColor(2);
   Ratio2->SetMarkerSize(1);
   Ratio2->SetLineColor(2);
   Ratio2->SetFillColor(0);
   Ratio2->Draw("E1");
   SaveCanvas(c1, outpath, Histo + "Rescale2");
   delete c1;
}
*/


void MakeHitSplit_Plot(string InputPattern){
   TCanvas* c1;
   TLegend* leg;
 
   GetSignalDefinition(signals);

   string Input = "Results/ANALYSE/" + InputPattern + "DumpHistos.root";
   string outpath = string("Results/PLOT/") + InputPattern;
   MakeDirectories(outpath);


   TFile* InputFile = new TFile(Input.c_str());
   TH1D* Data_05_I  = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_I_Data_SSHit05");
   TH1D* Data_10_I  = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_I_Data_SSHit10");
   TH1D* Data_15_I  = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_I_Data_SSHit15");
   TH1D* Data_20_I  = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_I_Data_SSHit20");
   Data_05_I = (TH1D*) Data_05_I->Rebin(100);
   Data_10_I = (TH1D*) Data_10_I->Rebin(100);
   Data_15_I = (TH1D*) Data_15_I->Rebin(100);
   Data_20_I = (TH1D*) Data_20_I->Rebin(100);
   Data_05_I->Scale(1.0/Data_05_I->Integral());
   Data_10_I->Scale(1.0/Data_10_I->Integral());
   Data_15_I->Scale(1.0/Data_15_I->Integral());
   Data_20_I->Scale(1.0/Data_20_I->Integral());

   TH1D* Data_05_Pt = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_Pt_Data_SSHit05");
   TH1D* Data_10_Pt = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_Pt_Data_SSHit10");
   TH1D* Data_15_Pt = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_Pt_Data_SSHit15");
   TH1D* Data_20_Pt = (TH1D*)GetObjectFromPath(InputFile, "CutFinder_Pt_Data_SSHit20");
   Data_05_Pt = (TH1D*) Data_05_Pt->Rebin(50);
   Data_10_Pt = (TH1D*) Data_10_Pt->Rebin(50);
   Data_15_Pt = (TH1D*) Data_15_Pt->Rebin(50);
   Data_20_Pt = (TH1D*) Data_20_Pt->Rebin(50);
   Data_05_Pt->Scale(1.0/Data_05_Pt->Integral());
   Data_10_Pt->Scale(1.0/Data_10_Pt->Integral());
   Data_15_Pt->Scale(1.0/Data_15_Pt->Integral());
   Data_20_Pt->Scale(1.0/Data_20_Pt->Integral());

   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetGridy(true);
   Data_05_I->SetTitle("");
   Data_05_I->SetStats(kFALSE);
   Data_05_I->GetXaxis()->SetTitle("dE/dx discriminator");
   Data_05_I->GetYaxis()->SetTitle("arbitrary units");
   Data_05_I->SetLineWidth(2);
   Data_05_I->SetLineColor(Color[0]);
   Data_05_I->SetMarkerColor(Color[0]);
   Data_05_I->SetMarkerStyle(Marker[0]);
   Data_05_I->Draw("E1");
   Data_05_I->Draw("E1 same");
   Data_10_I->SetLineWidth(2);
   Data_10_I->SetLineColor(Color[1]);
   Data_10_I->SetMarkerColor(Color[1]);
   Data_10_I->SetMarkerStyle(Marker[1]);
   Data_10_I->Draw("E1 same");
   Data_15_I->SetLineWidth(2);
   Data_15_I->SetLineColor(Color[2]);
   Data_15_I->SetMarkerColor(Color[2]);
   Data_15_I->SetMarkerStyle(Marker[2]);
   Data_15_I->Draw("E1 same");
   Data_20_I->SetLineWidth(2);
   Data_20_I->SetLineColor(Color[3]);
   Data_20_I->SetMarkerColor(Color[3]);
   Data_20_I->SetMarkerStyle(Marker[3]);
   Data_20_I->Draw("E1 same");
   c1->SetLogy(true);

   leg = new TLegend(0.80,0.93,0.80 - 0.20,0.93 - 6*0.03);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(Data_05_I, "05 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_10_I, "10 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_15_I, "15 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_20_I, "20 dE/dx Hits"   ,"P");
   leg->Draw();

   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "IDistribution");
   delete c1;


   c1 = new TCanvas("c1","c1", 800, 600);
   c1->SetGridy(true);
   Data_05_Pt->SetAxisRange(0,200,"X");
   Data_05_Pt->SetTitle("");
   Data_05_Pt->SetStats(kFALSE);
   Data_05_Pt->GetXaxis()->SetTitle("Pt (GeV/c)");
   Data_05_Pt->GetYaxis()->SetTitle("arbitrary units");
   Data_05_Pt->SetLineWidth(2);
   Data_05_Pt->SetLineColor(Color[0]);
   Data_05_Pt->SetMarkerColor(Color[0]);
   Data_05_Pt->SetMarkerStyle(Marker[0]);
   Data_05_Pt->Draw("E1");
   Data_05_Pt->Draw("E1 same");
   Data_10_Pt->SetLineWidth(2);
   Data_10_Pt->SetLineColor(Color[1]);
   Data_10_Pt->SetMarkerColor(Color[1]);
   Data_10_Pt->SetMarkerStyle(Marker[1]);
   Data_10_Pt->Draw("E1 same");
   Data_15_Pt->SetLineWidth(2);
   Data_15_Pt->SetLineColor(Color[2]);
   Data_15_Pt->SetMarkerColor(Color[2]);
   Data_15_Pt->SetMarkerStyle(Marker[2]);
   Data_15_Pt->Draw("E1 same");
   Data_20_Pt->SetLineWidth(2);
   Data_20_Pt->SetLineColor(Color[3]);
   Data_20_Pt->SetMarkerColor(Color[3]);
   Data_20_Pt->SetMarkerStyle(Marker[3]);
   Data_20_Pt->Draw("E1 same");
   c1->SetLogy(true);

   leg = new TLegend(0.80,0.93,0.80 - 0.20,0.93 - 6*0.03);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(Data_05_Pt, "05 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_10_Pt, "10 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_15_Pt, "15 dE/dx Hits"   ,"P");
   leg->AddEntry(Data_20_Pt, "20 dE/dx Hits"   ,"P");
   leg->Draw();

   DrawPreliminary(-1);
   SaveCanvas(c1, outpath, "PtDistribution");
   delete c1;

}


int JobIdToIndex(string JobId){
   for(unsigned int s=0;s<signals.size();s++){
      if(signals[s].Name==JobId)return s;
   }return -1;
}

TF1* GetMassLine(double M, bool MC)
{   
   double K = dEdxK_Data[dEdxMassIndex];
   double C = dEdxC_Data[dEdxMassIndex];
   if(MC){
      K = dEdxK_MC[dEdxMassIndex];
      C = dEdxC_MC[dEdxMassIndex];
   }
   
   double BetaMax = 0.9;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));
   
   double BetaMin = 0.2;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));

   TF1* MassLine = new TF1("MassLine","[2] + ([0]*[0]*[1])/(x*x)", PMin, PMax);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, K);
   MassLine->SetParameter(2, C);
   MassLine->SetLineWidth(2);
   return MassLine;
}
