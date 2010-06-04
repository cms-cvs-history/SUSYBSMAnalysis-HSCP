
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
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "tdrstyle.C"


namespace reco    { class Vertex; class Track; class GenParticle;}
namespace susybsm { class HSCParticle;}
namespace fwlite  { class ChainEvent;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;

#endif


#include "Analysis_Global.h"
#include "Analysis_Samples.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"


/////////////////////////// FUNCTION DECLARATION /////////////////////////////


void FillArray(int HitIndex, int EtaIndex, double* Array, double value);
void FillHisto(int HitIndex, int EtaIndex, TH1D**  Histo, double value                , double weight=1);
void FillHisto(int HitIndex, int EtaIndex, TH2D**  Histo, double value1, double value2, double weight=1);

double deltaR(double eta1, double phi1, double eta2, double phi2);

void DrawDEDXVsP(TH2* Histos, char* Title,  char* Xlegend, char* Ylegend, double xmin, double xmax, double ymin, double ymax, bool DrawMassLine);


double CutFromEfficiency(TH1* Histo, double Efficiency, bool DoesKeepLeft=false);
double Efficiency(TH1* Histo, double CutX);
double Efficiency(TH2* Histo, double CutX, double CutY);

double GetMass(double P, double I, bool MC=false);
TF1* GetMassLine(double M, bool MC=false);


void GetIndices(int NOM, double Eta, int& HitIndex, int& EtaIndex);
int GetCutIndex(int HitIndex, int EtaIndex);
void GetNameFromIndex(char* NameExt, int index);

//void PrintEventInRange(double min, double max,FILE* pFile=stdout);

void Find_WorkingPoint(char* SavePath);
void CompleteAnalysis(char* SavePath);
void Merge_Map(char* SavePath);
void Analysis_Step2();
void Analysis_Step3();
void Analysis_Step4(char* SavePath);
void Analysis_Step5(char* SavePath);

void InitHistos();
void CleanUpHistos();

double DistToHSCP     (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest);
bool   isGoodCandidate(const susybsm::HSCParticle& hscp, const reco::Vertex& vertex, double PtCut=0, double ICut=0, stPlots* st=NULL);
void DumpCandidateInfo(const susybsm::HSCParticle& hscp, const reco::Vertex& vertex, const fwlite::ChainEvent& ev, FILE* pFile);
bool PassTrigger      (const fwlite::ChainEvent& ev);

void SetWeight(double IntegratedLuminosityInPb=-1, double CrossSection=0, double MCEvents=0);


/////////////////////////// VARIABLE DECLARATION /////////////////////////////


TH1D*  Data_Pt[40*6];
TH1D*  Data_I [40*6];

TH1D*  Data_Mass[40*6];
TH1D*** Sign_Mass;
TH1D*  Pred_Mass[40*6];
TH1D*  MCTr_Mass[40*6];

double N_A[40*6];	double N_Aerr[40*6];
double N_B[40*6];	double N_Berr[40*6];
double N_C[40*6];	double N_Cerr[40*6];
double N_D[40*6];	double N_Derr[40*6];

TH1D*  Pred_P   [40*6];
TH1D*  Pred_I   [40*6];
TH2D*  Data_PI_A[40*6];
TH2D*  Data_PI_B[40*6];
TH2D*  Data_PI_C[40*6];
TH2D*  Data_PI_D[40*6];
TH2D*  Pred_PI  [40*6];

TH1D* Pred_Expected_Entries;
TH1D* Pred_Observed_Entries;
TH1D* Pred_Correlation_A;
TH1D* Pred_Correlation_B;
TH1D* Pred_Correlation_C;
TH1D* Pred_Correlation_D;

stPlots DataPlots;
std::vector<stPlots> SignPlots;
stPlots MCTrPlots;

std::vector<stSignal> signals;

/////////////////////////// CODE PARAMETERS /////////////////////////////

std::vector<string> DataFileName;
std::vector<string> MCTrFileName;
string TreeName;

float Event_Weight = 1;

void Analysis_Step2345(string MODE="COMPILE", double WP_Pt=-1.0, double WP_I=-1, int SplitMode_=2, int dEdxSel_=0, int dEdxMass_=0, int TypeMode_=0)
{
   if(MODE=="COMPILE")return;

   //////////////////////////////////////////////////     GLOBAL INIT
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.13);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);

   GetSignalDefinition(signals);

   char Buffer[2048];   
   char Command[2048];
   DataFileName.clear();
   GetInputFiles(DataFileName, "Data");
//   DataFileName.push_back("InputFiles/Data.root");

   MCTrFileName.clear();
//   MCTrFileName.push_back("InputFiles/MC.root");
   GetInputFiles(MCTrFileName, "MC");


   dEdxSeleIndex = dEdxSel_;
   dEdxMassIndex = dEdxMass_;

   TreeName = "HSCPTreeBuilder/MyTree";

   TypeMode  = TypeMode_;
   SplitMode = SplitMode_;
   if(SplitMode>0){    GlobalMinHit = 1;
   }else{         GlobalMinHit = 9;   }

   DefaultCutPt   = 0;
   DefaultCutI    = 0;
   SelectionCutPt = pow(10,WP_Pt);
   SelectionCutI  = pow(10,WP_I);

   sprintf(Buffer,"Results/"       );                                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%s%s/"         ,Buffer,MODE.c_str());               sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sSplitMode%i/",Buffer,SplitMode);                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sMinHit%02i/" ,Buffer,GlobalMinHit);               sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sSele_%s/"    ,Buffer,dEdxLabel[dEdxSeleIndex]);   sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sMass_%s/"    ,Buffer,dEdxLabel[dEdxMassIndex]);   sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sType%i/"     ,Buffer,TypeMode);                   sprintf(Command,"mkdir %s",Buffer); system(Command);

   printf("Running Mode = %s\n",MODE.c_str());  
   if(MODE==string("MAKE_MAP")){
      Find_WorkingPoint(Buffer);
   }else if(MODE==string("MERGE_MAP")){
      Merge_Map(Buffer);
   }else if(MODE==string("ANALYSE")){
      sprintf(Buffer,"%sWPPt%+03i/"  ,Buffer,(int)(10*log10(SelectionCutPt)));   sprintf(Command,"mkdir %s",Buffer); system(Command);
      sprintf(Buffer,"%sWPI%+03i/"   ,Buffer,(int)(10*log10(SelectionCutI)));    sprintf(Command,"mkdir %s",Buffer); system(Command);
      CompleteAnalysis(Buffer);
   }else{
      printf("UNKNOWN MODE\n");
   }
}


void CompleteAnalysis(char* SavePath)
{
   char Buffer[2048];
   sprintf(Buffer,"%s/DumpHistos.root",SavePath);
   TFile* tmp = new TFile(Buffer,"RECREATE");
   TH1::AddDirectory(kTRUE);

   InitHistos();
   Analysis_Step2();
   Analysis_Step3();
   Analysis_Step4(SavePath);
   Analysis_Step5(SavePath);
   tmp->Write();
   tmp->Close();
//   CleanUpHistos();
}

void Find_WorkingPoint(char* SavePath)
{
   InitHistos();
   Analysis_Step2();
   Analysis_Step4(NULL);

   char Buffer[2048];
   sprintf(Buffer,"%s/WP_Pt%+03i_I%+03i.tmp",SavePath,(int)(10*log10(SelectionCutPt)),(int)(10*log10(SelectionCutI)) );
   printf("Buffer = %s\n",Buffer);
   FILE* pFile = fopen(Buffer,"w");
   for(unsigned int s=0;s<signals.size();s++){
      double BEff = DataPlots   .WN_I/DataPlots   .WN_Total; //if(BEff<1E-10)BEff=1E-10;
      double SEff = SignPlots[s].WN_I/SignPlots[s].WN_Total; //if(SEff<1E-10)SEff=1E-10;

      fprintf(pFile,"Signal=%10s WP=(%f,%f) --> Numbers: B=%3.2E S=%3.2E S/B=%4.3E EFFIENCIES: B=%3.2E S=%3.2E S/B=%4.3E\n",signals[s].Name.c_str(),log10(SelectionCutPt),log10(SelectionCutI), DataPlots.WN_I,SignPlots[s].WN_I,SignPlots[s].WN_I/DataPlots.WN_I, BEff,SEff,SEff/BEff);
      fflush(pFile);
   }
   fclose(pFile); 
   CleanUpHistos();
}

void Merge_Map(char* SavePath)
{
   if(!SavePath)return;
   gStyle->SetNdivisions(520,"XYZ");

   string MapFilePath = SavePath;
   MapFilePath.replace(MapFilePath.find("MERGE_MAP"),9,"MAKE_MAP");

   TH2D** WP_B     = new TH2D*[signals.size()];
   TH2D** WP_BEff  = new TH2D*[signals.size()];
   TH2D** WP_S     = new TH2D*[signals.size()];
   TH2D** WP_SEff  = new TH2D*[signals.size()];
   TH2D** WP_SB    = new TH2D*[signals.size()];
   TH2D** WP_SBEff = new TH2D*[signals.size()];
   for(unsigned int s=0;s<signals.size();s++){
      string Name;
      Name = string("WP B ") + signals[s].Name;
      WP_B    [s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
      Name = string("WP B ") + signals[s].Name + "Eff";
      WP_BEff [s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
      Name  = string("WP S ") + signals[s].Name;
      WP_S    [s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
      Name = string("WP S ") + signals[s].Name + "Eff";
      WP_SEff [s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
      Name = string("WP S/B ") + signals[s].Name ;
      WP_SB   [s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
      Name = string("WP S/B ") + signals[s].Name + "Eff";
      WP_SBEff[s] = new TH2D(Name.c_str(), Name.c_str()    , 11,-5.25,0.25,11,-5.25,0.25);
   }

   char Buffer[2048];
   sprintf(Buffer,"%s/Map.txt",SavePath);
   FILE* pDump = fopen(Buffer,"w");

   for(float WP_Pt=0;WP_Pt>=-5;WP_Pt-=0.5f){
   for(float WP_I =0;WP_I >=-5;WP_I -=0.5f){

      sprintf(Buffer,"%s/WP_Pt%+03i_I%+03i.tmp",MapFilePath.c_str(),(int)(10*WP_Pt),(int)(10*WP_I));
      FILE* pFile = fopen(Buffer,"r");
      if(!pFile){
         printf("Not Found: %s\n",Buffer);
         continue;
      }

      for(unsigned int S=0;S<signals.size();S++){
         char  sname[256];
         float wp1, wp2, b,s,sb,be,se,sebe;

         fscanf (pFile,"Signal=%s WP=(%f,%f) --> Numbers: B=%E S=%E S/B=%E EFFIENCIES: B=%E S=%E S/B=%E\n",sname,&wp1,&wp2,&b,&s,&sb,&be,&se,&sebe);
         if(b==0)sb=0;      if(be==0)sebe=0;
         fprintf(pDump,"Signal=%10s WP=(%+6.2f,%+6.2f) --> Numbers: B=%3.2E S=%3.2E S/B=%4.3E EFFIENCIES: B=%3.2E S=%3.2E S/B=%4.3E\n",sname,wp1,wp2,b,s,sb,be,se,sebe);

         int Bin_Pt = WP_B[S]->GetXaxis()->FindBin(WP_Pt);
         int Bin_I  = WP_B[S]->GetYaxis()->FindBin(WP_I );

         WP_BEff[S]->SetBinContent(Bin_Pt,Bin_I,be);
         if(!(se!=se))WP_SEff[S]->SetBinContent(Bin_Pt,Bin_I,se);
         if(!(se!=se) && be!=0)WP_SBEff[S]->SetBinContent(Bin_Pt,Bin_I,sebe);
      
         WP_B[S]->SetBinContent(Bin_Pt,Bin_I,b);
         if(!(s!=s))WP_S[S]->SetBinContent(Bin_Pt,Bin_I,s);
         if(!(s!=s) && b>0)WP_SB[S]->SetBinContent(Bin_Pt,Bin_I,sb);
      }
      fclose(pFile);
   }}
   fclose(pDump);
   
   TCanvas* c1;
   for(unsigned int s=0;s<signals.size();s++){
      c1  = new TCanvas("BEFF", "BEFF", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      gStyle->SetPaintTextFormat("1.1E");
      WP_BEff[s]->SetMarkerSize(1.0);
      WP_BEff[s]->SetTitle("");
      WP_BEff[s]->SetStats(kFALSE);
      WP_BEff[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_BEff[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_BEff[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_BEff[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_BEff");
      delete c1;

      c1  = new TCanvas("SEFF", "SEFF", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      WP_SEff[s]->SetMarkerSize(1.0);
      WP_SEff[s]->SetTitle("");
      WP_SEff[s]->SetStats(kFALSE);
      WP_SEff[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_SEff[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_SEff[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_SEff[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_SEff");
      delete c1;

      c1  = new TCanvas("SBEFF", "SBEFF", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      WP_SBEff[s]->SetMarkerSize(1.0);
      WP_SBEff[s]->SetTitle("");
      WP_SBEff[s]->SetStats(kFALSE);
      WP_SBEff[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_SBEff[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_SBEff[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_SBEff[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_SBEff");
      delete c1;

      c1  = new TCanvas("B", "B", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      gStyle->SetPaintTextFormat("1.1E");
      WP_B[s]->SetMarkerSize(1.0);
      WP_B[s]->SetTitle("");
      WP_B[s]->SetStats(kFALSE);
      WP_B[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_B[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_B[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_B[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_B");

      WP_B[s]->SetAxisRange(1,1E6,"Z");
      WP_B[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("MapRanged_") + signals[s].Name + "_B");
      delete c1;

      c1  = new TCanvas("S", "S", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      WP_S[s]->SetMarkerSize(1.0);
      WP_S[s]->SetTitle("");
      WP_S[s]->SetStats(kFALSE);
      WP_S[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_S[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_S[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_S[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_S");
 
      WP_B[s]->SetAxisRange(1E-4,1E2,"Z");
      WP_S[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("MapRanged_") + signals[s].Name + "_S");
      delete c1;

      c1  = new TCanvas("SB", "SB", 600,600);
      c1->SetLogx(false);
      c1->SetLogy(false);
      c1->SetLogz(true);
      c1->SetGridx(true);
      c1->SetGridy(true);
      WP_SB[s]->SetMarkerSize(1.0);
      WP_SB[s]->SetTitle("");
      WP_SB[s]->SetStats(kFALSE);
      WP_SB[s]->GetXaxis()->SetTitle("Selection Efficiency on PT (log10)");
      WP_SB[s]->GetYaxis()->SetTitle("Selection Efficiency on I  (log10)");
      WP_SB[s]->GetYaxis()->SetTitleOffset(1.60);
      WP_SB[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("Map_") + signals[s].Name + "_SB");

      WP_SB[s]->SetAxisRange(1E-6,1E1,"Z");
      WP_SB[s]->Draw("COLZ TEXT45");
      SaveCanvas(c1, SavePath, string("MapRanged_") + signals[s].Name + "_SB" );
      delete c1;
   }
}


bool isGoodCandidate(const susybsm::HSCParticle& hscp, const reco::Vertex& vertex, double PtCut, double ICut, stPlots* st)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::matchedStandAloneMuon || hscp.type() == HSCParticleType::globalMuon))return false;

   reco::TrackRef   trackRef = hscp.trackRef(); if(trackRef.isNull())return false;
   reco::Track      track    = *trackRef;

   if(st){st->WN_Total+=Event_Weight;	st->UN_Total++;}

   double dz  = track.dz (vertex.position());
   double dxy = track.dxy(vertex.position());

   if(st){st->BS_Hits->Fill(track.found(),Event_Weight);}
   if(track.found()<GlobalMinHit)return false;
   if(st){st->AS_Hits->Fill(track.found(),Event_Weight);}
   if(st){st->WN_Hits  +=Event_Weight;   st->UN_Hits++;}

   if(st){st->BS_Qual->Fill(track.qualityMask(),Event_Weight);}
   if(track.qualityMask()<GlobalMinQual )return false;
   if(st){st->AS_Qual->Fill(track.qualityMask(),Event_Weight);}
   if(st){st->WN_Qual  +=Event_Weight;   st->UN_Qual++;}

   if(st){st->BS_Chi2->Fill(track.chi2()/track.ndof(),Event_Weight);}
   if(track.chi2()/track.ndof()>GlobalMaxChi2 )return false;
   if(st){st->AS_Chi2->Fill(track.chi2()/track.ndof(),Event_Weight);}
   if(st){st->WN_Chi2  +=Event_Weight;   st->UN_Chi2++;}

   if(st){st->BS_Pterr ->Fill(track.ptError()/track.pt(),Event_Weight);}
   if((track.ptError()/track.pt())>GlobalMaxPterr)return false;
   if(st){st->AS_Pterr ->Fill(track.ptError()/track.pt(),Event_Weight);}
   if(st){st->WN_Pterr   +=Event_Weight;   st->UN_Pterr ++;}

   if(st){st->BS_DZ->Fill(fabs(dz),Event_Weight);}
   if(fabs(dz)>GlobalMaxDZ )return false; 
   if(st){st->AS_DZ->Fill(fabs(dz),Event_Weight);} 
   if(st){st->WN_DZ   +=Event_Weight;   st->UN_DZ++;}

   if(st){st->BS_DXY->Fill(fabs(dxy),Event_Weight);}  
   if(fabs(dxy)>GlobalMaxDXY )return false; 
   if(st){st->AS_DXY->Fill(fabs(dxy),Event_Weight);}
   if(st){st->WN_DXY  +=Event_Weight;   st->UN_DXY++;}

   if(st){st->BS_MPt ->Fill(track.pt(),Event_Weight);}
   if(track.pt()<GlobalMinPt)return false;
   if(st){st->AS_MPt ->Fill(track.pt(),Event_Weight);}
   if(st){st->WN_MPt   +=Event_Weight;   st->UN_MPt ++;}

   if(st){st->BS_MI ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(hscp.dedx(dEdxSeleIndex).dEdx()<GlobalMinI)return false;
   if(st){st->AS_MI ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->WN_MI   +=Event_Weight;   st->UN_MI++;}

   if(st){st->BS_EtaP ->Fill(track.eta(),track.p(),Event_Weight);}
   if(st){st->BS_EtaPt->Fill(track.eta(),track.pt(),Event_Weight);}
   if(st){st->BS_PI   ->Fill(track.p(),hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}

   if(st){st->BS_Pt  ->Fill(track.pt(),Event_Weight);}
   if(track.pt()<PtCut)return false;
   if(st){st->AS_Pt  ->Fill(track.pt(),Event_Weight);}
   if(st){st->WN_Pt    +=Event_Weight;   st->UN_Pt ++;}

   if(st){st->BS_I  ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(hscp.dedx(dEdxSeleIndex).dEdx()<ICut)return false;
   if(st){st->AS_I  ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->WN_I    +=Event_Weight;   st->UN_I++;}

   if(st){st->AS_EtaP ->Fill(track.eta(),track.p(),Event_Weight);}
   if(st){st->AS_EtaPt->Fill(track.eta(),track.pt(),Event_Weight);}
   if(st){st->AS_PI   ->Fill(track.p(),hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}

   return true;
}

void DumpCandidateInfo(const susybsm::HSCParticle& hscp, const reco::Vertex& vertex, const fwlite::ChainEvent& ev, FILE* pFile)
{
   reco::MuonRef  muon  = hscp.muonRef();
   reco::TrackRef track = hscp.trackRef();
   if(track.isNull())return;

   double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
   if(Mass<MinCandidateMass || Mass>10000)return;
   double dz  = track->dz (vertex.position());
   double dxy = track->dxy(vertex.position());


   int HitIndex,EtaIndex;
   GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
   int CutIndex = GetCutIndex(HitIndex,EtaIndex);

   fprintf(pFile,"\n");
   fprintf(pFile,"---------------------------------------------------------------------------------------------------\n");
   fprintf(pFile,"Candidate Type = %i --> Mass: %7.2f GeV\n",hscp.type(),Mass);
   fprintf(pFile,"------------------------------------------ EVENT INFO ---------------------------------------------\n");
   fprintf(pFile,"Run=%i Lumi=%i Event=%i BX=%i  Orbit=%i Store=%i\n",ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().event(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().orbitNumber(),ev.eventAuxiliary().storeNumber());
   fprintf(pFile,"------------------------------------------ INNER TRACKER ------------------------------------------\n");
   fprintf(pFile,"Quality = %i Chi2/NDF=%6.2f dz=+%6.2f dxy=%+6.2f charge:%+i\n",track->qualityMask(), track->chi2()/track->ndof(), dz, dxy, track->charge());
   fprintf(pFile,"P=%7.2f  Pt=%7.2f+-%6.2f (Cut=%6.2f) Eta=%+6.2f  Phi=%+6.2f  NOH=%2i\n",track->p(),track->pt(), track->ptError(), CutPt[CutIndex], track->eta(), track->phi(), track->found() );
   fprintf(pFile,"dEdx for selection:%6.2f (Cut=%6.2f) NOM %2i NOS %2i\n",hscp.dedx(dEdxSeleIndex).dEdx(),CutI[CutIndex],hscp.dedx(dEdxSeleIndex).numberOfMeasurements(),hscp.dedx(dEdxSeleIndex).numberOfSaturatedMeasurements());
   fprintf(pFile,"dEdx for mass reco:%6.2f             NOM %2i NOS %2i\n",hscp.dedx(dEdxMassIndex).dEdx(),hscp.dedx(dEdxMassIndex).numberOfMeasurements(),hscp.dedx(dEdxMassIndex).numberOfSaturatedMeasurements());
   fprintf(pFile,"------------------------------------------ MUON INFO ----------------------------------------------\n");
   if(!muon.isNull()){
	   fprintf(pFile,"Quality=%i type=%i P=%7.2f  Pt=%7.2f Eta=%+6.2f Phi=%+6.2f #Chambers=%i\n",muon->isQualityValid(),muon->type(),muon->p(),muon->pt(),muon->eta(),muon->phi(),muon->numberOfChambers());
   fprintf(pFile,"muonTimeCombined: NDOF=%2i InvBeta=%6.2f+-%6.2f FreeInvBeta=%6.2f+-%6.2f\n",hscp.muonTimeCombined().nDof(),hscp.muonTimeCombined().inverseBeta(),hscp.muonTimeCombined().inverseBetaErr(),hscp.muonTimeCombined().freeInverseBeta(),hscp.muonTimeCombined().freeInverseBetaErr());
   fprintf(pFile,"muonTimeDT      : NDOF=%2i InvBeta=%6.2f+-%6.2f FreeInvBeta=%6.2f+-%6.2f\n",hscp.muonTimeDt().nDof(),hscp.muonTimeDt().inverseBeta(),hscp.muonTimeDt().inverseBetaErr(),hscp.muonTimeDt().freeInverseBeta(),hscp.muonTimeDt().freeInverseBetaErr());
   fprintf(pFile,"muonTimeCSC     : NDOF=%2i InvBeta=%6.2f+-%6.2f FreeInvBeta=%6.2f+-%6.2f\n",hscp.muonTimeCsc().nDof(),hscp.muonTimeCsc().inverseBeta(),hscp.muonTimeCsc().inverseBetaErr(),hscp.muonTimeCsc().freeInverseBeta(),hscp.muonTimeCsc().freeInverseBetaErr());
   }
   fprintf(pFile,"------------------------------------------ RPC INFO -----------------------------------------------\n");
   if(hscp.hasRpcInfo()){
   fprintf(pFile,"isCandidate %i Beta=%6.2f #Hits=%i\n",hscp.rpc().isCandidate,hscp.rpc().beta,hscp.rpc().hits.size());
   }
   fprintf(pFile,"------------------------------------------ CALO INFO ----------------------------------------------\n");
   if(hscp.hasCaloInfo()){
   fprintf(pFile,"HCAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f HO E=%6.2f\n",hscp.calo().hcalenergy,hscp.calo().hcal3by3dir, hscp.calo().hcal5by5dir, hscp.calo().hoenergy);
   fprintf(pFile,"ECAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f\n"           ,hscp.calo().ecalenergy,hscp.calo().ecal3by3dir, hscp.calo().ecal5by5dir);
   fprintf(pFile,"ECAL: time=%6.2f beta=%6.2f trkisodr=%6.2f\n"    ,hscp.calo().ecaltime  ,hscp.calo().ecalbeta   , hscp.calo().trkisodr);
   }
   fprintf(pFile,"---------------------------------------------------------------------------------------------------\n");
   fprintf(pFile,"\n");
}

bool PassTrigger(const fwlite::ChainEvent& ev)
{
      //Some run on data have a different trigger table... (they are 261Events on 871561)
      if(ev.eventAuxiliary().run()>=132440 && ev.eventAuxiliary().run()<=132528)return false;

      edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
      if(!tr.isValid())return false;
//      for(unsigned int i=0;i<tr.size();i++){
//         printf("Path %3i %50s --> %1i\n",i, tr.triggerName(i).c_str(),tr.accept(i));
//      }fflush(stdout);

      bool JetMetSD = false;
      if(TypeMode!=1){
      JetMetSD |= tr.accept("HLT_Jet15U");
      JetMetSD |= tr.accept("HLT_DiJetAve15U_8E29");
      JetMetSD |= tr.accept("HLT_FwdJet20U");
      JetMetSD |= tr.accept("HLT_Jet30U");
      JetMetSD |= tr.accept("HLT_Jet50U");
      JetMetSD |= tr.accept("HLT_DiJetAve30U_8E29");
      JetMetSD |= tr.accept("HLT_QuadJet15U");
      JetMetSD |= tr.accept("HLT_MET45");
      JetMetSD |= tr.accept("HLT_MET100");
      JetMetSD |= tr.accept("HLT_HT100U");
      JetMetSD |= tr.accept("HLT_SingleLooseIsoTau20");
      JetMetSD |= tr.accept("HLT_DoubleLooseIsoTau15");
      JetMetSD |= tr.accept("HLT_DoubleJet15U_ForwardBackward");
      JetMetSD |= tr.accept("HLT_BTagMu_Jet10U");
      JetMetSD |= tr.accept("HLT_BTagIP_Jet50U");
      JetMetSD |= tr.accept("HLT_StoppedHSCP_8E29");
//      printf("JetMetSD=%i\n",JetMetSD);
      }

      bool MuSD = false;
      MuSD |= tr.accept("HLT_L2Mu0");
      MuSD |= tr.accept("HLT_L2Mu3");
//    MuSD |= tr.accept("HLT_L2Mu5");
      MuSD |= tr.accept("HLT_L1Mu20");
      MuSD |= tr.accept("HLT_L2Mu9");
      MuSD |= tr.accept("HLT_L2Mu11");
      MuSD |= tr.accept("HLT_L1Mu14_L1SingleEG10");
      MuSD |= tr.accept("HLT_L1Mu14_L1SingleJet6U");
      MuSD |= tr.accept("HLT_L1Mu14_L1ETM30");
      MuSD |= tr.accept("HLT_L2DoubleMu0");
      MuSD |= tr.accept("HLT_L1DoubleMuOpen");
      MuSD |= tr.accept("HLT_DoubleMu0");
      MuSD |= tr.accept("HLT_DoubleMu3");
      MuSD |= tr.accept("HLT_Mu3");
      MuSD |= tr.accept("HLT_Mu5");
      MuSD |= tr.accept("HLT_Mu9");
      MuSD |= tr.accept("HLT_IsoMu3");
      MuSD |= tr.accept("HLT_Mu0_L1MuOpen");
      MuSD |= tr.accept("HLT_Mu0_Track0_Jpsi");
      MuSD |= tr.accept("HLT_Mu3_L1MuOpen");
      MuSD |= tr.accept("HLT_Mu3_Track0_Jpsi");
      MuSD |= tr.accept("HLT_Mu5_L1MuOpen");
      MuSD |= tr.accept("HLT_Mu5_Track0_Jpsi");
      MuSD |= tr.accept("HLT_Mu0_L2Mu0");
      MuSD |= tr.accept("HLT_Mu3_L2Mu0");
      MuSD |= tr.accept("HLT_Mu5_L2Mu0");

      return JetMetSD | MuSD;
}

void Analysis_Step2()
{
   printf("Step2: Optimizing Cuts\n");

   fwlite::ChainEvent tree(DataFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Finding Cuts                 :");
   int TreeStep = tree.size()/50;if(TreeStep==0)TreeStep=1;

   for(Long64_t ientry=0;ientry<tree.size();ientry++){
      tree.to(ientry);
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
      if(!PassTrigger(tree) )continue;

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(tree,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      std::vector<reco::Vertex> vertexColl = *vertexCollHandle;
      if(vertexColl.size()<1){printf("NO VERTEX\n"); continue;}

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(tree,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(!isGoodCandidate(hscp,vertexColl[0]))continue;

         int HitIndex, EtaIndex;
         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);

         FillHisto(HitIndex, EtaIndex, Data_Pt, track->pt()                     ,Event_Weight);
         FillHisto(HitIndex, EtaIndex, Data_I , hscp.dedx(dEdxSeleIndex).dEdx() ,Event_Weight);
      } // end of Track Loop
   }// end of Event Loop
   printf("\n");

   if(SelectionCutPt>=0){
      for(unsigned int i=0;i<40*6;i++){
         if(SplitMode==0 && i>0)continue;
         if(SplitMode==1 && (i==0 || i%6!=0))continue;
         if(SplitMode==2 && (i< 6 || i%6==0))continue;
         CutPt[i] = CutFromEfficiency(Data_Pt[i],SelectionCutPt);
         if(CutPt[i]>Data_Pt[i]->GetXaxis()->GetXmax())CutPt[i]=99999;
      }
   }else{
      for(unsigned int i=0;i<6*40;i++){ CutPt[i]=DefaultCutPt; }
   }  

   if(SelectionCutI>=0){
      for(unsigned int i=0;i<40*6;i++){
         if(SplitMode==0 && i>0)continue;
         if(SplitMode==1 && (i==0 || i%6!=0))continue;
         if(SplitMode==2 && (i< 6 || i%6==0))continue;
         CutI[i] = CutFromEfficiency(Data_I[i],SelectionCutI);
//         if(CutI[i]>Data_I[i]->GetXaxis()->GetXmax()){
//            printf("Interval %i --> %f entries --> Cut = %f\n",i,Data_I[i]->GetEntries(),CutI[i]);
//         }
         if(CutI[i]>Data_I[i]->GetXaxis()->GetXmax())CutI[i]=99999;
          
      }
   }else{
      for(unsigned int i=0;i<6*40;i++){ CutI[i]=DefaultCutI; }
   }
  
//   delete tree;
}


void Analysis_Step3()
{
   printf("Step3: Making the Background Prediction\n");

   //////////////////////////////////////////////////     TREE LOOP
   fwlite::ChainEvent tree(DataFileName);
   SetWeight(-1);
   printf("Predicting (Looping on Tree) :");
   int TreeStep = tree.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t ientry=0;ientry<tree.size();ientry++){
      tree.to(ientry);
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
      if(!PassTrigger(tree) )continue;

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(tree,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      std::vector<reco::Vertex> vertexColl = *vertexCollHandle;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(tree,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(!isGoodCandidate(hscp,vertexColl[0]))continue;

         int HitIndex, EtaIndex;
         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
         
         bool PassMinPt = GlobalMinPt;
         bool PassPtCut = track->pt()>=CutPt[CutIndex];
         bool PassMinI  = (hscp.dedx(dEdxSeleIndex).dEdx()>=GlobalMinI);
         bool PassICut  = (hscp.dedx(dEdxSeleIndex).dEdx()>=CutI[CutIndex]);

          //printf("Pt %6.2f|%6.2f  I %6.2f|%6.2f\n", track->pt(),CutPt[CutIndex], hscp.dedx(dEdxSeleIndex).dEdx(), CutI[CutIndex]);

          ///\ I
          // -----------------------------
          // |   |           |             |
          // |   |           |             |
          // |   |    B      |     D       |
          // |   |           |             |
          // ------------------------------
          // |   |           |             |
          // |   |   A       |    C        |
          // |3.5|           |             |
          // |---|-----------|-------------|
          // |   |           |             |
          // 0---10---------------------------> PT

          if(PassMinPt && PassMinI && !PassPtCut && !PassICut){        //Region A
             FillArray(HitIndex, EtaIndex, N_A             , Event_Weight);
             FillArray(HitIndex, EtaIndex, N_Aerr          , Event_Weight*Event_Weight);
             FillHisto(HitIndex, EtaIndex, Data_PI_A, track->p(),hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
          }else if(PassMinPt && PassMinI && !PassPtCut &&  PassICut){  //Region B
             FillArray(HitIndex, EtaIndex, N_B             , Event_Weight);
             FillArray(HitIndex, EtaIndex, N_Berr          , Event_Weight*Event_Weight);
             FillHisto(HitIndex, EtaIndex, Data_PI_B, track->p(),hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
             FillHisto(HitIndex, EtaIndex, Pred_I, hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
          }else if(PassMinPt && PassMinI && PassPtCut && !PassICut){   //Region C
             FillArray(HitIndex, EtaIndex, N_C             , Event_Weight);
             FillArray(HitIndex, EtaIndex, N_Cerr          , Event_Weight*Event_Weight);
             FillHisto(HitIndex, EtaIndex, Data_PI_C, track->p(),hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
             FillHisto(HitIndex, EtaIndex, Pred_P, track->p(), Event_Weight);
          }else if(PassMinPt && PassMinI && PassPtCut && PassICut){    //Region D
             FillArray(HitIndex, EtaIndex, N_D             , Event_Weight);
             FillArray(HitIndex, EtaIndex, N_Derr          , Event_Weight*Event_Weight);
             FillHisto(HitIndex, EtaIndex, Data_PI_D, track->p(),hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
          }        
      } // end of Track Loop
   }// end of Event Loop
   printf("\n");

   //////////////////////////////////////////////////      MAKING THE PREDICTION

   printf("Predicting (Finding Prob)    :");
   TreeStep = (40*6)/50;if(TreeStep==0)TreeStep=1;
   int CountStep = 0;
   for(unsigned int i=0;i<40*6;i++){
      if(i%TreeStep==0 && CountStep<=50){printf(".");fflush(stdout);CountStep++;}

      if(SplitMode==0 && i>0)continue;
      if(SplitMode==1 && (i==0 || i%6!=0))continue;
      if(SplitMode==2 && (i< 6 || i%6==0))continue;

      Pred_Correlation_A->SetBinContent(i, Data_PI_A[i]->GetCorrelationFactor()  );    Pred_Correlation_A->SetBinError(i,0);
      Pred_Correlation_B->SetBinContent(i, Data_PI_B[i]->GetCorrelationFactor()  );    Pred_Correlation_B->SetBinError(i,0);
      Pred_Correlation_C->SetBinContent(i, Data_PI_C[i]->GetCorrelationFactor()  );    Pred_Correlation_C->SetBinError(i,0);
      Pred_Correlation_D->SetBinContent(i, Data_PI_D[i]->GetCorrelationFactor()  );    Pred_Correlation_D->SetBinError(i,0);

      if(N_A[i]>0){
         Pred_Expected_Entries->SetBinContent(i, ((N_C[i]*N_B[i])/N_A[i])  );
         Pred_Expected_Entries->SetBinError  (i, sqrt((pow(N_B[i]/N_A[i],2)*N_Cerr[i]) + (pow(N_C[i]/N_A[i],2)*N_Berr[i]) + (pow((N_B[i]*(N_C[i])/(N_A[i]*N_A[i])),2)*N_Aerr[i])) );
	 Pred_Observed_Entries->SetBinContent(i, N_D[i]  );
         Pred_Observed_Entries->SetBinError  (i, sqrt(N_Derr[i])  );
      }

      double IntegralP = Pred_P[i]->Integral(0, Pred_P[i]->GetNbinsX()+1);
      double IntegralI = Pred_I[i]->Integral(0, Pred_I[i]->GetNbinsX()+1);
      if(IntegralP>0)Pred_P[i]->Scale(1.0/IntegralP);
      if(IntegralI>0)Pred_I[i]->Scale(1.0/IntegralI);

      double N_A_L = N_A[i] - sqrt(N_Aerr[i]); if(N_A_L<0)N_A_L=0;
      double N_A_C = N_A[i];                   if(N_A_C<0)N_A_C=0;
      double N_A_U = N_A[i] + sqrt(N_Aerr[i]); if(N_A_U<0)N_A_U=0;

      double N_B_L = N_B[i] - sqrt(N_Berr[i]); if(N_B_L<0)N_B_L=0;
      double N_B_C = N_B[i];                   if(N_B_C<0)N_B_C=0;
      double N_B_U = N_B[i] + sqrt(N_Berr[i]); if(N_B_U<0)N_B_U=0;

      double N_C_L = N_C[i] - sqrt(N_Cerr[i]); if(N_C_L<0)N_C_L=0;
      double N_C_C = N_C[i];                   if(N_C_C<0)N_C_C=0;
      double N_C_U = N_C[i] + sqrt(N_Cerr[i]); if(N_C_U<0)N_C_U=0;

      double NExpectedBckgEntriesC;
      double NExpectedBckgEntriesC2;

      if(AbsolutePredictiction){
         if(N_A[i]>0){
            NExpectedBckgEntriesC  = ((N_C[i]*N_B[i])/N_A[i]);
            NExpectedBckgEntriesC2 = sqrt((pow(N_B[i]/N_A[i],2)*N_Cerr[i]) + (pow(N_C[i]/N_A[i],2)*N_Berr[i]) + (pow((N_B[i]*(N_C[i])/(N_A[i]*N_A[i])),2)*N_Aerr[i]));
         }else{
            NExpectedBckgEntriesC  = 0;
            NExpectedBckgEntriesC2 = 0;
         }
      }else{
         NExpectedBckgEntriesC  = N_A[i]/N_A[0];
         NExpectedBckgEntriesC2 = NExpectedBckgEntriesC*NExpectedBckgEntriesC;
      }

      //Loop on Mass Line
      for(int m=0;m<Pred_Mass[i]->GetNbinsX()+1;m++){
         //Find which bins contributes to this particular mass bin
         std::vector<std::pair<int,int> > BinThatGivesThisMass;
         for(int x=1;x<Pred_P[i]->GetNbinsX()+1;x++){
         for(int y=1;y<Pred_I[i]->GetNbinsX()+1;y++){
            double Mass = GetMass( Pred_P[i]->GetXaxis()->GetBinCenter(x) , Pred_I[i]->GetXaxis()->GetBinCenter(y) );
            if(Mass>Pred_Mass[i]->GetXaxis()->GetBinLowEdge(m) && Mass<Pred_Mass[i]->GetXaxis()->GetBinUpEdge(m)){
               BinThatGivesThisMass.push_back(std::make_pair(x,y));
            }
         }}

         double MBinContent=0;
         double MBinError2 =0;

         //Loops on the bins that contribute to this mass bin.
         for(unsigned int b1=0;b1<BinThatGivesThisMass.size();b1++){
            double bx1 = BinThatGivesThisMass[b1].first;
            double by1 = BinThatGivesThisMass[b1].second;
            double vx1 = Pred_P[i]->GetBinContent(bx1);
            double vy1 = Pred_I[i]->GetBinContent(by1);
            double ex1 = Pred_P[i]->GetBinError(bx1);
            double ey1 = Pred_I[i]->GetBinError(by1);
            double vz1 = vx1*vy1;


            MBinContent += vz1*NExpectedBckgEntriesC;
            Pred_PI[i]->SetBinContent(bx1,by1,vz1*NExpectedBckgEntriesC);

            //Compute the errors with a covariance matrix (on the fly) --> Only vertical and horizontal term contributes.
            for(unsigned int b2=0;b2<BinThatGivesThisMass.size();b2++){
               double bx2 = BinThatGivesThisMass[b2].first;
               double by2 = BinThatGivesThisMass[b2].second;
               double vx2 = Pred_P[i]->GetBinContent(bx2);
               double vy2 = Pred_I[i]->GetBinContent(by2);
               //double ex2 = Pred_P[i]->GetBinError(bx2);
               //double ey2 = Pred_I[i]->GetBinError(by2);


               if(bx1==bx2 && by1==by2){
                  //Correlation with itself!
                  MBinError2 += NExpectedBckgEntriesC2*(ex1*ex1+ey1*ey1)*vz1*vz1;
               }else if(by1==by2){
                  //Vertical term
                  MBinError2 += NExpectedBckgEntriesC2*vx1*vx2*ey1;
               }else if(bx1==bx2){
                  //Horizontal term
                  MBinError2 += NExpectedBckgEntriesC2*vy1*vy2*ex1;
               }else{
                  //Diagonal term... do nothing
               }
            }
//            printf("Interval %i --> M = %i --> %f +- %f, %f+-%f\n",i,m,vx1,ex1,vy1,ey1);
//            printf("Interval %i --> M = %i --> %i Bins concerned --> %fEntries -->  %f +- %f\n",i,m,BinThatGivesThisMass.size(),NExpectedBckgEntriesC,MBinContent,NExpectedBckgEntriesC*sqrt(MBinError));

            Pred_Mass[i]->SetBinContent(m, MBinContent);
            Pred_Mass[i]->SetBinError  (m, sqrt(MBinError2));
         }
         BinThatGivesThisMass.clear();
      }

      if(SplitMode!=0)Pred_PI   [0]->Add(Pred_PI   [i],1);
      if(SplitMode!=0)Pred_Mass [0]->Add(Pred_Mass [i],1);     
   }
   printf("\n");


}

void Analysis_Step4(char* SavePath)
{
   printf("Step4: Building Mass Spectrum for B and S\n");

   int TreeStep;
   int HitIndex, EtaIndex;
   FILE* pFile = NULL;

   //////////////////////////////////////////////////     BUILD BACKGROUND MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_D_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");
   }


   fwlite::ChainEvent treeD(DataFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Building Mass Spectrum for D :");
   TreeStep = treeD.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t ientry=0;ientry<treeD.size();ientry++){
      treeD.to(ientry);
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}

      DataPlots.WN_TotalE+=Event_Weight;       DataPlots.UN_TotalE++;
      if(!PassTrigger(treeD) )continue;
      DataPlots.WN_TotalTE+=Event_Weight;      DataPlots.UN_TotalTE++;

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(treeD,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      std::vector<reco::Vertex> vertexColl = *vertexCollHandle;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeD,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;
 
         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
         if(!isGoodCandidate(hscp,vertexColl[0],CutPt[CutIndex], CutI[CutIndex], &DataPlots))continue;

         //DEBUG
         double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
         double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
         double Mass = GetMass(PBinned,IBinned);

//         double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx());
         FillHisto(HitIndex, EtaIndex,Data_Mass , Mass, Event_Weight);

         if(SavePath)DumpCandidateInfo(hscp, vertexColl[0], treeD, pFile);

      } // end of Track Loop
   }// end of Event Loop
   printf("\n");
   if(pFile){fclose(pFile);pFile=NULL;};

   //////////////////////////////////////////////////     BUILD MCTRUTH MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_M_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");
   }


   fwlite::ChainEvent treeM(MCTrFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Building Mass Spectrum for M :");
   TreeStep = treeM.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t ientry=0;ientry<treeM.size();ientry++){
      treeM.to(ientry);
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}

      MCTrPlots.WN_TotalE+=Event_Weight;       MCTrPlots.UN_TotalE++;
      if(!PassTrigger(treeM) )continue;
      MCTrPlots.WN_TotalTE+=Event_Weight;      MCTrPlots.UN_TotalTE++;

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(treeM,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      std::vector<reco::Vertex> vertexColl = *vertexCollHandle;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeM,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;

         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
         if(!isGoodCandidate(hscp,vertexColl[0],CutPt[CutIndex], CutI[CutIndex], &MCTrPlots))continue;

         //DEBUG
         double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
         double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
         double Mass = GetMass(PBinned,IBinned);

//         double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
         FillHisto(HitIndex, EtaIndex, MCTr_Mass, Mass, Event_Weight);

         if(SavePath)DumpCandidateInfo(hscp, vertexColl[0], treeM, pFile);
      } // end of Track Loop 
   }// end of Event Loop
   printf("\n");
   if(pFile){fclose(pFile);pFile=NULL;};


   //////////////////////////////////////////////////     BUILD SIGNAL MASS SPECTRUM


   for(unsigned int s=0;s<signals.size();s++){
      if(SavePath){
         char Buffer[2048];
         sprintf(Buffer,"%s/Candidate_%s_Dump.txt",SavePath,signals[s].Name.c_str());
         pFile = fopen(Buffer,"w");
      }

      std::vector<string> SignFileName;
      GetInputFiles(SignFileName, signals[s].Name);
//      SignFileName.push_back(string("InputFiles/")+signals[s].Name+".root");

      fwlite::ChainEvent treeS(SignFileName);
      SetWeight(IntegratedLuminosity,signals[s].XSec,(double)treeS.size());
      printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
      printf("Building Mass for %10s :",signals[s].Name.c_str());
      TreeStep = treeS.size()/50;if(TreeStep==0)TreeStep=1;
      for(Long64_t ientry=0;ientry<treeS.size();ientry++){
         treeS.to(ientry);
         if(ientry%TreeStep==0){printf(".");fflush(stdout);}

         SignPlots[s].WN_TotalE+=Event_Weight;       SignPlots[s].UN_TotalE++;
         if(!PassTrigger(treeS) )continue;
         SignPlots[s].WN_TotalTE+=Event_Weight;      SignPlots[s].UN_TotalTE++;

         fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
         vertexCollHandle.getByLabel(treeS,"offlinePrimaryVertices");
         if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
         std::vector<reco::Vertex> vertexColl = *vertexCollHandle;

         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeS,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         genCollHandle.getByLabel(treeS, "genParticles");
         if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
         std::vector<reco::GenParticle> genColl = *genCollHandle;


         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

            int ClosestGen;
            if(DistToHSCP(hscp, genColl, ClosestGen)>0.03)continue;

            GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
            int CutIndex = GetCutIndex(HitIndex,EtaIndex);
            if(!isGoodCandidate(hscp,vertexColl[0],CutPt[CutIndex], CutI[CutIndex], &SignPlots[s]))continue;         

            //DEBUG
            double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
            double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
            double Mass = GetMass(PBinned,IBinned);

//         double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
           FillHisto(HitIndex, EtaIndex, Sign_Mass[s], Mass, Event_Weight);

           if(SavePath)DumpCandidateInfo(hscp, vertexColl[0], treeS, pFile);

         } // end of Track Loop 
       }// end of Event Loop
      printf("\n");
      if(pFile){fclose(pFile);pFile=NULL;};
   }// end of signal Type loop
}

void Analysis_Step5(char* SavePath)
{
   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   printf("Step5: Creating Plots and Dumping information\n");

//   if(!AbsolutePredictiction){
      //SCALE DATA AND PREDICTION
//      Pred_PI   [0]->Scale((Data_Mass[0]->Integral()-Sign_Mass[0]->Integral())/Pred_Mass[0]->Integral());
//      Pred_Mass [0]->Scale((Data_Mass[0]->Integral()-Sign_Mass[0]->Integral())/Pred_Mass[0]->Integral());
//   }

   //////////////////////////////////////////////////     DUMP USEFUL INFORMATION

   for(unsigned int s=0;s<signals.size();s++){
      double BEff = DataPlots   .WN_I/DataPlots   .WN_Total; //if(BEff<1E-10)BEff=1E-10;
      double SEff = SignPlots[s].WN_I/SignPlots[s].WN_Total; //if(SEff<1E-10)SEff=1E-10;
      printf("%10s --> SUMMARY EFFIENCIES: B=%3.2E S=%3.2E >> S/B=%4.3E\n",signals[s].Name.c_str(),BEff,SEff,SEff/BEff);
   }

   char Buffer[2048];
   sprintf(Buffer,"%s/CUT_Dump.txt",SavePath);
   FILE* pFile = fopen(Buffer,"w");
   fprintf(pFile,"MODE          = %i\n",SplitMode);
   fprintf(pFile,"Selection     = %s\n",dEdxLabel[dEdxSeleIndex]);
   fprintf(pFile,"Mass          = %s\n",dEdxLabel[dEdxMassIndex]);
   fprintf(pFile,"WP PT         = %4.3E\n",SelectionCutPt);
   fprintf(pFile,"WP I          = %4.3E\n",SelectionCutI);
   fprintf(pFile,"GlobalMinHit  = %02i\n",GlobalMinHit);
   fprintf(pFile,"GlobalMaxChi2 = %6.2f\n",GlobalMaxChi2);
   fprintf(pFile,"--------------------\n");

   double CutMin_I  = 9999;   double CutMax_I  = 0;
   double CutMin_Pt = 9999;   double CutMax_Pt = 0;
   for(unsigned int i=0;i<40*6;i++){
      if(SplitMode==0 && i>0)continue;
      if(SplitMode==1 && (i==0 || i%6!=0))continue;
      if(SplitMode==2 && (i< 6 || i%6==0))continue;

      if(CutI [i]<CutMin_I                                                         )CutMin_I =CutI [i];
      if(CutI [i]>CutMax_I  && CutI [i]<Data_I [0]->GetXaxis()->GetXmax())CutMax_I =CutI [i];
      if(CutPt[i]<CutMin_Pt                                                        )CutMin_Pt=CutPt[i];
      if(CutPt[i]>CutMax_Pt && CutPt[i]<Data_Pt[0]->GetXaxis()->GetXmax())CutMax_Pt=CutPt[i];
      
      fprintf(pFile,"CutIndex=%03i  PtCut=%14.5f   ICut=%14.5f\n",i,CutPt[i],CutI[i]);
   }

   fprintf(pFile,"--------------------\n");
   fprintf(pFile,"PtCut Range =[%10.5f,%10.5f]\n",CutMin_Pt,CutMax_Pt);
   fprintf(pFile,"ICut  Range =[%10.5f,%10.5f]\n",CutMin_I ,CutMax_I);
   fprintf(pFile,"--------------------\n");


   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"DATA SELECTION DETAILS\n");
   fprintf(pFile,"--------------------\n");
   stPlots_Dump(DataPlots, pFile);

   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"MC TRUTH SELECTION DETAILS\n");
   fprintf(pFile,"--------------------\n");
   stPlots_Dump(MCTrPlots, pFile);

   
   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"SIGNAL SELECTION DETAILS\n");
   fprintf(pFile,"--------------------\n");
   for(unsigned int s=0;s<signals.size();s++){
      fprintf(pFile,"##### ##### %10s ##### #####\n",signals[s].Name.c_str());
      stPlots_Dump(SignPlots[s], pFile);
   }

   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"PREDICTION OF THE MASS DISTRIBUTION\n");
   fprintf(pFile,"--------------------\n");
   for(unsigned int i=0;i<40*6;i++){
      if(SplitMode==0 && i>0)continue;
      if(SplitMode==1 && (i==0 || i%6!=0))continue;
      if(SplitMode==2 && (i< 6 || i%6==0))continue;

      fprintf(pFile,"CutIndex=%03i --> N_A=%E, N_B=%E, N_C=%E, N_D=%E <--> %E\n",i,N_A[i],N_B[i],N_C[i],N_D[i],(N_C[i]*N_B[i])/N_A[i] );      
   }
   fprintf(pFile,"--------------------\n");


//   PrintEventInRange(100,300,pFile);
//   PrintEventInRange(0,10000,pFile);
//   PrintEventInRange(0,20000,pFile);
   fprintf(pFile,"--------------------\n");

   //////////////////////////////////////////////////     CREATE PLOTS OF SELECTION


   if(DataPlots.WN_Total>0){
      sprintf(Buffer,"%s/Selection_Data",SavePath);
      stPlots_Draw(DataPlots, Buffer);
   }

   if(MCTrPlots.WN_Total>0){
      sprintf(Buffer,"%s/Selection_MCTr",SavePath);
      stPlots_Draw(MCTrPlots, Buffer);
   }


   for(unsigned int s=0;s<signals.size();s++){
      if(SignPlots[s].WN_Total>0){
         sprintf(Buffer,"%s/Selection_%s",SavePath,signals[s].Name.c_str());
         stPlots_Draw(SignPlots[s], Buffer);
      }

      sprintf(Buffer,"%s/Selection_Comp_%s",SavePath,signals[s].Name.c_str());
      stPlots_DrawComparison(SignPlots[s], MCTrPlots, DataPlots, Buffer);
   }

   //////////////////////////////////////////////////     CREATE PLOTS OF PREDICTION

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = Pred_Correlation_A;                legend.push_back("Region A");
   Histos[1] = Pred_Correlation_B;                legend.push_back("Region B");
   Histos[2] = Pred_Correlation_C;                legend.push_back("Region C");
   Histos[3] = Pred_Correlation_D;                legend.push_back("Region D");
   DrawSuperposedHistos((TH1**)Histos, legend, "P",  "Interval Index", "Correlation Factor", 0,0, 0,0);
   DrawLegend(Histos,legend,"","P");
   SaveCanvas(c1,SavePath,"Correlation");
   delete c1;


   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = Pred_Expected_Entries;             legend.push_back("Predicted");
   Histos[1] = Pred_Observed_Entries;             legend.push_back("Observed");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Interval Index", "#Tracks", 0,0, 0,0);
   DrawLegend(Histos,legend,"","P");
   SaveCanvas(c1,SavePath,"Prediction_Entries_Absolute");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   TH1D* Diff = (TH1D*)Pred_Observed_Entries->Clone();
   Diff->Add(Pred_Expected_Entries,-1);
   Histos[0] = Diff;                              legend.push_back("Observed-Predicted");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Interval Index", "#Tracks", 0,0, 0,0);
   DrawLegend(Histos,legend,"","P");
   SaveCanvas(c1,SavePath,"Prediction_Entries_Difference");
   delete Diff;
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = Pred_P[0];                         legend.push_back("Predicted");
   DrawSuperposedHistos((TH1**)Histos, legend, "Hist E1",  "P (Gev/c)", "#Tracks", 0,0, 0,0);
   DrawLegend(Histos,legend,"","P");
   SaveCanvas(c1,SavePath,"Prediction_PSpectrum");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = Pred_I[0];                         legend.push_back("Predicted");
   DrawSuperposedHistos((TH1**)Histos, legend, "Hist E1",  "Ionization Variable", "#Tracks", 0,0, 0,0);
   DrawLegend(Histos,legend,"","P");
   SaveCanvas(c1,SavePath,"Prediction_ISpectrum");
   delete c1;

   //////////////////////////////////////////////////     CREATE PLOTS OF MASS DISTRIBUTION

   for(unsigned int s=0;s<signals.size();s++){

   TH1D* DATA = (TH1D*)Data_Mass[0]->Clone();
   TH1D* PRED = (TH1D*)Pred_Mass[0]->Clone();
   TH1D* MCTR = (TH1D*)MCTr_Mass[0]->Clone();
   TH1D* SIGN = (TH1D*)Sign_Mass[s][0]->Clone();
   TH1D* MCSI = (TH1D*)MCTR->Clone();MCSI->Add(SIGN,1);

   PRED->Scale(DATA->Integral()/PRED->Integral());
   MCTR->Scale(DATA->Integral()/MCTR->Integral());

   DATA->Rebin(4);
   PRED->Rebin(4);
   MCTR->Rebin(4);
   SIGN->Rebin(4);
   MCSI->Rebin(4);

   double Max = std::max(DATA->GetMaximum(), PRED->GetMaximum());
   Max        = std::max(MCTR->GetMaximum(), Max);
   Max        = std::max(MCSI->GetMaximum(), Max);
   Max       *= 1.5;

   c1 = new TCanvas("c1","c1,",600,600);

   MCSI->GetXaxis()->SetNdivisions(505);
   MCSI->SetTitle("");
   MCSI->SetStats(kFALSE);
   MCSI->GetXaxis()->SetTitle("Reconstructed Mass (GeV)");
   MCSI->GetYaxis()->SetTitle("#Tracks");
   MCSI->GetYaxis()->SetTitleOffset(1.50);
   MCSI->SetLineColor(46);
   MCSI->SetFillColor(2);
   MCSI->SetMarkerStyle(1);
   MCSI->SetMarkerColor(46);
   MCSI->SetMaximum(Max);
   MCSI->SetMinimum(0.01);
   MCSI->Draw("HIST");
   TH1D* MCSIErr = (TH1D*)MCSI->Clone();
   MCSIErr->SetLineColor(46);
   MCSIErr->Draw("E1 same");

   MCTR->SetLineColor(39);
   MCTR->SetFillColor(64);
   MCTR->SetMarkerStyle(1);
   MCTR->SetMarkerColor(39);
   MCTR->Draw("HIST same");
   TH1D* MCTRErr = (TH1D*)MCTR->Clone();
   MCTRErr->SetLineColor(39);
   MCTRErr->Draw("E1 same");

   PRED->SetMarkerStyle(21);
   PRED->SetMarkerColor(8);
   PRED->SetMarkerSize(1);
   PRED->SetLineColor(8);
   PRED->SetFillColor(0);
   PRED->Draw("E1 same");

   DATA->SetMarkerStyle(20);
   DATA->SetMarkerColor(1);
   DATA->SetMarkerSize(1);
   DATA->SetLineColor(1);
   DATA->SetFillColor(0);
   DATA->Draw("E1 same");

   TLegend* leg = new TLegend(0.79,0.93,0.59,0.73);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(MCSI, "MC+Sign"    ,"F");
   leg->AddEntry(MCTR, "MC"         ,"F");
   leg->AddEntry(PRED, "Prediction" ,"P");
   leg->AddEntry(DATA, "Data"       ,"P");
   leg->Draw();

   SaveCanvas(c1, SavePath, signals[s].Name + "_MassLinear");
   c1->SetLogy(true);
   SaveCanvas(c1, SavePath, signals[s].Name + "_Mass");
   }

   sprintf(Buffer,"%s/Histos.root"  ,SavePath);   
   TFile * outHistos = new TFile(Buffer, "RECREATE");
   Data_Mass[0]->SetName("Mass_Data");
   Data_Mass[0]->Write();
   Pred_Mass[0]->SetName("Mass_Pred");
   Pred_Mass[0]->Write();
   MCTr_Mass[0]->SetName("Mass_MCTr");
   MCTr_Mass[0]->Write();

   for(unsigned int s=0;s<signals.size();s++){
   string tmp = string("Mass_") + signals[s].Name;
   Sign_Mass[s][0]->SetName(tmp.c_str());
   Sign_Mass[s][0]->Write();
   }
   outHistos->Write();
   outHistos->Close();
   delete outHistos;

}


void InitHistos(){
   stPlots_Init(DataPlots,"Data");
   stPlots_Init(MCTrPlots,"MCTr");
   for(unsigned int s=0;s<signals.size();s++){
      stPlots tmp;
      stPlots_Init(tmp,signals[s].Name);
      SignPlots.push_back(tmp);
   }

   Pred_Expected_Entries = new TH1D("Pred_Expected_Entries","Pred_Expected_Entries",40*6,0,40*6);
   Pred_Observed_Entries = new TH1D("Pred_Observed_Entries","Pred_Observed_Entries",40*6,0,40*6);

   Pred_Correlation_A = new TH1D("Pred_Correlation_A","Pred_Correlation_A",40*6,0,40*6);
   Pred_Correlation_B = new TH1D("Pred_Correlation_B","Pred_Correlation_B",40*6,0,40*6);
   Pred_Correlation_C = new TH1D("Pred_Correlation_C","Pred_Correlation_C",40*6,0,40*6);
   Pred_Correlation_D = new TH1D("Pred_Correlation_D","Pred_Correlation_D",40*6,0,40*6);


   Sign_Mass = new TH1D**[signals.size()];
   for(unsigned int s=0;s<signals.size();s++){
      Sign_Mass[s] = new TH1D*[40*6];
   }


   char BckgExt[1024];
   char SignExt[1024];
   char DataExt[1024];
   char MCTrExt[1024];
   char Name   [1024];
   for(unsigned int i=0;i<40*6;i++){
      sprintf(BckgExt,"Bckg");
      GetNameFromIndex(BckgExt, i);
      sprintf(SignExt,"Sign");
      GetNameFromIndex(SignExt, i);
      sprintf(DataExt,"Data");
      GetNameFromIndex(DataExt, i);
      sprintf(MCTrExt,"MCTr");
      GetNameFromIndex(MCTrExt, i);


      sprintf(Name,"CutFinder_I_%s",DataExt);
      Data_I[i]         = new TH1D(Name,Name, 10000,0,dEdxUpLim[dEdxSeleIndex]);
      Data_I[i]->Sumw2();

      sprintf(Name,"CutFinder_Pt_%s",DataExt);
      Data_Pt[i]       = new TH1D(Name,Name,10000,0,PtHistoUpperBound);
      Data_Pt[i]->Sumw2();

      sprintf(Name,"Mass_%s",BckgExt);
      Data_Mass[i]         = new TH1D(Name,Name,200,0,MassHistoUpperBound);
      Data_Mass[i]->Sumw2();

      for(unsigned int s=0;s<signals.size();s++){
         sprintf(SignExt,"%s",signals[s].Name.c_str());
         GetNameFromIndex(SignExt, i);
         sprintf(Name,"Mass_%s",SignExt);
         Sign_Mass[s][i] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
         Sign_Mass[s][i]->Sumw2();
      }

      sprintf(Name,"Mass_%s",DataExt);
      Pred_Mass[i] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
      Pred_Mass[i]->Sumw2();

      sprintf(Name,"Mass_%s",MCTrExt);
      MCTr_Mass[i] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
      MCTr_Mass[i]->Sumw2();

      sprintf(Name,"Pred_I_%s",DataExt);
      Pred_I[i]  = new TH1D(Name,Name,200,0,dEdxUpLim[dEdxMassIndex]);
      Pred_I[i]->Sumw2();

      sprintf(Name,"Pred_P_%s",DataExt);
      Pred_P[i]  = new TH1D(Name,Name,200,0,PtHistoUpperBound);
      Pred_P[i]->Sumw2();

      sprintf(Name,"Data_PI_%s",DataExt);
      Pred_PI[i] = new TH2D(Name,Name,200,0,PtHistoUpperBound, 200, 0, dEdxUpLim[dEdxMassIndex]);
      Pred_PI[i]->Sumw2();

      sprintf(Name,"Data_PI_A_%s",DataExt);
      Data_PI_A[i] = new TH2D(Name,Name,200,0,PtHistoUpperBound, 200, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_A[i]->Sumw2();

      sprintf(Name,"Data_PI_B_%s",DataExt);
      Data_PI_B[i] = new TH2D(Name,Name,200,0,PtHistoUpperBound, 200, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_B[i]->Sumw2();

      sprintf(Name,"Data_PI_C_%s",DataExt);
      Data_PI_C[i] = new TH2D(Name,Name,200,0,PtHistoUpperBound, 200, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_C[i]->Sumw2();

      sprintf(Name,"Data_PI_D_%s",DataExt);
      Data_PI_D[i] = new TH2D(Name,Name,200,0,PtHistoUpperBound, 200, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_D[i]->Sumw2();

      N_A[i] = 0;
      N_B[i] = 0;
      N_C[i] = 0;
      N_D[i] = 0;

      N_Aerr[i] = 0;
      N_Berr[i] = 0;
      N_Cerr[i] = 0;
      N_Derr[i] = 0;
   }
}


void CleanUpHistos(){
   stPlots_Clear(DataPlots);
   stPlots_Clear(MCTrPlots);
   for(unsigned int s=0;s<signals.size();s++){
      stPlots_Clear(SignPlots[s]);
   }SignPlots.clear();

   for(unsigned int i=0;i<40*6;i++){
      delete Data_Pt  [i];
      delete Data_I   [i];   

      delete Pred_Mass   [i];

      delete Data_Mass   [i];
      delete MCTr_Mass   [i];

      for(unsigned int s=0;s<signals.size();s++){
         delete Sign_Mass[s][i];
      }
   }
}

void DrawDEDXVsP(TH2* Histos, char* Title,  char* Xlegend, char* Ylegend, double xmin, double xmax, double ymin, double ymax, bool DrawMassLine)
{
   Histos->SetTitle(Title);
   Histos->SetStats(kFALSE);
   Histos->GetXaxis()->SetTitle(Xlegend);
   Histos->GetYaxis()->SetTitle(Ylegend);
   Histos->GetYaxis()->SetTitleOffset(1.20);
   if(xmin!=xmax)Histos->SetAxisRange(xmin,xmax,"X");
   if(ymin!=ymax)Histos->SetAxisRange(ymin,ymax,"Y");

   Histos->Draw("COLZ");

   if(DrawMassLine){
      for(double i=1;i<=15;i++){
         GetMassLine(i*100)->Draw("same");
      }
   }
}

/*
void PrintEventInRange(double min, double max, FILE* pfile){
  fprintf(pfile,"###########################################\n");
  fprintf(pfile,"Counting events in range [%f,%f]\n",min,max);

  int binMin = Sign_Mass[0]->GetXaxis()->FindBin(min);
  int binMax = Sign_Mass[0]->GetXaxis()->FindBin(max);
  double S = Sign_Mass[0]->Integral(binMin,binMax);
  double B = Data_Mass[0]->Integral(binMin,binMax);


  binMin = Pred_Mass[0]->GetXaxis()->FindBin(min);
  binMax = Pred_Mass[0]->GetXaxis()->FindBin(max);
  double P=0, PL=0, PH=0;
  for(int i=binMin;i<=binMax;i++){
     double pval =  Pred_Mass[0]->GetBinContent(i);
     double perr =  Pred_Mass[0]->GetBinError(i);
     P  += pval;
     PL += pval - perr;
     PH += pval + perr;
  }
//  double P    = Pred_Mass [0]->Integral(binMin,binMax);

  fprintf(pfile,"S = %f  & B = %f & P = %f < %f < %f\n",S,B,PL,P,PH);
  fprintf(pfile,"###########################################\n");
}
*/

double GetMass(double P, double I, bool MC){
   double K = dEdxK_Data[dEdxMassIndex];
   double C = dEdxC_Data[dEdxMassIndex];
   if(MC){
      K = dEdxK_MC[dEdxMassIndex];
      C = dEdxC_MC[dEdxMassIndex];
   }

   return sqrt((I-C)/K)*P;
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

double CutFromEfficiency(TH1* Histo, double Efficiency, bool DoesKeepLeft)
{
   if(DoesKeepLeft){  Efficiency = 1 - Efficiency;  }

   char Buffer[1024];
   sprintf(Buffer,"%s_NTracks",Histo->GetName());
   TH1D* Temp = new TH1D(Buffer,Buffer, Histo->GetXaxis()->GetNbins(), Histo->GetXaxis()->GetXmin(), Histo->GetXaxis()->GetXmax());

   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1);
   Temp->SetBinContent(0,Entries);

   double Integral = 0;
   for(int i=0;i<=Histo->GetXaxis()->GetNbins()+1;i++){
      Integral += Histo->GetBinContent(i);
      if(Integral>Entries)Integral = Entries;
       Temp->SetBinContent(i,   Entries - Integral);      
   }

   unsigned int CutPosition = Temp->GetXaxis()->GetNbins()+1;
   for(int c=0;c<=Temp->GetXaxis()->GetNbins()+1;c++){
      if(Temp->GetBinContent(c)/Entries < Efficiency){ CutPosition = c;  break; }
   }
   delete Temp;

   if(DoesKeepLeft){
      return Histo->GetXaxis()->GetBinLowEdge(CutPosition);
   }else{
      return Histo->GetXaxis()->GetBinUpEdge(CutPosition);
   }
}

double Efficiency(TH1* Histo, double CutX){
   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1);
   double Integral = Histo->Integral(Histo->GetXaxis()->FindBin(CutX),Histo->GetNbinsX()+1);
   return Integral/Entries;
}

double Efficiency(TH2* Histo, double CutX, double CutY){
   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1, 0,Histo->GetNbinsY()+1);
   double Integral = Histo->Integral(Histo->GetXaxis()->FindBin(CutX),Histo->GetNbinsX()+1, Histo->GetYaxis()->FindBin(CutY),Histo->GetNbinsY()+1);
   return Integral/Entries;
}

void GetIndices(int NOM, double Eta, int& HitIndex, int& EtaIndex)
{
   HitIndex = NOM*6;
   EtaIndex = 0;   

         if(fabs(Eta)<0.5)EtaIndex = 1;
   else  if(fabs(Eta)<1.0)EtaIndex = 2;
   else  if(fabs(Eta)<1.5)EtaIndex = 3;
   else  if(fabs(Eta)<2.0)EtaIndex = 4;
   else                   EtaIndex = 5;
}

void FillArray(int HitIndex, int EtaIndex, double* Array, double value){
   Array[ 0                 ] +=  value;
   Array[           EtaIndex] +=  value;
   Array[HitIndex           ] +=  value;
   Array[HitIndex + EtaIndex] +=  value;
}

void FillHisto(int HitIndex, int EtaIndex, TH1D** Histo, double value, double weight){
   Histo[ 0                 ]->Fill(value,weight);
   Histo[           EtaIndex]->Fill(value,weight);
   Histo[HitIndex           ]->Fill(value,weight);
   Histo[HitIndex + EtaIndex]->Fill(value,weight);
}

void FillHisto(int HitIndex, int EtaIndex, TH2D** Histo, double value1, double value2, double weight){
   Histo[ 0                 ]->Fill(value1,value2,weight);
   Histo[           EtaIndex]->Fill(value1,value2,weight);
   Histo[HitIndex           ]->Fill(value1,value2,weight);
   Histo[HitIndex + EtaIndex]->Fill(value1,value2,weight);
}

int GetCutIndex(int HitIndex, int EtaIndex){
   int CutIndex;
   if(SplitMode==0){
      CutIndex = 0;
   }else if(SplitMode==1){
      CutIndex = HitIndex;
   }else{
      CutIndex = HitIndex + EtaIndex;
   }
   return CutIndex;
}

void GetNameFromIndex(char* NameExt, int index)
{
      unsigned int Hit = index/6;
      unsigned int Eta = index%6;
      if(Hit>=1)sprintf(NameExt,"%s_SSHit%02i",NameExt,Hit);
      if(Eta==1)sprintf(NameExt,"%s_Eta00to05",NameExt);
      if(Eta==2)sprintf(NameExt,"%s_Eta05to10",NameExt);
      if(Eta==3)sprintf(NameExt,"%s_Eta10to15",NameExt);
      if(Eta==4)sprintf(NameExt,"%s_Eta15to20",NameExt);
      if(Eta==5)sprintf(NameExt,"%s_Eta20to25",NameExt);
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
   double deta = eta1 - eta2;
   double dphi = phi1 - phi2;
   while (dphi >   M_PI) dphi -= 2*M_PI;
   while (dphi <= -M_PI) dphi += 2*M_PI;
   return sqrt(deta*deta + dphi*dphi);
}

double DistToHSCP (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest){
   reco::TrackRef   trackRef = hscp.trackRef(); if(trackRef.isNull())return false;
   reco::Track      track    = *trackRef;

   double RMin = 9999; IndexOfClosest=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(abs(genColl[g].pdgId())<1000000)continue;
      if(genColl[g].status()!=1)continue;
      if(genColl[g].pt()<5)continue;

      double dR = deltaR(track.eta(), track.phi(), genColl[g].eta(), genColl[g].phi());
      if(dR<RMin){RMin=dR;IndexOfClosest=g;}
   }
   return RMin;
}

void SetWeight(double IntegratedLuminosityInPb, double CrossSection, double MCEvents){
   if(IntegratedLuminosityInPb>0){
      Event_Weight = (CrossSection * IntegratedLuminosityInPb) / MCEvents;
   }else{
      Event_Weight=1;
   }
}




