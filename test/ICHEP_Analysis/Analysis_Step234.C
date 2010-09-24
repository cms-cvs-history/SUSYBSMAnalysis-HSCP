
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
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}


#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;
#endif


#include "Analysis_Global.h"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"
#include "Analysis_Samples.h"


/////////////////////////// FUNCTION DECLARATION /////////////////////////////


//void FillArray(int HitIndex, int EtaIndex, double* Array, double value);
//void FillHisto(int HitIndex, int EtaIndex, TH1D**  Histo, double value                , double weight=1);
//void FillHisto(int HitIndex, int EtaIndex, TH2D**  Histo, double value1, double value2, double weight=1);

//double CutFromEfficiency(TH1* Histo, double Efficiency, bool DoesKeepLeft=false);
//double Efficiency(TH1* Histo, double CutX);
//double Efficiency(TH2* Histo, double CutX, double CutY);

void CompleteAnalysis(char* SavePath);
void Analysis_Step2();
void Analysis_Step3();
void Analysis_Step4(char* SavePath);

void InitHistos();
void CleanUpHistos();

double DistToHSCP     (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest);
int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl);
bool   isGoodCandidate(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double PtCut=0, double ICut=0, stPlots* st=NULL, double PtRescale=1.0, double IRescale=1.0);
void DumpCandidateInfo(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, FILE* pFile);
bool PassTrigger      (const fwlite::ChainEvent& ev);
bool hasGoodPtHat     (const fwlite::ChainEvent& ev, double PtMax);

void SetWeight(double IntegratedLuminosityInPb=-1, double CrossSection=0, double MCEvents=0);
void SetWeightMC(double IntegratedLuminosityInPb, double SampleEquivalentLumi, double SampleSize, double MaxEvent);

bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold=false);

/////////////////////////// VARIABLE DECLARATION /////////////////////////////

TH1D*  Data_Pt[40*6];
TH1D*  Data_I [40*6];

TH1D*  Data_Mass[40*6];
TH1D*** Sign_Mass;
TH1D** Sign_Mass_Syst_PtLow;
TH1D** Sign_Mass_Syst_ILow;
TH1D*  Pred_Mass[40*6];
TH1D*  MCTr_Mass[40*6];

double N_A[40*6];	double N_Aerr[40*6];
double N_B[40*6];	double N_Berr[40*6];
double N_C[40*6];	double N_Cerr[40*6];
double N_D[40*6];	double N_Derr[40*6];

TH1D*  Pred_P    [40*6];
TH1D*  Pred_I    [40*6];
TH2D*  Data_PI_A [40*6];
TH2D*  Data_PI_B [40*6];
TH2D*  Data_PI_C [40*6];
TH2D*  Data_PI_D [40*6];
TH2D*  Pred_PI   [40*6];
TH1D*  Ctrl_BckgP;
TH1D*  CtrlPt_BckgIs;
TH1D*  CtrlPt_BckgIm;
TH1D*  CtrlP_BckgIs;
TH1D*  CtrlP_BckgIm;
TH1D*  Ctrl_SignP;
TH1D*  CtrlPt_SignIs;
TH1D*  CtrlPt_SignIm;
TH1D*  CtrlP_SignIs;
TH1D*  CtrlP_SignIm;


TH1D* Pred_Expected_Entries;
TH1D* Pred_Observed_Entries;
TH1D* Pred_Correlation_A;
TH1D* Pred_Correlation_B;
TH1D* Pred_Correlation_C;
TH1D* Pred_Correlation_D;


stPlots DataPlots;
std::vector<stPlots> SignPlots;
std::vector<stPlots> MCPlots;
stPlots MCTrPlots;

std::vector<stSignal> signals;
std::vector<stMC>     MCsample;

string LegendTitle;

/////////////////////////// CODE PARAMETERS /////////////////////////////

std::vector<string> DataFileName;
string TreeName;

float Event_Weight = 1;
int MaxEntry = -1;

void Analysis_Step234(string MODE="COMPILE", double WP_Pt=-1.0, double WP_I=-1, int SplitMode_=2, int dEdxSel_=0, int dEdxMass_=0, int TypeMode_=0)
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
   gStyle->SetNdivisions(505);

   GetSignalDefinition(signals);
   GetMCDefinition(MCsample);

   char Buffer[2048];   
   char Command[2048];
   DataFileName.clear();
   GetInputFiles(DataFileName, "Data");
//   DataFileName.push_back("InputFiles/Data.root");

   dEdxSeleIndex = dEdxSel_;
   dEdxMassIndex = dEdxMass_;

   TreeName = "HSCPTreeBuilder/MyTree";

   TypeMode  = TypeMode_;
   SplitMode = SplitMode_;
   if(SplitMode>0){    GlobalMinNOH = 1;
   }else{         GlobalMinNOH = 9;   }

   if(TypeMode==0){
      LegendTitle = "Tracker - Only";
   }else{
      LegendTitle = "Tracker + Muon";
   }

   DefaultCutPt   = 0;
   DefaultCutI    = 0;
   SelectionCutPt = pow(10,WP_Pt);
   SelectionCutI  = pow(10,WP_I);

   sprintf(Buffer,"Results/"       );                                  sprintf(Command,"mkdir %s",Buffer); system(Command);
//   sprintf(Buffer,"%s%s/"         ,Buffer,MODE.c_str());               sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sSplitMode%i/",Buffer,SplitMode);                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sMinHit%02i/" ,Buffer,GlobalMinNOH);               sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sSele_%s/"    ,Buffer,dEdxLabel[dEdxSeleIndex]);   sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sMass_%s/"    ,Buffer,dEdxLabel[dEdxMassIndex]);   sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sType%i/"     ,Buffer,TypeMode);                   sprintf(Command,"mkdir %s",Buffer); system(Command);


/*   printf("Running Mode = %s\n",MODE.c_str());  
   if(MODE==string("MERGE_MAP")){
      Merge_Map(Buffer,  0,1000);
      Merge_Map(Buffer, 50,1000);
      Merge_Map(Buffer, 75,1000);
      Merge_Map(Buffer,100,1000);
      Merge_Map(Buffer,125,1000);
   }else if(MODE==string("ANALYSE")){*/
      sprintf(Buffer,"%sWPPt%+03i/"  ,Buffer,(int)(10*log10(SelectionCutPt)));   sprintf(Command,"mkdir %s",Buffer); system(Command);
      sprintf(Buffer,"%sWPI%+03i/"   ,Buffer,(int)(10*log10(SelectionCutI)));    sprintf(Command,"mkdir %s",Buffer); system(Command);
      CompleteAnalysis(Buffer);
//   }else{
//      printf("UNKNOWN MODE\n");
//   }

   return;
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

   CleanUpHistos();


   tmp->Write();
   tmp->Close();

}

bool hasGoodPtHat(const fwlite::ChainEvent& ev, double PtMax){
   if(PtMax<0)return true;
   fwlite::Handle< GenEventInfoProduct > genInfo;
   genInfo.getByLabel(ev, "generator");
   if(!genInfo.isValid()){printf("genInfo NotFound\n");return false;}
   if((genInfo->binningValues()[0])<PtMax)return true;
   return false;
}


bool isGoodCandidate(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double PtCut, double ICut, stPlots* st, double PtRescale, double IRescale)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;

   reco::TrackRef   trackRef = hscp.trackRef(); if(trackRef.isNull())return false;
   reco::Track      track    = *trackRef;

   if(st){st->WN_Total+=Event_Weight;	st->UN_Total++;}

   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
   std::vector<reco::Vertex> vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

   double dz  = track.dz (vertexColl[0].position());
   double dxy = track.dxy(vertexColl[0].position());
   for(unsigned int i=1;i<vertexColl.size();i++){
      if(fabs(track.dz (vertexColl[i].position())) < fabs(dz) ){
         dz  = track.dz (vertexColl[i].position());
         dxy = track.dxy(vertexColl[i].position());
      }
   }

   if(st){st->BS_Hits->Fill(track.found(),Event_Weight);}
   if(track.found()<GlobalMinNOH)return false;
   if(hscp.dedx(dEdxSeleIndex).numberOfMeasurements()<GlobalMinNOM)return false;
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
   if(track.pt()*PtRescale<GlobalMinPt)return false;
   if(st){st->AS_MPt ->Fill(track.pt(),Event_Weight);}
   if(st){st->WN_MPt   +=Event_Weight;   st->UN_MPt ++;}

   if(st){st->BS_MIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->BS_MIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}
   if(hscp.dedx(dEdxSeleIndex).dEdx()*IRescale<GlobalMinI)return false;
   if(st){st->AS_MIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->AS_MIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}
   if(st){st->WN_MI   +=Event_Weight;   st->UN_MI++;}

   if(st){st->BS_Pt  ->Fill(track.pt(),Event_Weight);}
   if(st){st->BS_Is ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->BS_Im ->Fill(hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}

   if(st){st->BS_EtaP ->Fill(track.eta(),track.p(),Event_Weight);}
   if(st){st->BS_EtaPt->Fill(track.eta(),track.pt(),Event_Weight);}
   if(st){st->BS_PIs  ->Fill(track.p(),hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->BS_PIm  ->Fill(track.p(),hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}

   if(track.pt()*PtRescale<PtCut)return false;
   if(st){st->WN_Pt    +=Event_Weight;   st->UN_Pt ++;}
   if(hscp.dedx(dEdxSeleIndex).dEdx()*IRescale<ICut)return false;
   if(st){st->WN_I    +=Event_Weight;   st->UN_I++;}

   if(st){st->AS_Pt  ->Fill(track.pt(),Event_Weight);}
   if(st){st->AS_Is ->Fill(hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->AS_Im ->Fill(hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}

   if(st){st->AS_EtaP ->Fill(track.eta(),track.p(),Event_Weight);}
   if(st){st->AS_EtaPt->Fill(track.eta(),track.pt(),Event_Weight);}
   if(st){st->AS_PIs  ->Fill(track.p(),hscp.dedx(dEdxSeleIndex).dEdx(),Event_Weight);}
   if(st){st->AS_PIm  ->Fill(track.p(),hscp.dedx(dEdxMassIndex).dEdx(),Event_Weight);}

   return true;
}

void DumpCandidateInfo(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, FILE* pFile)
{
   reco::MuonRef  muon  = hscp.muonRef();
   reco::TrackRef track = hscp.trackRef();
   if(track.isNull())return;

   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return;}
   std::vector<reco::Vertex> vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return;}
   const reco::Vertex& vertex = vertexColl[0];

   double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
   double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
   double Mass = GetMass(PBinned,IBinned);   
   double MassExact = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
//   if(Mass<MinCandidateMass)return;
   double dz  = track->dz (vertex.position());
   double dxy = track->dxy(vertex.position());

   int HitIndex,EtaIndex;
   GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
   int CutIndex = GetCutIndex(HitIndex,EtaIndex);

   fprintf(pFile,"\n");
   fprintf(pFile,"---------------------------------------------------------------------------------------------------\n");
   fprintf(pFile,"Candidate Type = %i --> Mass (Binned): %7.2f GeV  Mass (UnBinned): %7.2f\n",hscp.type(),Mass, MassExact);
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
   fprintf(pFile,"\n");
}

bool PassTrigger(const fwlite::ChainEvent& ev)
{
      //THIS IS ONLY USE FOR SIGNAL SAMPLE (TRIGGER IS ALREADY APPLIED IN BOTH DATA AND MC!)

      //Some run on data have a different trigger table... (they are 261Events on 871561)
      if(ev.eventAuxiliary().run()>=132440 && ev.eventAuxiliary().run()<=132528){
         std::cout << "Using special event range, maybe crashing\n";
      }//return false;

      edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
      if(!tr.isValid()){
         tr = ev.triggerResultsByName("REDIGI36X");
         if(!tr.isValid())return false;
      }

      if(tr.accept("HLT_DoubleMu3")                                                    ) return true;
      if(tr.accept("HLT_Mu9")                                                          ) return true;
      if(tr.accept("HLT_MET100")                                                       ) return true;

      fwlite::Handle< trigger::TriggerEvent > trEvHandle;
      trEvHandle.getByLabel(ev,"hltTriggerSummaryAOD");
      trigger::TriggerEvent trEv = *trEvHandle;

      if(IncreasedTreshold(trEv, InputTag("hlt1jet50U"        ,"","HLT"), 100, 1      )) return true;
      if(IncreasedTreshold(trEv, InputTag("hltDiJetAve30U8E29","","HLT"),  70, 2, true)) return true;
      if(IncreasedTreshold(trEv, InputTag("hlt4jet15U"        ,"","HLT"),  25, 4      )) return true;
      return false;
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
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
//      if(!PassTrigger(tree) )continue;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(tree,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(!isGoodCandidate(hscp,tree))continue;

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
         if(CutPt[i]<GlobalMinPt)CutPt[i]=GlobalMinPt;
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
         if(CutI[i]<GlobalMinI)CutI[i]=GlobalMinI;
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
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
//      if(!PassTrigger(tree) )continue;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(tree,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(!isGoodCandidate(hscp,tree))continue;

         int HitIndex, EtaIndex;
         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
         
         bool PassMinPt = GlobalMinPt;
         bool PassPtCut = track->pt()>=CutPt[CutIndex];
         bool PassMinI  = (hscp.dedx(dEdxSeleIndex).dEdx()>=GlobalMinI);
         bool PassICut  = (hscp.dedx(dEdxSeleIndex).dEdx()>=CutI[CutIndex]);


         if(PassMinPt && track->pt()<25){
            CtrlPt_BckgIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(), Event_Weight);
            CtrlPt_BckgIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
         }
         if(PassMinPt &&  track->pt()>25){
            CtrlPt_SignIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(), Event_Weight);
            CtrlPt_SignIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
         }

         if(track->p()>10 && track->p()<25){
            CtrlP_BckgIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(), Event_Weight);
            CtrlP_BckgIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
         }
         if(track->p()>10 &&  track->p()>25){
            CtrlP_SignIs->Fill(hscp.dedx(dEdxSeleIndex).dEdx(), Event_Weight);
            CtrlP_SignIm->Fill(hscp.dedx(dEdxMassIndex).dEdx(), Event_Weight);
         }


         if(PassMinI && hscp.dedx(dEdxSeleIndex).dEdx()<0.5){
            Ctrl_BckgP->Fill(track->p(), Event_Weight);
         }
         if(PassMinI && hscp.dedx(dEdxSeleIndex).dEdx()>0.5){
            Ctrl_SignP->Fill(track->p(), Event_Weight);
         }




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
   FILE* pFileTrg = NULL;

   //////////////////////////////////////////////////     BUILD BACKGROUND MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_D_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");

      sprintf(Buffer,"%s/Candidate_D_Trigger.txt",SavePath);
      pFileTrg = fopen(Buffer,"w");
   }


   unsigned int COUNT_PRINTED_DATA=0;
   fwlite::ChainEvent treeD(DataFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Building Mass Spectrum for D :");
   TreeStep = treeD.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t ientry=0;ientry<treeD.size();ientry++){
      treeD.to(ientry);
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}

      DataPlots.WN_TotalE+=Event_Weight;       DataPlots.UN_TotalE++;
//      if(!PassTrigger(treeD) )continue;
      DataPlots.WN_TotalTE+=Event_Weight;      DataPlots.UN_TotalTE++;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeD,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      bool HSCPTk = false;
      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;
 
         GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
         if(!isGoodCandidate(hscp,treeD,CutPt[CutIndex], CutI[CutIndex], &DataPlots))continue;

         //DEBUG
         double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
         double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
         double Mass = GetMass(PBinned,IBinned);
         if(!isnan(Mass)){
 
            HSCPTk = true;
//          double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx());
            FillHisto(HitIndex, EtaIndex,Data_Mass , Mass, Event_Weight);
            if(SavePath && (Mass>=MinCandidateMass || COUNT_PRINTED_DATA<1000) ){DumpCandidateInfo(hscp, treeD, pFile);   COUNT_PRINTED_DATA++;}
            if(SavePath)fprintf(pFileTrg,"Run=%i Lumi=%i Event=%i\n",treeD.eventAuxiliary().run(),treeD.eventAuxiliary().luminosityBlock(),treeD.eventAuxiliary().event());
         }
      } // end of Track Loop
      if(HSCPTk){DataPlots.WN_HSCPE+=Event_Weight;  DataPlots.UN_HSCPE++;          }
   }// end of Event Loop
   printf("\n");
   if(pFile){fclose(pFile);pFile=NULL;};
   if(pFileTrg){fclose(pFileTrg);pFileTrg=NULL;};

   //////////////////////////////////////////////////     BUILD MCTRUTH MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_M_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");
   }

   for(unsigned int m=0;m<MCsample.size();m++){
      unsigned int COUNT_PRINTED_SIGNAL=0;
      fprintf(pFile,"\nXXXXXXXXXXXXXXXXXXXXXXX %10s XXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n",MCsample[m].Name.c_str());

      if(SavePath){
         char Buffer[2048];
         sprintf(Buffer,"%s/Candidate_M_%s_Trigger.txt",SavePath,MCsample[m].Name.c_str());
         pFileTrg = fopen(Buffer,"w");
      }

      std::vector<string> FileName;
      GetInputFiles(FileName, MCsample[m].Name);

      fwlite::ChainEvent treeM(FileName);
      SetWeightMC(IntegratedLuminosity,MCsample[m].ILumi, treeM.size(), MCsample[m].MaxEvent);
      printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
      printf("Building Mass for %10s :",MCsample[m].Name.c_str());
      TreeStep = treeM.size()/50;if(TreeStep==0)TreeStep=1;
      for(Long64_t ientry=0;ientry<treeM.size();ientry++){       
          treeM.to(ientry);
         if(MaxEntry>0 && ientry>MaxEntry)break;
         if(MCsample[m].MaxEvent>0 && ientry>MCsample[m].MaxEvent)break;
         if(ientry%TreeStep==0){printf(".");fflush(stdout);}

         if(!hasGoodPtHat(treeM, MCsample[m].MaxPtHat)){continue;}


         MCTrPlots .WN_TotalE+=Event_Weight;       MCTrPlots.UN_TotalE++;
         MCPlots[m].WN_TotalE+=Event_Weight;       MCPlots[m].UN_TotalE++;
//         if(!PassTrigger(treeM) )continue;
         MCTrPlots .WN_TotalTE+=Event_Weight;      MCTrPlots.UN_TotalTE++;
         MCPlots[m].WN_TotalTE+=Event_Weight;      MCPlots[m].UN_TotalTE++;

         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeM,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

         bool HSCPTk = false;
         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

            GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
            int CutIndex = GetCutIndex(HitIndex,EtaIndex);
                isGoodCandidate(hscp,treeM,CutPt[CutIndex], CutI[CutIndex], &MCPlots[m]);
            if(!isGoodCandidate(hscp,treeM,CutPt[CutIndex], CutI[CutIndex], &MCTrPlots))continue;

            //DEBUG
            double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
            double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
            double Mass = GetMass(PBinned,IBinned);
            if(!isnan(Mass)){
               HSCPTk = true;
//             double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
               FillHisto(HitIndex, EtaIndex, MCTr_Mass, Mass, Event_Weight);


               if(SavePath)fprintf(pFileTrg,"Run=%i Lumi=%i Event=%i\n",treeM.eventAuxiliary().run(),treeM.eventAuxiliary().luminosityBlock(),treeM.eventAuxiliary().event());
               if(SavePath && COUNT_PRINTED_SIGNAL<10){
                  if(SavePath && Mass>=MinCandidateMass)DumpCandidateInfo(hscp, treeM, pFile);
                  COUNT_PRINTED_SIGNAL++;
               }
            }
         } // end of Track Loop 
         if(HSCPTk){MCTrPlots .WN_HSCPE+=Event_Weight; MCTrPlots .UN_HSCPE++;          }
         if(HSCPTk){MCPlots[m].WN_HSCPE+=Event_Weight; MCPlots[m].UN_HSCPE++;          }
      }// end of Event Loop
      printf("\n");
      if(pFileTrg){fclose(pFileTrg);pFileTrg=NULL;};
   }
   if(pFile){fclose(pFile);pFile=NULL;};


   //////////////////////////////////////////////////     BUILD SIGNAL MASS SPECTRUM


   for(unsigned int s=0;s<signals.size();s++){
      unsigned int COUNT_PRINTED_SIGNAL=0;
      if(SavePath){
         char Buffer[2048];
         sprintf(Buffer,"%s/Candidate_%s_Dump.txt",SavePath,signals[s].Name.c_str());
         pFile = fopen(Buffer,"w");

         sprintf(Buffer,"%s/Candidate_%s_Trigger.txt",SavePath,signals[s].Name.c_str());
         pFileTrg = fopen(Buffer,"w");
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
         if(MaxEntry>0 && ientry>MaxEntry)break;
         if(ientry%TreeStep==0){printf(".");fflush(stdout);}

         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         genCollHandle.getByLabel(treeS, "genParticles");
         if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
         std::vector<reco::GenParticle> genColl = *genCollHandle;
         if(signals[s].NChargedHSCP>=0 && HowManyChargedHSCP(genColl)!=signals[s].NChargedHSCP)continue;

         SignPlots[s].WN_TotalE+=Event_Weight;       SignPlots[s].UN_TotalE++;
         if(!PassTrigger(treeS) )continue;
         SignPlots[s].WN_TotalTE+=Event_Weight;      SignPlots[s].UN_TotalTE++;

         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeS,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;


         bool HSCPTk = false;
         bool HSCPTkSystA = false;
         bool HSCPTkSystB = false;
         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

            int ClosestGen;
            if(DistToHSCP(hscp, genColl, ClosestGen)>0.03)continue;

            GetIndices(hscp.dedx(dEdxSeleIndex).numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
            int CutIndex = GetCutIndex(HitIndex,EtaIndex);

            //FOR SYSTEMATIC COMPUTATION (START)
            //A Signal Pt(&P) -->0.95*Pt(&P)
            if(isGoodCandidate(hscp,treeS,CutPt[CutIndex], CutI[CutIndex], NULL, 0.95, 1.0)){
               SignPlots[s].WN_I_SYSTA+=Event_Weight; SignPlots[s].UN_I_SYSTA++;
               double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()*0.95));
               double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
               double Mass    = GetMass(PBinned,IBinned);
               if(!isnan(Mass)){
                  HSCPTkSystA    = true;
                  Sign_Mass_Syst_PtLow[s]->Fill(Mass, Event_Weight);
               }
            }

            //B Signal I -->0.95*I
            if(isGoodCandidate(hscp,treeS,CutPt[CutIndex], CutI[CutIndex], NULL, 1.0, 0.95)){
               SignPlots[s].WN_I_SYSTB+=Event_Weight; SignPlots[s].UN_I_SYSTB++;
               double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
               double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()*0.95));
               double Mass    = GetMass(PBinned,IBinned);
               if(!isnan(Mass)){
                  HSCPTkSystB    = true;
                  Sign_Mass_Syst_ILow[s]->Fill(Mass, Event_Weight);
               }
            }
            //FOR SYSTEMATIC COMPUTATION (END)


            if(!isGoodCandidate(hscp,treeS,CutPt[CutIndex], CutI[CutIndex], &SignPlots[s]))continue;         

            //DEBUG
            double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
            double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(hscp.dedx(dEdxMassIndex).dEdx()));
            double Mass = GetMass(PBinned,IBinned);
            if(!isnan(Mass)){

               HSCPTk = true;

//             double Mass = GetMass(track->p(),hscp.dedx(dEdxMassIndex).dEdx(), true);
               FillHisto(HitIndex, EtaIndex, Sign_Mass[s], Mass, Event_Weight);

               if(SavePath)fprintf(pFileTrg,"Run=%i Lumi=%i Event=%i\n",treeS.eventAuxiliary().run(),treeS.eventAuxiliary().luminosityBlock(),treeS.eventAuxiliary().event());
               if(SavePath && (Mass>=MinCandidateMass || COUNT_PRINTED_SIGNAL<10)){
                  DumpCandidateInfo(hscp, treeS, pFile);
                  COUNT_PRINTED_SIGNAL++;
               }
            }
         } // end of Track Loop 
         if(HSCPTk     ){SignPlots[s].WN_HSCPE      +=Event_Weight;  SignPlots[s].UN_HSCPE++;          }
         if(HSCPTkSystA){SignPlots[s].WN_HSCPE_SYSTA+=Event_Weight;  SignPlots[s].UN_HSCPE_SYSTA++;    }
         if(HSCPTkSystB){SignPlots[s].WN_HSCPE_SYSTB+=Event_Weight;  SignPlots[s].UN_HSCPE_SYSTB++;    }
       }// end of Event Loop
      printf("\n");
      if(pFile){fclose(pFile);pFile=NULL;};
      if(pFileTrg){fclose(pFileTrg);pFileTrg=NULL;};
   }// end of signal Type loop




   //////////////////////////////////////////////////     DUMP USEFUL INFORMATION

   char Buffer[2048];
   sprintf(Buffer,"%s/CUT_Dump.txt",SavePath);
   pFile = fopen(Buffer,"w");
   fprintf(pFile,"MODE          = %i\n",SplitMode);
   fprintf(pFile,"Selection     = %s\n",dEdxLabel[dEdxSeleIndex]);
   fprintf(pFile,"Mass          = %s\n",dEdxLabel[dEdxMassIndex]);
   fprintf(pFile,"WP PT         = %4.3E\n",SelectionCutPt);
   fprintf(pFile,"WP I          = %4.3E\n",SelectionCutI);
   fprintf(pFile,"GlobalMinNOH  = %02i\n",GlobalMinNOH);
   fprintf(pFile,"GlobalMinNOM  = %02i\n",GlobalMinNOM);
   fprintf(pFile,"GlobalMaxChi2 = %6.2f\n",GlobalMaxChi2);
   fprintf(pFile,"--------------------\n");

   double CutMin_I  = 9999;   double CutMax_I  = 0;  double CutMean_I  = 0;
   double CutMin_Pt = 9999;   double CutMax_Pt = 0;  double CutMean_Pt = 0;
   int NCutsI=0;  int NCutsPt=0;
   for(unsigned int i=0;i<40*6;i++){
      if(SplitMode==0 && i>0)continue;
      if(SplitMode==1 && (i==0 || i%6!=0))continue;
      if(SplitMode==2 && (i< 6 || i%6==0))continue;
      

      if(CutI [i]<CutMin_I                                                         )CutMin_I =CutI [i];
      if(CutI [i]>CutMax_I  && CutI [i]<Data_I [0]->GetXaxis()->GetXmax())CutMax_I =CutI [i];
      if(CutI [i]<Data_I [0]->GetXaxis()->GetXmax()){CutMean_I+=CutI [i];NCutsI++;}
      if(CutPt[i]<CutMin_Pt                                                        )CutMin_Pt=CutPt[i];
      if(CutPt[i]>CutMax_Pt && CutPt[i]<Data_Pt[0]->GetXaxis()->GetXmax())CutMax_Pt=CutPt[i];
      if(CutPt[i]<Data_Pt[0]->GetXaxis()->GetXmax()){CutMean_Pt+=CutPt[i];NCutsPt++;}
   
      char IntervalName[2048];
      sprintf(IntervalName," ");//Just here to initialize the char*
      GetNameFromIndex(IntervalName,i);
      fprintf(pFile,"CutIndex=%03i  %20s PtCut=%14.5f   ICut=%14.5f\n",i,IntervalName,CutPt[i],CutI[i]);
   }
   CutMean_I /= NCutsI;   CutMean_Pt /= NCutsPt;   

   DataPlots.MeanICut = CutMean_I;
   MCTrPlots.MeanICut = CutMean_I;
   for(unsigned int s=0;s<signals.size();s++){
      SignPlots[s].MeanICut = CutMean_I;
   }

   DataPlots.MeanPtCut = CutMean_Pt;
   MCTrPlots.MeanPtCut = CutMean_Pt;
   for(unsigned int s=0;s<signals.size();s++){
      SignPlots[s].MeanPtCut = CutMean_Pt;
   }

   fprintf(pFile,"--------------------\n");
   fprintf(pFile,"PtCut Range =[%10.5f,%10.5f]  Mean = %10.5f\n",CutMin_Pt,CutMax_Pt, CutMean_Pt);
   fprintf(pFile,"ICut  Range =[%10.5f,%10.5f]  Mean = %10.5f\n",CutMin_I ,CutMax_I,  CutMean_I);
   fprintf(pFile,"--------------------\n");

   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"DATA SELECTION DETAILS\n");
   fprintf(pFile,"--------------------\n");
   stPlots_Dump(DataPlots, pFile);

   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"MC TRUTH SELECTION DETAILS\n");
   fprintf(pFile,"--------------------\n");
   stPlots_Dump(MCTrPlots, pFile);

   for(unsigned int m=0;m<MCsample.size();m++){
      fprintf(pFile,"##### ##### %10s ##### #####\n",MCsample[m].Name.c_str());
      stPlots_Dump(MCPlots[m], pFile);
   }
   
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

   fprintf(pFile,"\nIntegral in range [0,1000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D",GetEventInRange(0,1000,Data_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","P",GetEventInRange(0,1000,Pred_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","M",GetEventInRange(0,1000,MCTr_Mass[0]));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(0,1000,Sign_Mass[s][0]));
   }
   fprintf(pFile,"\nIntegral in range [75,1000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D",GetEventInRange(75,1000,Data_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","P",GetEventInRange(75,1000,Pred_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","M",GetEventInRange(75,1000,MCTr_Mass[0]));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(75,1000,Sign_Mass[s][0]));
   }
   fprintf(pFile,"\nIntegral in range [100,1000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","S",GetEventInRange(100,1000,Data_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","P",GetEventInRange(100,1000,Pred_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","M",GetEventInRange(100,1000,MCTr_Mass[0]));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(100,1000,Sign_Mass[s][0]));
   }
   fprintf(pFile,"\nIntegral in range [125,1000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D",GetEventInRange(125,1000,Data_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","P",GetEventInRange(125,1000,Pred_Mass[0]));
   fprintf(pFile,"%15s = %5.3E\n","M",GetEventInRange(125,1000,MCTr_Mass[0]));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(125,1000,Sign_Mass[s][0]));
   }
   fprintf(pFile,"--------------------\n");
   fclose(pFile);

   //////////////////////////////////////////////////     CREATE EFFICIENCY FILE

   sprintf(Buffer,"%s/Aeff.tmp",SavePath);
   pFile = fopen(Buffer,"w");
   for(unsigned int s=0;s<signals.size();s++){
      fprintf(pFile,"%15s Eff=%4.3E (%4.3E)",signals[s].Name.c_str(),SignPlots[s].WN_I        /(  SignPlots[s].WN_HSCPE),      SignPlots[s].UN_I       /(  SignPlots[s].UN_HSCPE  ));
      fprintf(pFile,"SYSTA: Eff=%4.3E (%4.3E)",                      SignPlots[s].WN_I_SYSTA  /(  SignPlots[s].WN_HSCPE_SYSTA),SignPlots[s].UN_I_SYSTA /(  SignPlots[s].UN_HSCPE_SYSTA  ));
      fprintf(pFile,"SYSTB: Eff=%4.3E (%4.3E)",                      SignPlots[s].WN_I_SYSTB  /(  SignPlots[s].WN_HSCPE_SYSTB),SignPlots[s].UN_I_SYSTB /(  SignPlots[s].UN_HSCPE_SYSTB  ));
      fprintf(pFile,"\n");
   }
   fclose(pFile);

   fflush(stdout);

}


void InitHistos(){
   stPlots_Init(DataPlots,"Data");
   stPlots_Init(MCTrPlots,"MCTr");
   for(unsigned int m=0;m<MCsample.size();m++){
      stPlots tmp;
      stPlots_Init(tmp,MCsample[m].Name);
      MCPlots.push_back(tmp);
   }

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

   Ctrl_BckgP    = new TH1D("Ctrl_BckgP"   ,"Ctrl_BckgP"   ,200,0,PtHistoUpperBound);	      Ctrl_BckgP->Sumw2();
   CtrlPt_BckgIs = new TH1D("CtrlPt_BckgIs","CtrlPt_BckgIs",200,0,dEdxUpLim[dEdxSeleIndex]);  CtrlPt_BckgIs->Sumw2();
   CtrlPt_BckgIm = new TH1D("CtrlPt_BckgIm","CtrlPt_BckgIm",200,0,dEdxUpLim[dEdxMassIndex]);  CtrlPt_BckgIm->Sumw2();
   CtrlP_BckgIs  = new TH1D("CtrlP_BckgIs" ,"CtrlP_BckgIs" ,200,0,dEdxUpLim[dEdxSeleIndex]);  CtrlP_BckgIs->Sumw2();
   CtrlP_BckgIm  = new TH1D("CtrlP_BckgIm" ,"CtrlP_BckgIm" ,200,0,dEdxUpLim[dEdxMassIndex]);  CtrlP_BckgIm->Sumw2();
   Ctrl_SignP    = new TH1D("Ctrl_SignP"   ,"Ctrl_SignP"   ,200,0,PtHistoUpperBound);         Ctrl_SignP->Sumw2();
   CtrlPt_SignIs = new TH1D("CtrlPt_SignIs","CtrlPt_SignIs",200,0,dEdxUpLim[dEdxSeleIndex]);  CtrlPt_SignIs->Sumw2();
   CtrlPt_SignIm = new TH1D("CtrlPt_SignIm","CtrlPt_SignIm",200,0,dEdxUpLim[dEdxMassIndex]);  CtrlPt_SignIm->Sumw2();
   CtrlP_SignIs  = new TH1D("CtrlP_SignIs" ,"CtrlP_SignIs" ,200,0,dEdxUpLim[dEdxSeleIndex]);  CtrlP_SignIs->Sumw2();
   CtrlP_SignIm  = new TH1D("CtrlP_SignIm" ,"CtrlP_SignIm" ,200,0,dEdxUpLim[dEdxMassIndex]);  CtrlP_SignIm->Sumw2();

   Sign_Mass = new TH1D**[signals.size()];
   for(unsigned int s=0;s<signals.size();s++){
      Sign_Mass[s] = new TH1D*[40*6];
   }
   Sign_Mass_Syst_PtLow = new TH1D*[signals.size()];
   Sign_Mass_Syst_ILow  = new TH1D*[signals.size()];


   char PredExt[1024];
   char SignExt[1024];
   char DataExt[1024];
   char MCTrExt[1024];
   char Name   [1024];
   for(unsigned int i=0;i<40*6;i++){
      sprintf(PredExt,"Pred");
      GetNameFromIndex(PredExt, i);
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

      sprintf(Name,"Mass_%s",DataExt);
      Data_Mass[i]         = new TH1D(Name,Name,400,0,MassHistoUpperBound);
      Data_Mass[i]->Sumw2();

      sprintf(Name,"Mass_%s",PredExt);
      Pred_Mass[i] = new TH1D(Name,Name,400,0,MassHistoUpperBound);
      Pred_Mass[i]->Sumw2();

      sprintf(Name,"Mass_%s",MCTrExt);
      MCTr_Mass[i] = new TH1D(Name,Name,400,0,MassHistoUpperBound);
      MCTr_Mass[i]->Sumw2();

      for(unsigned int s=0;s<signals.size();s++){
         sprintf(SignExt,"%s",signals[s].Name.c_str());
         GetNameFromIndex(SignExt, i);
         sprintf(Name,"Mass_%s",SignExt);
         Sign_Mass[s][i] = new TH1D(Name,Name,400,0,MassHistoUpperBound);
         Sign_Mass[s][i]->Sumw2();
      }

      sprintf(Name,"I_%s",PredExt);
      Pred_I[i]  = new TH1D(Name,Name,200,0,dEdxUpLim[dEdxMassIndex]);
      Pred_I[i]->Sumw2();

      sprintf(Name,"P_%s",PredExt);
      Pred_P[i]  = new TH1D(Name,Name,400,0,PtHistoUpperBound);
      Pred_P[i]->Sumw2();

      sprintf(Name,"PI_%s",PredExt);
      Pred_PI[i] = new TH2D(Name,Name,400,0,PtHistoUpperBound, 400, 0, dEdxUpLim[dEdxMassIndex]);
      Pred_PI[i]->Sumw2();

      sprintf(Name,"PI_A_%s",DataExt);
      Data_PI_A[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_A[i]->Sumw2();

      sprintf(Name,"PI_B_%s",DataExt);
      Data_PI_B[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_B[i]->Sumw2();

      sprintf(Name,"PI_C_%s",DataExt);
      Data_PI_C[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxUpLim[dEdxMassIndex]);
      Data_PI_C[i]->Sumw2();

      sprintf(Name,"PI_D_%s",DataExt);
      Data_PI_D[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxUpLim[dEdxMassIndex]);
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

   for(unsigned int s=0;s<signals.size();s++){
      sprintf(Name,"Mass_%s_Syst_PtLow", signals[s].Name.c_str());
      Sign_Mass_Syst_PtLow[s] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
      Sign_Mass_Syst_PtLow[s]->Sumw2();

      sprintf(Name,"Mass_%s_Syst_ILow", signals[s].Name.c_str());
      Sign_Mass_Syst_ILow[s] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
      Sign_Mass_Syst_ILow[s]->Sumw2();
   }


}


void CleanUpHistos(){
    //Crashing for unknown reason, is not a high priority problem, so just switch off
/*
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
*/
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

void SetWeightMC(double IntegratedLuminosityInPb, double SampleEquivalentLumi, double SampleSize, double MaxEvent){
   if(MaxEvent<0)MaxEvent=SampleSize;
   printf("SetWeight MC: IntLumi = %6.2E  SampleLumi = %6.2E --> EventWeight = %6.2E\n",IntegratedLuminosityInPb,SampleEquivalentLumi, IntegratedLuminosityInPb/SampleEquivalentLumi);
   printf("Sample NEvent = %6.2E   SampleEventUsed = %6.2E --> Weight Rescale = %6.2E\n",SampleSize, MaxEvent, SampleSize/MaxEvent);
   Event_Weight = (IntegratedLuminosityInPb/SampleEquivalentLumi) * (SampleSize/MaxEvent);
   printf("FinalWeight = %6.2f\n",Event_Weight);
}







int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl){
   int toReturn = 0;
   for(unsigned int g=0;g<genColl.size();g++){
      if(abs(genColl[g].pdgId())<1000000)continue;
      if(genColl[g].status()!=1)continue;

      //skip neutral RHadrons (gluino)
      if(abs(genColl[g].pdgId())==1000993)continue;
      if(abs(genColl[g].pdgId())==1009313)continue;
      if(abs(genColl[g].pdgId())==1009113)continue;
      if(abs(genColl[g].pdgId())==1009223)continue;
      if(abs(genColl[g].pdgId())==1009333)continue;
      if(abs(genColl[g].pdgId())==1092114)continue;
      if(abs(genColl[g].pdgId())==1093214)continue;
      if(abs(genColl[g].pdgId())==1093324)continue;

      //skip neutral RHadrons (stop)
      if(abs(genColl[g].pdgId())==1000622)continue;
      if(abs(genColl[g].pdgId())==1000642)continue;
      if(abs(genColl[g].pdgId())==1006113)continue;
      if(abs(genColl[g].pdgId())==1006311)continue;
      if(abs(genColl[g].pdgId())==1006313)continue;
      if(abs(genColl[g].pdgId())==1006333)continue;


      toReturn++;
   }
   return toReturn;
}



bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold)
{
   unsigned int filterIndex = trEv.filterIndex(InputPath);
   //if(filterIndex<trEv.sizeFilters())printf("SELECTED INDEX =%i --> %s    XXX   %s\n",filterIndex,trEv.filterTag(filterIndex).label().c_str(), trEv.filterTag(filterIndex).process().c_str());

   if (filterIndex<trEv.sizeFilters()){
      const trigger::Vids& VIDS(trEv.filterIds(filterIndex));
      const trigger::Keys& KEYS(trEv.filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      const trigger::TriggerObjectCollection& TOC(trEv.getObjects());


      if(!averageThreshold){
         int NObjectAboveThresholdObserved = 0;
         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            if(TO.pt()> NewThreshold) NObjectAboveThresholdObserved++;
            //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }
         if(NObjectAboveThresholdObserved>=NObjectAboveThreshold)return true;

      }else{
         std::vector<double> ObjPt;

         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            ObjPt.push_back(TO.pt());
            //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }
         if((int)(ObjPt.size())<NObjectAboveThreshold)return false;
         std::sort(ObjPt.begin(), ObjPt.end());

         double Average = 0;
         for(int i=0; i<NObjectAboveThreshold;i++){
            Average+= ObjPt[ObjPt.size()-1-i];
         }Average/=NObjectAboveThreshold;
         //cout << "AVERAGE = " << Average << endl;

         if(Average>NewThreshold)return true;
      }
   }
   return false;
}

