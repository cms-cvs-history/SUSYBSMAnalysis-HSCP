#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDCacheFile.h"
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
#include "TRandom3.h"
#include "TProfile.h"
#include "TDirectory.h"

namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle; class HSCPIsolation; class MuonSegment;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag; class LumiReWeighting;}
namespace reweight{class PoissonMeanShifter;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;
using namespace reweight;

#endif

#include "Analysis_Global.h"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"
#include "Analysis_Samples.h"

/////////////////////////// FUNCTION DECLARATION /////////////////////////////

void Analysis_Step3(char* SavePath);

void InitHistos();

double DistToHSCP      (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest);
int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl);
void  GetGenHSCPBeta   (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, bool onlyCharged=true);
bool   PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st=NULL, const double& GenBeta=-1, bool RescaleP=false, const double& RescaleI=0.0, const double& RescaleT=0.0);
bool   PassSAPreselection(const susybsm::HSCParticle& hscp, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, int& DzType, bool Control=false, const double& GenBeta=-1, const double& GenPt=-1, const double& GenCharge=-1);
bool PassSelection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const fwlite::ChainEvent& ev, const int& CutIndex=0, stPlots* st=NULL, const double& GenBeta=-1, bool RescaleP=false, const double& RescaleI=0.0, const double& RescaleT=0.0);
bool PassTrigger      (const fwlite::ChainEvent& ev);
bool hasGoodPtHat     (const fwlite::ChainEvent& ev, const double& PtMax);
double SegSep(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double& minPhi, double& minEta);

double GetPUWeight(const fwlite::ChainEvent& ev, const bool& Iss4pileup, double &PUSystFactor);
double GetSampleWeight(const double& IntegratedLuminosityInPb=-1, const double& IntegratedLuminosityInPbBeforeTriggerChange=-1, const double& CrossSection=0, const double& MCEvents=0, int period=0);
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent);
double RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge);
unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames);
/////////////////////////// VARIABLE DECLARATION /////////////////////////////

class DuplicatesClass{
   private :
      typedef std::map<std::pair<unsigned int, unsigned int>, bool > RunEventHashMap;
      RunEventHashMap map;
   public :
        DuplicatesClass(){}
        ~DuplicatesClass(){}
        void Clear(){map.clear();}
        bool isDuplicate(unsigned int Run, unsigned int Event){
	   RunEventHashMap::iterator it = map.find(std::make_pair(Run,Event));
           if(it==map.end()){
   	      map[std::make_pair(Run,Event)] = true;
              return false;
           }
           return true;
        }
};


TFile* HistoFile;
/*
TH1D*  Hist_Pt ;
TH1D*  Hist_Is  ;
TH1D*  Hist_TOF;

TH2D*  Pred_Mass;
TH2D*  Pred_MassTOF;
TH2D*  Pred_MassComb;

TH1D* H_A;
TH1D* H_B;
TH1D* H_C;
TH1D* H_D;
TH1D* H_E;
TH1D* H_F;
TH1D* H_G;
TH1D* H_H;
TH1D* H_P;
*/
std::vector<double>  CutPt ;
std::vector<double>  CutI  ;
std::vector<double>  CutTOF;

TH1D*  HCuts_Pt;
TH1D*  HCuts_I;
TH1D*  HCuts_TOF;
/*
TH3D*  Pred_EtaP ;
TH2D*  Pred_I    ;
TH2D*  Pred_TOF  ;
TH2D*  Pred_EtaB;
TH2D*  Pred_EtaS;
TH2D*  Pred_EtaS2;

TH2D*  RegionD_P   ;
TH2D*  RegionD_I   ;
TH2D*  RegionD_TOF  ;
*/
std::vector<stSignal> signals;
std::vector<stMC>     MCsample;
std::vector<string>   DataFileName;

stPlots              DataPlots;  
stPlots              DataPlotsControl;
std::vector<stPlots> SignPlots; 
std::vector<stPlots> MCPlots;  
stPlots              MCTrPlots;
//for initializing PileUpReweighting utility.
const   float TrueDist2011_f[35] = {0.00285942, 0.0125603, 0.0299631, 0.051313, 0.0709713, 0.0847864, 0.0914627, 0.0919255, 0.0879994, 0.0814127, 0.0733995, 0.0647191, 0.0558327, 0.0470663, 0.0386988, 0.0309811, 0.0241175, 0.018241, 0.0133997, 0.00956071, 0.00662814, 0.00446735, 0.00292946, 0.00187057, 0.00116414, 0.000706805, 0.000419059, 0.000242856, 0.0001377, 7.64582e-05, 4.16101e-05, 2.22135e-05, 1.16416e-05, 5.9937e-06, 5.95542e-06};//from 2011 Full dataset

const   float Pileup_MC[35]= {1.45346E-01, 6.42802E-02, 6.95255E-02, 6.96747E-02, 6.92955E-02, 6.84997E-02, 6.69528E-02, 6.45515E-02, 6.09865E-02, 5.63323E-02, 5.07322E-02, 4.44681E-02, 3.79205E-02, 3.15131E-02, 2.54220E-02, 2.00184E-02, 1.53776E-02, 1.15387E-02, 8.47608E-03, 6.08715E-03, 4.28255E-03, 2.97185E-03, 2.01918E-03, 1.34490E-03, 8.81587E-04, 5.69954E-04, 3.61493E-04, 2.28692E-04, 1.40791E-04, 8.44606E-05, 5.10204E-05, 3.07802E-05, 1.81401E-05, 1.00201E-05, 5.80004E-06};

edm::LumiReWeighting LumiWeightsMC_;
std::vector< float > BgLumiMC; //MC                                           
std::vector< float > TrueDist2011;                                    
reweight::PoissonMeanShifter PShift_(0.6);//0.6 for upshift, -0.6 for downshift
/////////////////////////// CODE PARAMETERS /////////////////////////////

void Analysis_Step23(string MODE_="COMPILE", int TypeMode_=0, string dEdxSel_="dedxASmi", string dEdxMass_="dedxHarm2", string TOF_Label_="combined", double CutPt_=-1.0, double CutI_=-1, double CutTOF_=-1, float MinPt_=GlobalMinPt, float MaxEta_=GlobalMaxEta)
{
  MODE=MODE_;
   if(MODE=="COMPILE")return;

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
   TH1::AddDirectory(kTRUE);

   bool TkOnly=true;
   if(TypeMode_==2)  TkOnly=false;
   GetSignalDefinition(signals, TkOnly);
   GetMCDefinition(MCsample);

   char Buffer[2048];   
   char Command[2048];
   DataFileName.clear();
   GetInputFiles(DataFileName, "Data");

   dEdxS_Label = dEdxSel_;
   dEdxM_Label = dEdxMass_;
   TOF_Label   = TOF_Label_;
   InitdEdx(dEdxS_Label);

   TypeMode  = TypeMode_;
   GlobalMaxEta = MaxEta_;
   GlobalMinPt    = MinPt_;

   //for 2012 running with dEdx triggers
   if(GlobalMinPt>=50){GlobalMinIm   =   3.5;}

   if(TypeMode==0){
      GlobalMinNDOF   = 0; 
      GlobalMinTOF    = 0;
   }else{
      GlobalMaxTIsol *= 2;
      GlobalMaxEIsol *= 2;
   }

   if(TypeMode!=3) {CutPt .push_back(GlobalMinPt);   CutI  .push_back(GlobalMinIs);  CutTOF.push_back(GlobalMinTOF);}
   else CutPt .push_back(SAMinPt);   CutI  .push_back(-1);  CutTOF.push_back(GlobalMinTOF);

   if(TypeMode==0){   
      for(double Pt =GlobalMinPt+5 ; Pt <200;Pt+=5){
      for(double I  =GlobalMinIs+0.025  ; I  <0.45 ;I+=0.025){
         CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(-1);
      }}
   }else if(TypeMode==2){
      for(double Pt =GlobalMinPt+5 ; Pt <120;  Pt+=5){
      for(double I  =GlobalMinIs +0.025; I  <0.40;  I+=0.025){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.35;TOF+=0.025){
         CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(TOF);
      }}}
   }else if(TypeMode==3){
     for(double Pt =SAMinPt+30 ; Pt <550;  Pt+=30){
       for(double TOF=GlobalMinTOF+0.01; TOF<1.4;TOF+=0.01){
	 CutPt .push_back(Pt);   CutI  .push_back(-1);  CutTOF.push_back(TOF);
       }}
   }
   printf("%i Different Final Selection will be tested\n",(int)CutPt.size());

   //initialize LumiReWeighting
   for(int i=0; i<35; ++i)   BgLumiMC.push_back(Pileup_MC[i]);
   for(int i=0; i<35; ++i)    TrueDist2011.push_back(TrueDist2011_f[i]);
   LumiWeightsMC_ = edm::LumiReWeighting(BgLumiMC, TrueDist2011);

   sprintf(Buffer,"Results/"       );                                          sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%s%s/"         ,Buffer,dEdxS_Label.c_str());                sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%s%s/"         ,Buffer,TOF_Label.c_str());                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sEta%02.0f/"  ,Buffer,10.0*GlobalMaxEta);                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sPtMin%02.0f/",Buffer,GlobalMinPt);                        sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sType%i/"     ,Buffer,TypeMode);                           sprintf(Command,"mkdir %s",Buffer); system(Command);
//   sprintf(Buffer,"%sPt%03.0f/"   ,Buffer,CutPt [0]);		               sprintf(Command,"mkdir %s",Buffer); system(Command);
//   sprintf(Buffer,"%sI%05.2f/"    ,Buffer,CutI  [0]);                          sprintf(Command,"mkdir %s",Buffer); system(Command);
//   sprintf(Buffer,"%sTOF%05.2f/"  ,Buffer,CutTOF[0]);                          sprintf(Command,"mkdir %s",Buffer); system(Command);

   time_t start = time(NULL);
   if(MODE=="ANALYSE_DATA"){
      signals.clear();  //Remove all signal samples
      MCsample.clear();
      HistoFile = new TFile((string(Buffer) + "/Histos_Data.root").c_str(),"RECREATE");
   }else if(MODE=="ANALYSE_SIGNAL"){
      DataFileName.clear();  //Remove all data files
      MCsample.clear();
      HistoFile = new TFile((string(Buffer) + "/Histos.root").c_str(),"RECREATE");
   }else if(MODE=="ANALYSE_MC"){
      DataFileName.clear();  //Remove all data files
      signals.clear();  //Remove all signal samples
      HistoFile = new TFile((string(Buffer) + "/Histos_MC.root").c_str(),"RECREATE");
   }else if(MODE=="ANALYSE_COSMIC"){
     signals.clear();  //Remove all signal samples                                                                                                                                  
     MCsample.clear();
     HistoFile = new TFile((string(Buffer) + "/Histos_Cosmic.root").c_str(),"RECREATE");
   }else{
      printf("You must select a MODE:\n");
      printf("MODE='ANALYSE_DATA'   : Will run the analysis on Data\n"); 
      printf("MODE='ANALYSE_SIGNAL' : Will run the analysis on Signal MC\n");
      printf("MODE='ANALYSE_MC'     : Will run the analysis on Background MC\n");
      return;
   }
   InitHistos();
   Analysis_Step3(Buffer);
   HistoFile->Write();
   HistoFile->Close();
   time_t end = time(NULL);
   printf("RUN TIME = %i sec\n",(int)(end-start));
   return;
}

bool hasGoodPtHat(const fwlite::ChainEvent& ev, const double& PtMax){
   if(PtMax<0)return true;
   fwlite::Handle< GenEventInfoProduct > genInfo;
   genInfo.getByLabel(ev, "generator");
   if(!genInfo.isValid()){printf("genInfo NotFound\n");return false;}
   if((genInfo->binningValues()[0])<PtMax)return true;
   return false;
}


#ifdef ANALYSIS2011
bool PassTrigger(const fwlite::ChainEvent& ev)
{
      edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");
      if(!tr.isValid())return false;

      if(TypeMode!=3) {
      if(tr.accept(tr.triggerIndex("HscpPathSingleMu")))return true;
//      if(tr.accept(tr.triggerIndex("HscpPathDoubleMu")))return true;
      else if(tr.accept(tr.triggerIndex("HscpPathPFMet"))) {
	if(MODE!="ANALYSE_DATA") Event_Weight=Event_Weight*0.96;
	return true;
      }}

      else {
      fwlite::Handle<reco::PFMETCollection> pfMETCollection;
      pfMETCollection.getByLabel(ev,"pfMet");
      if(tr.accept(tr.triggerIndex("HSCPPathSAMU")) && pfMETCollection->begin()->et()>60) return true;}

//      if(tr.accept(tr.triggerIndex("HscpPathCaloMet")))return true;
      return false;
}
#else
bool PassTrigger(const fwlite::ChainEvent& ev)
{
    edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");
    if(!tr.isValid())return false;

    if(TypeMode!=3) {
      if(tr.accept("HSCPHLTTriggerMetDeDxFilter"))return true;
      if(tr.accept("HSCPHLTTriggerMuDeDxFilter"))return true;
      if(tr.accept("HSCPHLTTriggerMuFilter"))return true;
    }
    else {
      if(MODE!="ANALYSE_COSMIC") {if(tr.accept(tr.triggerIndex("HSCPHLTTriggerL2MuFilter"))) return true;}
      else {if(tr.accept(tr.triggerIndex("HSCPHLTTriggerCosmicFilter"))) return true;}
    }
    if(tr.accept(tr.triggerIndex("HSCPHLTTriggerL2MuFilter"))) return true;

      return false;
}
#endif

bool PassSAPreselection(const susybsm::HSCParticle& hscp, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, int& DzType, bool Control, const double& GenBeta, const double& GenPt, const double &GenCharge)
{   

  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return false;  

   reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull())return false;
   reco::TrackRef innertrack = hscp.trackRef();

   bool isGlobal=(muon->isGlobalMuon() && muon->isTrackerMuon());

   if(!tof) return false;

   //Make distributions without any cuts
   if(st) {
     st->BS_Pt_All->Fill(track->pt(), Event_Weight);
     if(tof->nDof()>0)st->BS_TOF_All->Fill(tof->inverseBeta(), Event_Weight);
     st->Total->Fill(0.0,Event_Weight);
     if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);
   }

   //Require track to match trigger object
   //st->DistTrigger->Fill(DistToTrigger(hscp, ev),Event_Weight);
   //if(DistToTrigger(hscp, ev)>MaxDistTrigger) return false;
   st->TriggerMatch->Fill(0.0, Event_Weight);

   //Match to a SA track without vertex constraint for IP cuts
   fwlite::Handle< std::vector<reco::Track> > noVertexTrackCollHandle;
   noVertexTrackCollHandle.getByLabel(ev,"RefitSAMuons", "");
   //noVertexTrackCollHandle.getByLabel(ev,"RefitMTSAMuons", "");
   if(!noVertexTrackCollHandle.isValid()){
     noVertexTrackCollHandle.getByLabel(ev,"refittedStandAloneMuons", "");
     if(!noVertexTrackCollHandle.isValid()){
       noVertexTrackCollHandle.getByLabel(ev,"RefitMTSAMuons", "");
       if(!noVertexTrackCollHandle.isValid()){printf("No Vertex Track Collection Not Found\n");return false;}
     }
   }

   const std::vector<reco::Track>& noVertexTrackColl = *noVertexTrackCollHandle;

   reco::Track NVTrack;
   double minDr=15;
   for(unsigned int i=0;i<noVertexTrackColl.size();i++){
     double dEta = track->eta()-noVertexTrackColl[i].eta();
     double dPhi = track->phi()-noVertexTrackColl[i].phi();
     if(dPhi>3.14159) dPhi=6.29-dPhi;
     double dR=sqrt(dEta*dEta+dPhi*dPhi);
     if(dR<minDr) {minDr=dR;
       NVTrack=noVertexTrackColl[i];}
   }

   if(st) st->BS_dR_NVTrack->Fill(minDr,Event_Weight);
   if(minDr>0.4) return false;
   st->NVTrack->Fill(0.0,Event_Weight);
   //Cut on max eta
   if(fabs(track->eta())>SAMaxEta) return false;
   if(st) st->Eta->Fill(0.0, Event_Weight);

   //Cut on min pt
   if(track->pt()<SAMinPt)return false;
   if(st) st->MinPt->Fill(0.0, Event_Weight);

   //Cut on number of matched muon stations
   int count=track->hitPattern().muonStationsWithValidHits();
   if(st) {
     st->BS_MatchedStations->Fill(count, Event_Weight);  ;
   }
   if(count<2) return false;
   if(st) st->Stations->Fill(0.0, Event_Weight);

   //Cut on number of TOF degrees of freedom
   if(st){st->BS_nDof->Fill(tof->nDof(),Event_Weight);}
   if(tof->nDof()<GlobalMinNDOF) return false;
   if(st)st->nDof  ->Fill(0.0,Event_Weight);

   //Cut on inverse beta error
   if(st){st->BS_InvBetaErr->Fill(tof->inverseBetaErr(), Event_Weight);}
   if(tof->inverseBetaErr()>GlobalMaxTOFErr)return false;
   if(st) st->tofError->Fill(0.0, Event_Weight);

   //Plots on number hits overall, barrel only, and endcap
   if(st) {st->BS_TNOH->Fill(track->hitPattern().numberOfValidMuonHits(),Event_Weight);
     if(fabs(track->eta())<CSCRegion) st->BS_TNOH_Barrel->Fill(track->hitPattern().numberOfValidMuonHits(),Event_Weight);
     if(fabs(track->eta())>CSCRegion) st->BS_TNOH_Endcap->Fill(track->hitPattern().numberOfValidMuonHits(),Event_Weight);
   }

   //Get vertex collection and require at least one good primary vertex, not required for cosmics
   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
   const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

   bool goodVertex=false;
   int goodVerts=0;
   for(unsigned int i=0;i<vertexColl.size();i++){
     if(fabs(vertexColl[i].z())<15 && sqrt(vertexColl[i].x()*vertexColl[i].x()+vertexColl[i].y()*vertexColl[i].y())<2 && vertexColl[i].ndof()>3) {goodVertex=true; goodVerts++;}
   }

   if(st) st->BS_PV->Fill(goodVerts,Event_Weight);

   if(MODE!="ANALYSE_COSMIC" && !goodVertex) return false;

   //Find displacement of tracks with respect to beam spot
   fwlite::Handle<reco::BeamSpot> beamSpotCollHandle;
   beamSpotCollHandle.getByLabel(ev,"offlineBeamSpot");
   if(!beamSpotCollHandle.isValid()){printf("Beam Spot Collection NotFound\n");return false;}
   const reco::BeamSpot& beamSpotColl = *beamSpotCollHandle;

   double dz  = NVTrack.dz (beamSpotColl.position());
   double dxy = NVTrack.dxy(beamSpotColl.position());
   double v3d = sqrt(dz*dz+dxy*dxy);

   //Find number of DT segments with a ZED projection, returns 6 if any CSC hits are used in track
   //int ZedSegs=6;
   //ZedSegs=Zed(hscp, ev);

   //Find distance to nearest segment on opposite side of detector
   double minPhi=10, minEta=10;
   double segSep=SegSep(hscp, ev, minPhi, minEta);

   if(st){st->BS_SegSep->Fill(segSep, Event_Weight);
     if(fabs(track->eta())>CSCRegion) st->BS_SegMinEtaSep_CSC->Fill(minEta, Event_Weight);
     else if(fabs(track->eta())<DTRegion) st->BS_SegMinEtaSep_DT->Fill(minEta, Event_Weight);
     st->BS_SegMinPhiSep->Fill(minPhi, Event_Weight);
     st->BS_SegMinEtaSep->Fill(minEta, Event_Weight);
     st->BS_SegEta_MinEtaSep->Fill(track->eta(), minEta, Event_Weight);

     //Plots for tracking failing Eta Sep cut
     if(fabs(minEta)<minSegEtaSep) {
       st->BS_Dz_FailPhi->Fill(dz);
       st->BS_Pt_FailPhi->Fill(track->pt(), Event_Weight);
       st->BS_TOF_FailPhi->Fill(tof->inverseBeta(), Event_Weight);
     }
     //Plots for tracks passing eta separation cut
     else {
       st->BS_Dz_PassPhi->Fill(dz);
       st->BS_Pt_PassPhi->Fill(track->pt(), Event_Weight);
       st->BS_TOF_PassPhi->Fill(tof->inverseBeta(), Event_Weight);
     }

     //Plotting segment separation depending on whether track passed dz cut
     if(fabs(dz)>SAMaxDz) {
       if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) {
         st->BS_SegMinEtaSep_FailDz->Fill(minEta, Event_Weight);
       }
     }
     else {
       st->BS_SegMinEtaSep_PassDz->Fill(minEta, Event_Weight);
     }
   }

   //Now cut Eta separation
   if(fabs(minEta)<minSegEtaSep) return false;
   if(st){st->SegSep->Fill(0.0,Event_Weight);}

   if(st){
     //Various plots involving dxy
     st->BS_V3D->Fill(v3d,Event_Weight);
     st->BS_Dxy->Fill(dxy,Event_Weight);
     st->BS_Dxy_Dz->Fill(dxy, dz, Event_Weight);

     //Plots for tracks passing or failing dxy cut
     if(fabs(dxy)>10) {
       st->BS_Pt_FailDxy->Fill(track->pt(), Event_Weight);
       st->BS_TOF_FailDxy->Fill(tof->inverseBeta(), Event_Weight);
       st->BS_Eta_FailDxy->Fill(track->eta(), Event_Weight);
       if(isGlobal)st->BS_GlobalPt_FailDxy->Fill(muon->pt(), Event_Weight);
       st->BS_Dz_FailDxy->Fill(dz,Event_Weight);
     }
     else {
       st->BS_Pt_PassDxy->Fill(track->pt(), Event_Weight);
       st->BS_TOF_PassDxy->Fill(tof->inverseBeta(), Event_Weight);
       st->BS_Eta_PassDxy->Fill(track->eta(), Event_Weight);
       if(isGlobal)st->BS_GlobalPt_PassDxy->Fill(muon->pt(), Event_Weight);
       st->BS_Dz_PassDxy->Fill(dz,Event_Weight);
     }
   }
   //Cut on dxy only for non global tracks
   if(fabs(dxy)>SAMaxDxy) return false;
   if(st) st->Dxy->Fill(0.0,Event_Weight);

   if(st) {
     //Plots for tracks in dz control region
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) {
       st->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);
       st->BS_TOF_FailDz->Fill(tof->inverseBeta(), Event_Weight);
       st->BS_Time_FailDz->Fill(tof->timeAtIpInOut(), Event_Weight);
       st->BS_Eta_FailDz->Fill(track->eta(), Event_Weight);
       if(fabs(track->eta())>CSCRegion) {
	 st->BS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Time_FailDz_CSC->Fill(tof->timeAtIpInOut(), Event_Weight);
	 st->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);
       }
       else if(fabs(track->eta())<DTRegion) {
	 st->BS_TOF_FailDz_DT->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Time_FailDz_DT->Fill(tof->timeAtIpInOut(), Event_Weight);
	 st->BS_Pt_FailDz_DT->Fill(track->pt(), Event_Weight);
       }
     }
     //Plots for tracks with relatively small dz
     else if(fabs(dz)<25){
       st->BS_Pt_PassDz->Fill(track->pt(), Event_Weight);
       st->BS_TOF_PassDz->Fill(tof->inverseBeta(), Event_Weight);
       st->BS_Time_PassDz->Fill(tof->timeAtIpInOut(), Event_Weight);
       if(fabs(track->eta())>CSCRegion) {
	 st->BS_TOF_PassDz_CSC->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Time_PassDz_CSC->Fill(tof->timeAtIpInOut(), Event_Weight);
         st->BS_Pt_PassDz_CSC->Fill(track->pt(), Event_Weight);
       }
       else if(fabs(track->eta())<DTRegion) {
	 st->BS_TOF_PassDz_DT->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Time_PassDz_DT->Fill(tof->timeAtIpInOut(), Event_Weight);
	 st->BS_Pt_PassDz_DT->Fill(track->pt(), Event_Weight);
       }
       st->BS_Eta_PassDz->Fill(track->eta(), Event_Weight);
     }

     //Plots of dz
     st->BS_Dz->Fill(dz,Event_Weight);
     if(fabs(track->eta())>CSCRegion) st->BS_Dz_CSC->Fill(dz,Event_Weight);
     else if(fabs(track->eta())<DTRegion) st->BS_Dz_DT->Fill(dz,Event_Weight);
     st->BS_EtaDz->Fill(track->eta(),dz,Event_Weight);
     st->BS_DzTime->Fill(dz,tof->timeAtIpInOut(),Event_Weight);
     if(fabs(track->eta())<DTRegion) st->BS_DzTime_DT->Fill(dz,tof->timeAtIpInOut(),Event_Weight);
     if(fabs(track->eta())>CSCRegion) st->BS_DzTime_CSC->Fill(dz,tof->timeAtIpInOut(),Event_Weight); 
     st->BS_DzPt->Fill(dz,track->pt(),Event_Weight);
   }

   //Count number of SA only tracks falling in dz control region
   if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz && (!isGlobal || MODE=="ANALYSE_DATA")) {
     st->FailDz->Fill(0.0,Event_Weight);
     if(fabs(track->eta())<DTRegion) st->FailDz_DT->Fill(0.0,Event_Weight);
     else st->FailDz_CSC->Fill(0.0,Event_Weight);
   }

   //Split into different dz regions, each different region used to predict cosmic background and find systematic
   if(Control && !isGlobal) {
     if(fabs(dz)<SAMaxDz) DzType=0;
     else if(fabs(dz)<30) DzType=1;
     else if(fabs(dz)<50) DzType=2;
     else if(fabs(dz)<70) DzType=3;
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) DzType=4;
     if(fabs(dz)>CosmicMaxDz) DzType=5;
     //Count number of tracks in dz sidebands passing the TOF cut
     //The pt cut is not applied to increase statistics
     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
       if(tof->inverseBeta()>=CutTOF[CutIndex]) {
	 st->H_D_DzSidebands->Fill(CutIndex, DzType);
	 if(fabs(track->eta())<DTRegion) st->H_D_DzSidebands_DT->Fill(CutIndex, DzType);
         else st->H_D_DzSidebands_CSC->Fill(CutIndex, DzType);
       }
     }
   }

   //Cut on dz for SA only tracks but not if this for the control region
   if(fabs(dz)>SAMaxDz && !Control) return false;

   //Require control region cuts
   if(Control && (fabs(dz)<CosmicMinDz || fabs(dz)>CosmicMaxDz)) return false;

   //Fill a bunch of plots for tracks passing preselection
   if(st){
     st->Dz  ->Fill(0.0,Event_Weight);
     //if(ZedSegs==0) st->BS_Dz_NoZed->Fill(dz,Event_Weight);
     st->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);
     st->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);
     st->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);
     st->BS_EtaTOF->Fill(track->eta(),tof->inverseBeta(),Event_Weight);
     st->BS_EtaTime->Fill(track->eta(),tof->timeAtIpInOut(),Event_Weight);
     st->BS_Eta->Fill(track->eta(),Event_Weight);
     st->BS_Phi->Fill(track->phi(),Event_Weight);
     st->BS_PhiTime->Fill(track->phi(),tof->timeAtIpInOut(),Event_Weight);
     st->BS_PterrSq ->Fill(track->ptError()/(track->pt()*track->pt()),Event_Weight);
     st->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);

     if(GenPt>0 && GenCharge!=0) {
       st->BS_QoverPt->Fill(1-(track->charge()*GenPt)/(GenCharge*track->pt()),Event_Weight);
     }

     if(GenBeta>=0)st->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
     if(GenBeta>=0)st->BS_BetaDiff->Fill(((1./GenBeta)-tof->inverseBeta())/(1./GenBeta));
     st->Preselected ->Fill(0.0,Event_Weight);
     if(fabs(track->eta())<DTRegion) st->Preselected_DT->Fill(0.0,Event_Weight);
     else st->Preselected_CSC->Fill(0.0,Event_Weight);
     st->BS_P  ->Fill(track->p(),Event_Weight);
     st->BS_Pt ->Fill(track->pt(),Event_Weight);
     st->BS_PtTOF->Fill(track->pt(),tof->inverseBeta(), Event_Weight);
     st->BS_TOF->Fill(tof->inverseBeta(),Event_Weight);

     if(fabs(track->eta())<0.9) {
       st->BS_TOF_Bar->Fill(tof->inverseBeta(),Event_Weight);
       st->BS_Pt_Bar->Fill(track->pt(),Event_Weight);
     }
     else {
       st->BS_TOF_For->Fill(tof->inverseBeta(),Event_Weight);
       st->BS_Pt_For->Fill(track->pt(),Event_Weight);
     }

     if(dttof->nDof()>6) st->BS_TOF_DT->Fill(dttof->inverseBeta(),Event_Weight);
     if(csctof->nDof()>6) st->BS_TOF_CSC->Fill(csctof->inverseBeta(),Event_Weight);

     st->BS_VertexTime->Fill(tof->timeAtIpInOut(),Event_Weight);
     if(innertrack.isNull()) st->BS_IsTracker->Fill(0.0, Event_Weight);
     else st->BS_IsTracker->Fill(1.0, Event_Weight);
     if(!muon->combinedQuality().updatedSta) st->BS_IsUpdated->Fill(0.0, Event_Weight);
     else st->BS_IsUpdated->Fill(1.0, Event_Weight);
     if(isGlobal) st->BS_Pt_Global->Fill(muon->pt(),Event_Weight);
   }
   return true;
}

bool PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT)
{

   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
   if(TypeMode==2 && hscp.type() != HSCParticleType::globalMuon)return false;
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

   if(st){st->Total->Fill(0.0,Event_Weight);
     if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);
     st->BS_Eta->Fill(track->eta(),Event_Weight);
   }

   if(fabs(track->eta())>GlobalMaxEta) return false;

   if(st){st->BS_TNOH->Fill(track->found(),Event_Weight);
          st->BS_TNOHFraction->Fill(track->validFraction(),Event_Weight);
   }

   if(track->found()<GlobalMinNOH)return false;
   if(track->validFraction()<0.80)return false;
   if(track->hitPattern().numberOfValidPixelHits()<2)return false;

   if(st){st->TNOH  ->Fill(0.0,Event_Weight);
          st->BS_TNOM->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
   }
   if(dedxSObj->numberOfMeasurements()<GlobalMinNOM)return false;
   if(st){st->TNOM  ->Fill(0.0,Event_Weight);}

   if(tof){
   if(st){st->BS_nDof->Fill(tof->nDof(),Event_Weight);}
   if(TypeMode==2 && tof->nDof()<GlobalMinNDOF && (dttof->nDof()<GlobalMinNDOFDT || csctof->nDof()<GlobalMinNDOFCSC) )return false;
   }

   if(st){st->nDof  ->Fill(0.0,Event_Weight);
          st->BS_Qual->Fill(track->qualityMask(),Event_Weight);
   }

   if(track->qualityMask()<GlobalMinQual )return false;
   if(st){st->Qual  ->Fill(0.0,Event_Weight);
          st->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);
   }
   if(track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(st){st->Chi2  ->Fill(0.0,Event_Weight);}

   if(st && GenBeta>=0)st->Beta_PreselectedA->Fill(GenBeta, Event_Weight);

   if(st){st->BS_MPt ->Fill(track->pt(),Event_Weight);}
   if(RescaleP){ if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<GlobalMinPt)return false;
   }else{        if(track->pt()<GlobalMinPt)return false;   }

   if(st){st->MPt   ->Fill(0.0,Event_Weight);
          st->BS_MIs->Fill(dedxSObj->dEdx(),Event_Weight);
          st->BS_MIm->Fill(dedxMObj->dEdx(),Event_Weight);
   }
   if(dedxSObj->dEdx()+RescaleI<GlobalMinIs)return false;
   if(dedxMObj->dEdx()<GlobalMinIm)return false;
   if(st){st->MI   ->Fill(0.0,Event_Weight);}
   if(tof){
   if(st){st->BS_MTOF ->Fill(tof->inverseBeta(),Event_Weight);}
   //if(TypeMode==2 && tof->inverseBeta()+RescaleT<GlobalMinTOF)return false;
   if(TypeMode==2 && tof->inverseBetaErr()>GlobalMaxTOFErr)return false;
   }
   if(st){st->MTOF ->Fill(0.0,Event_Weight);
      if(GenBeta>=0)st->Beta_PreselectedB->Fill(GenBeta, Event_Weight);
   }

   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
   const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

   double dz  = track->dz (vertexColl[0].position());
   double dxy = track->dxy(vertexColl[0].position());
   for(unsigned int i=1;i<vertexColl.size();i++){
      if(fabs(track->dz (vertexColl[i].position())) < fabs(dz) ){
         dz  = track->dz (vertexColl[i].position());
         dxy = track->dxy(vertexColl[i].position());
      }
   }
   double v3d = sqrt(dz*dz+dxy*dxy);

   if(st){st->BS_V3D->Fill(v3d,Event_Weight);}
   if(v3d>GlobalMaxV3D )return false;
   if(st){st->V3D  ->Fill(0.0,Event_Weight);}

   fwlite::Handle<HSCPIsolationValueMap> IsolationH;
   IsolationH.getByLabel(ev, "HSCPIsolation03");
   if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
   const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();

   HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());

   if(st){st->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
    if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;
   if(st){st->TIsol   ->Fill(0.0,Event_Weight);}

   double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
   if(st){st->BS_EIsol ->Fill(EoP,Event_Weight);}
   if(EoP>GlobalMaxEIsol)return false;
   if(st){st->EIsol   ->Fill(0.0,Event_Weight);}

   if(st){st->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if((track->ptError()/track->pt())>GlobalMaxPterr)return false;

   if(std::max(0.0,track->pt())<GlobalMinPt)return false;
   if(st){st->Pterr   ->Fill(0.0,Event_Weight);}

   if(st){st->BS_EtaIs->Fill(track->eta(),dedxSObj->dEdx(),Event_Weight);
          st->BS_EtaIm->Fill(track->eta(),dedxMObj->dEdx(),Event_Weight);
          st->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);
          st->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);
          if(tof)st->BS_EtaTOF->Fill(track->eta(),tof->inverseBeta(),Event_Weight);
   }
   if(fabs(track->eta())>GlobalMaxEta) return false;

   if(st){if(GenBeta>=0)st->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
          st->BS_P  ->Fill(track->p(),Event_Weight);
          st->BS_Pt ->Fill(track->pt(),Event_Weight);
          st->BS_Is ->Fill(dedxSObj->dEdx(),Event_Weight);
          st->BS_Im ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(tof) {
	    st->BS_TOF->Fill(tof->inverseBeta(),Event_Weight);
	    if(dttof->nDof()>6) st->BS_TOF_DT->Fill(dttof->inverseBeta(),Event_Weight);
            if(csctof->nDof()>6) st->BS_TOF_CSC->Fill(csctof->inverseBeta(),Event_Weight);
	  }
          st->BS_PIs  ->Fill(track->p()  ,dedxSObj->dEdx(),Event_Weight);
          st->BS_PIm  ->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
          st->BS_PtIs ->Fill(track->pt() ,dedxSObj->dEdx(),Event_Weight);
          st->BS_PtIm ->Fill(track->pt() ,dedxMObj->dEdx(),Event_Weight);
          if(tof)st->BS_TOFIs->Fill(tof->inverseBeta(),dedxSObj->dEdx(),Event_Weight);
          if(tof)st->BS_TOFIm->Fill(tof->inverseBeta(),dedxMObj->dEdx(),Event_Weight);
   }

   return true;
}

bool PassSelection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const fwlite::ChainEvent& ev, const int& CutIndex, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT){
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

   double MuonTOF = GlobalMinTOF;

   if(tof){
      MuonTOF = tof->inverseBeta();
   }

   double Is=0;
   if(dedxSObj) Is=dedxSObj->dEdx();
   double Ih=0;
   if(dedxMObj) Ih=dedxMObj->dEdx();

   if(RescaleP)
   {
     if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<CutPt[CutIndex])return false;
     //if(std::max(0.0,RescaledPt(track->pt() - track->ptError(),track->eta(),track->phi(),track->charge()))<CutPt[CutIndex])return false;
   }
   else
   {
     if(track->pt()<CutPt[CutIndex])return false;
     //if(std::max(0.0,(track->pt() - track->ptError()))<CutPt[CutIndex])return false;
   } 
   if(st){st->Pt    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedP->Fill(CutIndex,GenBeta, Event_Weight);
   }

   if(Is+RescaleI<CutI[CutIndex])return false;
   if(st){st->I    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedI->Fill(CutIndex, GenBeta, Event_Weight);
   }

   if(TypeMode==2 && MuonTOF+RescaleT<CutTOF[CutIndex])return false;
   if(st){st->TOF  ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedT->Fill(CutIndex, GenBeta, Event_Weight);
          st->AS_P  ->Fill(CutIndex,track->p(),Event_Weight);
          st->AS_Pt ->Fill(CutIndex,track->pt(),Event_Weight);
          st->AS_Is ->Fill(CutIndex,Is,Event_Weight);
          st->AS_Im ->Fill(CutIndex,Ih,Event_Weight);
          st->AS_TOF->Fill(CutIndex,MuonTOF,Event_Weight);
//        st->AS_EtaIs->Fill(CutIndex,track->eta(),Is,Event_Weight);
//        st->AS_EtaIm->Fill(CutIndex,track->eta(),Ih,Event_Weight);
//        st->AS_EtaP ->Fill(CutIndex,track->eta(),track->p(),Event_Weight);
//        st->AS_EtaPt->Fill(CutIndex,track->eta(),track->pt(),Event_Weight);
          st->AS_PIs  ->Fill(CutIndex,track->p()  ,Is,Event_Weight);
          st->AS_PIm  ->Fill(CutIndex,track->p()  ,Ih,Event_Weight);
          st->AS_PtIs ->Fill(CutIndex,track->pt() ,Is,Event_Weight);
          st->AS_PtIm ->Fill(CutIndex,track->pt() ,Ih,Event_Weight);
          st->AS_TOFIs->Fill(CutIndex,MuonTOF     ,Is,Event_Weight);
          st->AS_TOFIm->Fill(CutIndex,MuonTOF     ,Ih,Event_Weight);
   }

   return true;
}

void Analysis_FillControlAndPredictionHist(const susybsm::HSCParticle& hscp, const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, stPlots* st=NULL, int DzType=-1){
         reco::TrackRef   track;
         if(TypeMode!=3) track = hscp.trackRef();
	 else {
	   reco::MuonRef muon = hscp.muonRef();
	   if(muon.isNull()) return;
	   track = muon->standAloneMuon();
	 }
	 if(track.isNull())return;

         double MuonTOF = GlobalMinTOF;
         if(tof){MuonTOF = tof->inverseBeta(); }

	 double Is=0;
	 if(dedxSObj) Is=dedxSObj->dEdx();
	 double Ih=0;
	 if(dedxMObj) Ih=dedxMObj->dEdx();

	 st->Hist_Pt->Fill(track->pt(),Event_Weight);
         st->Hist_Is->Fill(Is,Event_Weight);
         st->Hist_TOF->Fill(MuonTOF,Event_Weight);

//          /\ I
//       /\  |----------------------------
//        |  |   |           |             |
//        |  |   |           |             |
//        |  |   |    B      |     D       |
//        |  |   |           |             |
//        |  ------------------------------
//        |  |   |           |             |
//        |  |   |   A       |    C        |
//        |  |   |           |             |
//        |  |---|-----------|-------------|
//        |  |   |           |             |
//        |  /---15---------------------------> PT
//        | /
//         /------------------------------->
//        /
//      TOF

	 if(TypeMode!=3) {
            if(track->pt()>100){
               st->CtrlPt_S4_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S4_Im->Fill(Ih, Event_Weight);
               if(tof){
		 st->CtrlPt_S4_TOF->Fill(MuonTOF, Event_Weight);
		 if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);
		 else st->CtrlFor_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);}
            }else if(track->pt()>80){
               st->CtrlPt_S3_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S3_Im->Fill(Ih, Event_Weight);
               if(tof){
		 st->CtrlPt_S3_TOF->Fill(MuonTOF, Event_Weight);
		 if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S3_TOF->Fill(MuonTOF, Event_Weight);
		 else st->CtrlFor_Pt_S3_TOF->Fill(MuonTOF, Event_Weight);}
            }else if(track->pt()>60){
               st->CtrlPt_S2_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S2_Im->Fill(Ih, Event_Weight);
               if(tof){
		 st->CtrlPt_S2_TOF->Fill(MuonTOF, Event_Weight);
		 if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S2_TOF->Fill(MuonTOF, Event_Weight);
		 else st->CtrlFor_Pt_S2_TOF->Fill(MuonTOF, Event_Weight);}
            }else{
               st->CtrlPt_S1_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S1_Im->Fill(Ih, Event_Weight);
               if(tof){
		 st->CtrlPt_S1_TOF->Fill(MuonTOF, Event_Weight);
		 if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);
		 else st->CtrlFor_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);}
            }

            if(Is>0.2){           if(tof)st->CtrlIs_S4_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Is>0.1){     if(tof)st->CtrlIs_S3_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Is>0.05){     if(tof)st->CtrlIs_S2_TOF->Fill(MuonTOF, Event_Weight);
            }else{                             if(tof)st->CtrlIs_S1_TOF->Fill(MuonTOF, Event_Weight);
            }

            if(Ih>4.4){           if(tof)st->CtrlIm_S4_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Ih>4.1){     if(tof)st->CtrlIm_S3_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Ih>3.8){     if(tof)st->CtrlIm_S2_TOF->Fill(MuonTOF, Event_Weight);
            }else{                             if(tof)st->CtrlIm_S1_TOF->Fill(MuonTOF, Event_Weight);
            }}

	 //Fill Control plots for muon only analysis.  Uses different pt regions
	 else {
	 //Fill control plots of TOF in different pt regions
	 if(track->pt()>250){
	   st->CtrlPt_S4_TOF->Fill(MuonTOF, Event_Weight);
	   if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);
	   else st->CtrlFor_Pt_S4_TOF->Fill(MuonTOF, Event_Weight);
	 }else if(track->pt()>130){
	   st->CtrlPt_S3_TOF->Fill(MuonTOF, Event_Weight);
	   if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S3_TOF->Fill(MuonTOF, Event_Weight);
	   else st->CtrlFor_Pt_S3_TOF->Fill(MuonTOF, Event_Weight);
	 }else if(track->pt()>90){
	   st->CtrlPt_S2_TOF->Fill(MuonTOF, Event_Weight);
	   if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S2_TOF->Fill(MuonTOF, Event_Weight);
	   else st->CtrlFor_Pt_S2_TOF->Fill(MuonTOF, Event_Weight);
	 }else{
	   st->CtrlPt_S1_TOF->Fill(MuonTOF, Event_Weight);
	   if(fabs(track->eta())<DTRegion) st->CtrlCen_Pt_S1_TOF->Fill(MuonTOF, Event_Weight);
	   else st->CtrlFor_Pt_S1_TOF->Fill(MuonTOF, Event_Weight);
	 }}

	 //Now do the same for different TOF regions
         if(MuonTOF<0.9){
           st->CtrlTOF_S4_Pt->Fill(track->pt(), Event_Weight);
           if(fabs(track->eta())<DTRegion) st->CtrlCen_TOF_S4_Pt->Fill(track->pt(), Event_Weight);
           else st->CtrlFor_TOF_S4_Pt->Fill(track->pt(), Event_Weight);
         }else if(MuonTOF<1.0){
           st->CtrlTOF_S3_Pt->Fill(track->pt(), Event_Weight);
           if(fabs(track->eta())<DTRegion) st->CtrlCen_TOF_S3_Pt->Fill(track->pt(), Event_Weight);
           else st->CtrlFor_TOF_S3_Pt->Fill(track->pt(), Event_Weight);
         }else if(MuonTOF<1.1){
           st->CtrlTOF_S2_Pt->Fill(track->pt(), Event_Weight);
           if(fabs(track->eta())<DTRegion) st->CtrlCen_TOF_S2_Pt->Fill(track->pt(), Event_Weight);
           else st->CtrlFor_TOF_S2_Pt->Fill(track->pt(), Event_Weight);
         }else{
           st->CtrlTOF_S1_Pt->Fill(track->pt(), Event_Weight);
           if(fabs(track->eta())<DTRegion) st->CtrlCen_TOF_S1_Pt->Fill(track->pt(), Event_Weight);
           else st->CtrlFor_TOF_S1_Pt->Fill(track->pt(), Event_Weight);
         }


         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
	   if(MuonTOF<GlobalMinTOF) continue;
            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassICut   = (Is>=CutI[CutIndex]);
            bool PassTOFCut = MuonTOF>=CutTOF[CutIndex];

            if(       PassTOFCut &&  PassPtCut &&  PassICut){   //Region D
               st->H_D      ->Fill(CutIndex,                Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_D_Cen->Fill(CutIndex, Event_Weight);
               else st->H_D_For->Fill(CutIndex, Event_Weight);
               st->RegionD_P  ->Fill(CutIndex,track->p(),     Event_Weight);
               st->RegionD_I  ->Fill(CutIndex,Ih,Event_Weight);
               st->RegionD_TOF->Fill(CutIndex,MuonTOF,        Event_Weight);
	       st->AS_Eta_RegionD->Fill(CutIndex,track->eta());
            }else if( PassTOFCut &&  PassPtCut && !PassICut){   //Region C
               st->H_C     ->Fill(CutIndex,                 Event_Weight);
	       if(fabs(track->eta())<DTRegion) st->H_C_Cen->Fill(CutIndex, Event_Weight);
	       else st->H_C_For->Fill(CutIndex, Event_Weight);
               if(TypeMode!=2)st->Pred_EtaP  ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight);
//               Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               st->AS_Eta_RegionC->Fill(CutIndex,track->eta());
            }else if( PassTOFCut && !PassPtCut &&  PassICut){   //Region B
               st->H_B     ->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_B_Cen->Fill(CutIndex, Event_Weight);
               else st->H_B_For->Fill(CutIndex, Event_Weight);
               if(TypeMode!=2)st->Pred_I  ->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode!=2)st->Pred_EtaS->Fill(CutIndex,track->eta(),         Event_Weight);
//               Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               st->AS_Eta_RegionB->Fill(CutIndex,track->eta());
            }else if( PassTOFCut && !PassPtCut && !PassICut){   //Region A
               st->H_A     ->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_A_Cen->Fill(CutIndex, Event_Weight);
               else st->H_A_For->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               if(TypeMode!=2)st->Pred_EtaB->Fill(CutIndex,track->eta(),         Event_Weight);
               if(TypeMode==2)st->Pred_EtaS2->Fill(CutIndex,track->eta(),        Event_Weight);
               st->AS_Eta_RegionA->Fill(CutIndex,track->eta());
            }else if(!PassTOFCut &&  PassPtCut &&  PassICut){   //Region H
               st->H_H   ->Fill(CutIndex,          Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_H_Cen->Fill(CutIndex, Event_Weight);
               else st->H_H_For->Fill(CutIndex, Event_Weight);
//               Pred_P->Fill(CutIndex,track->p(),        Event_Weight);
//               Pred_I->Fill(CutIndex,Ih,   Event_Weight);
               st->AS_Eta_RegionH->Fill(CutIndex,track->eta());
            }else if(!PassTOFCut &&  PassPtCut && !PassICut){   //Region G
               st->H_G     ->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_G_Cen->Fill(CutIndex, Event_Weight);
               else st->H_G_For->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_EtaP  ->Fill(CutIndex,track->eta(),track->p(),     Event_Weight);
               st->AS_Eta_RegionG->Fill(CutIndex,track->eta());
            }else if(!PassTOFCut && !PassPtCut &&  PassICut){   //Region F
               st->H_F     ->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_F_Cen->Fill(CutIndex, Event_Weight);
               else st->H_F_For->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_I  ->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode==2)st->Pred_EtaS->Fill(CutIndex,track->eta(),         Event_Weight);
               st->AS_Eta_RegionF->Fill(CutIndex,track->eta());
            }else if(!PassTOFCut && !PassPtCut && !PassICut){   //Region E
               st->H_E     ->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_E_Cen->Fill(CutIndex, Event_Weight);
               else st->H_E_For->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_EtaB->Fill(CutIndex,track->eta(),         Event_Weight);
               st->AS_Eta_RegionE->Fill(CutIndex,track->eta());
            }

	    /*
	    //Plots for cosmic background prediction and systematic for muon only search
	    //Fill number of events in sideband regions
	    if(DzType>-1) {
            if(       PassTOFCut &&  PassPtCut &&  PassICut){   //Region D
               st->H_D_Syst[DzType]->Fill(CutIndex,                Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_D_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_D_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if( PassTOFCut &&  PassPtCut && !PassICut){   //Region C
               st->H_C_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
	       if(fabs(track->eta())<DTRegion) st->H_C_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
	       else st->H_C_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if( PassTOFCut && !PassPtCut &&  PassICut){   //Region B
               st->H_B_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_B_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_B_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if( PassTOFCut && !PassPtCut && !PassICut){   //Region A
               st->H_A_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_A_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_A_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if(!PassTOFCut &&  PassPtCut &&  PassICut){   //Region H
               st->H_H_Syst[DzType]->Fill(CutIndex,          Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_H_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_H_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if(!PassTOFCut &&  PassPtCut && !PassICut){   //Region G
               st->H_G_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_G_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_G_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if(!PassTOFCut && !PassPtCut &&  PassICut){   //Region F
               st->H_F_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_F_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_F_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
            }else if(!PassTOFCut && !PassPtCut && !PassICut){   //Region E
               st->H_E_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_E_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
               else st->H_E_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
	    }
	    }
	    */
         }

	 //Use events with low TOF to check accuracy of background prediction
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
	   if(MuonTOF>GlobalMinTOF) continue;
            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassICut   = (Is>=CutI[CutIndex]);
            bool PassTOFCut = MuonTOF<=(2-CutTOF[CutIndex]);

            if(       PassTOFCut &&  PassPtCut &&  PassICut){   //Region D
               st->H_D_Flip->Fill(CutIndex,                Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_D_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_D_For_Flip->Fill(CutIndex, Event_Weight);
            }else if( PassTOFCut &&  PassPtCut && !PassICut){   //Region C
               st->H_C_Flip->Fill(CutIndex,                 Event_Weight);
	       if(fabs(track->eta())<DTRegion) st->H_C_Cen_Flip->Fill(CutIndex, Event_Weight);
	       else st->H_C_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode!=2)st->Pred_EtaP_Flip->Fill(CutIndex,track->eta(), track->p(),     Event_Weight);
//               Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
            }else if( PassTOFCut && !PassPtCut &&  PassICut){   //Region B
               st->H_B_Flip->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_B_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_B_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode!=2)st->Pred_I_Flip->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode!=2)st->Pred_EtaS_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
//               Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
            }else if( PassTOFCut && !PassPtCut && !PassICut){   //Region A
               st->H_A_Flip->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_A_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_A_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
               if(TypeMode!=2)st->Pred_EtaB_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
               if(TypeMode==2)st->Pred_EtaS2_Flip->Fill(CutIndex,track->eta(),        Event_Weight);
            }else if(!PassTOFCut &&  PassPtCut &&  PassICut){   //Region H
               st->H_H_Flip->Fill(CutIndex,          Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_H_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_H_For_Flip->Fill(CutIndex, Event_Weight);
//               Pred_P_Flip->Fill(CutIndex,track->p(),        Event_Weight);
//               Pred_I_Flip->Fill(CutIndex,Ih,   Event_Weight);
            }else if(!PassTOFCut &&  PassPtCut && !PassICut){   //Region G
               st->H_G_Flip->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_G_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_G_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_EtaP_Flip->Fill(CutIndex,track->eta(),track->p(),     Event_Weight);
            }else if(!PassTOFCut && !PassPtCut &&  PassICut){   //Region F
               st->H_F_Flip->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_F_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_F_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_I_Flip->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode==2)st->Pred_EtaS_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
            }else if(!PassTOFCut && !PassPtCut && !PassICut){   //Region E
               st->H_E_Flip->Fill(CutIndex,                 Event_Weight);
               if(fabs(track->eta())<DTRegion) st->H_E_Cen_Flip->Fill(CutIndex, Event_Weight);
               else st->H_E_For_Flip->Fill(CutIndex, Event_Weight);
               if(TypeMode==2)st->Pred_EtaB_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
            }
         }
	 

         //if(fabs(track->eta())<DTRegion) st->H_DzCounts_DT->Fill(DzType, Event_Weight);
         //else st->H_DzCounts_CSC->Fill(DzType, Event_Weight);
}





void Analysis_Step3(char* SavePath)
{
   printf("Step3: Building Mass Spectrum for B and S\n");

   int TreeStep;
   //////////////////////////////////////////////////     BUILD BACKGROUND MASS SPECTRUM

   if(DataFileName.size()) {
     stPlots_Init(HistoFile, DataPlots,"Data", CutPt.size());
     if(TypeMode==3) stPlots_Init(HistoFile, DataPlotsControl,"Data_Control", CutPt.size());
   }
   HistoFile->cd();

   DuplicatesClass Duplicates;
   Duplicates.Clear();

   fwlite::ChainEvent treeD(DataFileName);
   double SampleWeight = GetSampleWeight(-1);
   Event_Weight = SampleWeight;
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Building Mass Spectrum for D :");
   TreeStep = treeD.size()/50;if(TreeStep==0)TreeStep=1;

   bool* HSCPTk = new bool[CutPt.size()]; 
   double* MaxMass = new double[CutPt.size()];

   for(Long64_t ientry=0;ientry<treeD.size();ientry++){
      treeD.to(ientry);
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
      if(treeD.eventAuxiliary().run()>193092 && treeD.eventAuxiliary().run() < 194619) continue;
      if(Duplicates.isDuplicate(treeD.eventAuxiliary().run(),treeD.eventAuxiliary().event())){continue;}

      DataPlots.TotalE->Fill(0.0,Event_Weight);  
      if(!PassTrigger(treeD) ) continue;
      DataPlots.TotalTE->Fill(0.0,Event_Weight);

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeD,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

      fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
      dEdxSCollH.getByLabel(treeD, dEdxS_Label.c_str());
      if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

      fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
      dEdxMCollH.getByLabel(treeD, dEdxM_Label.c_str());
      if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(treeD, "muontiming",TOF_Label.c_str());
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

      fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
      TOFDTCollH.getByLabel(treeD, "muontiming",TOFdt_Label.c_str());
      if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");return;}

      fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
      TOFCSCCollH.getByLabel(treeD, "muontiming",TOFcsc_Label.c_str());
      if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");return;}


      for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk[CutIndex] = false;   }
      for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass[CutIndex] = -1; }
      for(unsigned int c=0;c<hscpColl.size();c++){

         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();

         if(track.isNull() && TypeMode!=3)continue;
	 if(TypeMode!=0 && muon.isNull()) continue;

	 const DeDxData* dedxSObj = NULL;
         const DeDxData* dedxMObj = NULL;
	 if(TypeMode!=3 && !track.isNull()) {
	   dedxSObj  = &dEdxSCollH->get(track.key());
	   dedxMObj  = &dEdxMCollH->get(track.key());
	 }

         const reco::MuonTimeExtra* tof = NULL;
         const reco::MuonTimeExtra* dttof = NULL;
         const reco::MuonTimeExtra* csctof = NULL;
        if(TypeMode!=0 && !hscp.muonRef().isNull()){ tof  = &TOFCollH->get(hscp.muonRef().key()); dttof = &TOFDTCollH->get(hscp.muonRef().key());  csctof = &TOFCSCCollH->get(hscp.muonRef().key());}

	//double MuonTOF = GlobalMinTOF;
	//if(tof){MuonTOF = tof->inverseBeta(); }
	bool Preselected=false;

	//if(TypeMode!=3 && !PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, treeD, &DataPlots))continue;

	 //Dz type is used to determine in which control region in the dz distribution a track falls
	 //The different regions are used to predict the cosmic background and its uncertainty
	 int DzType=-1;
         if(TypeMode==3 && !muon->isGlobalMuon()) PassSAPreselection(hscp, tof, dttof, csctof, treeD, &DataPlotsControl, DzType, true);
	 if(TypeMode==3 && !PassSAPreselection(hscp, tof, dttof, csctof, treeD, &DataPlots, DzType)) continue;

         Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, &DataPlots);

	 double Mass=-1, MassTOF=-1, MassComb=-1;

         if(dedxMObj) Mass     = GetMass(track->p(),dedxMObj->dEdx());
	 double p=-1;
	 if(!track.isNull()) p=track->p();

         if(tof) MassTOF=GetTOFMass(p,tof->inverseBeta());
	 if(tof && dedxMObj) MassComb=GetMassFromBeta(p, (GetIBeta(dedxMObj->dEdx()) + (1/tof->inverseBeta()))*0.5 );
	 else if(dedxMObj) MassComb = Mass;
	 else if(tof) MassComb = MassTOF;
         bool PassNonTrivialSelection=false;

         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
            //Full Selection
            if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, treeD, CutIndex, &DataPlots))continue;
            if(CutIndex!=0)PassNonTrivialSelection=true;
            HSCPTk[CutIndex] = true;
	    if(Mass>MaxMass[CutIndex]) MaxMass[CutIndex]=Mass;

      	    DataPlots.Mass->Fill(CutIndex, Mass,Event_Weight);
            if(tof){
               DataPlots.MassTOF ->Fill(CutIndex, MassTOF , Event_Weight);
            }
            DataPlots.MassComb->Fill(CutIndex, MassComb, Event_Weight);
         } //end of Cut loop

//         if(track->pt()>40 && Mass>75)stPlots_FillTree(DataPlots, treeD.eventAuxiliary().run(),treeD.eventAuxiliary().event(), c, track->pt(), dedxSObj->dEdx(), tof ? tof->inverseBeta() : -1);
         if (PassNonTrivialSelection) stPlots_FillTree(DataPlots, treeD.eventAuxiliary().run(),treeD.eventAuxiliary().event(), c, track->pt(), dedxSObj ? dedxSObj->dEdx() : -1, tof ? tof->inverseBeta() : -1, Mass, -1);

      } // end of Track Loop
      for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  if(HSCPTk[CutIndex]){DataPlots.HSCPE->Fill(CutIndex,Event_Weight); DataPlots.MaxEventMass->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);} }
   }// end of Event Loop
   delete [] HSCPTk;
   delete [] MaxMass;
   printf("\n");

   if(DataFileName.size())stPlots_Clear(DataPlots, true);


   //////////////////////////////////////////////////     BUILD MCTRUTH MASS SPECTRUM
   if(MCsample.size())stPlots_Init(HistoFile, MCTrPlots,"MCTr", CutPt.size());

   for(unsigned int m=0;m<MCsample.size();m++){
      stPlots_Init(HistoFile,MCPlots[m],MCsample[m].Name, CutPt.size());

      std::vector<string> FileName;

      GetInputFiles(FileName, MCsample[m].Name);

      fwlite::ChainEvent treeM(FileName);

      //get PU reweighted total # MC events.
      double NMCevents=0;
      for(Long64_t ientry=0;ientry<treeM.size();ientry++){
	treeM.to(ientry);
	if(MaxEntry>0 && ientry>MaxEntry)break;
	double PUSystFactor;
	double puwt= GetPUWeight(treeM, MCsample[m].IsS4PileUp, PUSystFactor);
	NMCevents+=puwt;
      }

      double SampleWeight = GetSampleWeightMC(IntegratedLuminosity,FileName, MCsample[m].XSection, treeM.size(), NMCevents);

      printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
      printf("Building Mass for %10s :",MCsample[m].Name.c_str());
      TreeStep = treeM.size()/50;if(TreeStep==0)TreeStep=1;

      bool* HSCPTk = new bool[CutPt.size()]; 
      double* MaxMass = new double[CutPt.size()];
      for(Long64_t ientry=0;ientry<treeM.size();ientry++){       
          treeM.to(ientry);
         if(MaxEntry>0 && ientry>MaxEntry)break;
         if(MCsample[m].MaxEvent>0 && ientry>MCsample[m].MaxEvent)break;
         if(ientry%TreeStep==0){printf(".");fflush(stdout);}

         if(!hasGoodPtHat(treeM, MCsample[m].MaxPtHat)){continue;}
	 double PUSystFactor;
         Event_Weight = SampleWeight * GetPUWeight(treeM, MCsample[m].IsS4PileUp, PUSystFactor);

         MCTrPlots .TotalE->Fill(0.0,Event_Weight);
         MCPlots[m].TotalE->Fill(0.0,Event_Weight);
         if(!PassTrigger(treeM) )continue;
         MCTrPlots .TotalTE->Fill(0.0,Event_Weight);
         MCPlots[m].TotalTE->Fill(0.0,Event_Weight);

         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeM,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

         fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
         dEdxSCollH.getByLabel(treeM, dEdxS_Label.c_str());
         if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

         fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
         dEdxMCollH.getByLabel(treeM, dEdxM_Label.c_str());
         if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFCollH;
         TOFCollH.getByLabel(treeM, "muontiming",TOF_Label.c_str());
         if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");continue;}
         
         fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
         TOFDTCollH.getByLabel(treeM, "muontiming",TOFdt_Label.c_str());
         if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
         TOFCSCCollH.getByLabel(treeM, "muontiming",TOFcsc_Label.c_str());
         if(!TOFCSCCollH.isValid()){printf("Invalid CSCTOF collection\n");continue;}

         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk[CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass[CutIndex] = -1; }
         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

	    const DeDxData* dedxSObj = NULL;
	    const DeDxData* dedxMObj = NULL;
	    if(TypeMode!=3 && !track.isNull()) {
	      dedxSObj  = &dEdxSCollH->get(track.key());
	      dedxMObj  = &dEdxMCollH->get(track.key());
	    }

            const reco::MuonTimeExtra* tof = NULL;
            const reco::MuonTimeExtra* dttof = NULL;
            const reco::MuonTimeExtra* csctof = NULL;
            if(TypeMode==2 && !hscp.muonRef().isNull()){ tof  = &TOFCollH->get(hscp.muonRef().key()); dttof  = &TOFDTCollH->get(hscp.muonRef().key()); csctof  = &TOFCSCCollH->get(hscp.muonRef().key());}

                PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, treeM,           &MCPlots[m]);
            if(!PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, treeM,           &MCTrPlots))continue;
            Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, &MCTrPlots);

            double Mass     = GetMass(track->p(),dedxMObj->dEdx());
            double MassTOF  = -1;   if(tof)MassTOF  = GetTOFMass(track->p(),tof->inverseBeta());
            double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()) + (1/tof->inverseBeta()))*0.5 ) ;


            for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){

                   PassSelection   (hscp, dedxSObj, dedxMObj, tof, treeM, CutIndex, &MCPlots[m]);
               if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, treeM, CutIndex, &MCTrPlots))continue;
               HSCPTk[CutIndex] = true;
	       if(Mass>MaxMass[CutIndex]) MaxMass[CutIndex]=Mass;

               MCTrPlots .Mass->Fill(CutIndex , Mass,Event_Weight);
               MCPlots[m].Mass->Fill(CutIndex, Mass,Event_Weight);

               if(tof){
                  MCTrPlots .MassTOF ->Fill(CutIndex, MassTOF , Event_Weight);
                  MCPlots[m].MassTOF ->Fill(CutIndex, MassTOF , Event_Weight);
               }
               MCTrPlots .MassComb->Fill(CutIndex, MassComb, Event_Weight);
               MCPlots[m].MassComb->Fill(CutIndex, MassComb, Event_Weight);
         } //end of Cut loo
	    if(track->pt()>35)stPlots_FillTree(MCTrPlots , treeM.eventAuxiliary().run(),treeM.eventAuxiliary().event(), c, track->pt(), dedxSObj->dEdx(), tof ? tof->inverseBeta() : -1, Mass);
	    if(track->pt()>35)stPlots_FillTree(MCPlots[m], treeM.eventAuxiliary().run(),treeM.eventAuxiliary().event(), c, track->pt(), dedxSObj->dEdx(), tof ? tof->inverseBeta() : -1, Mass);

         } // end of Track Loop 
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  if(HSCPTk[CutIndex]){
	     MCTrPlots .HSCPE->Fill(CutIndex,Event_Weight);MCPlots[m].HSCPE->Fill(CutIndex,Event_Weight);
	     MCTrPlots.MaxEventMass->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);MCPlots[m].MaxEventMass->Fill(CutIndex,MaxMass[CutIndex], Event_Weight); 
	   } }
      }// end of Event Loop
      delete [] HSCPTk;
      delete [] MaxMass;
      stPlots_Clear(MCPlots[m], true);
      printf("\n");
   }
   if(MCsample.size())stPlots_Clear(MCTrPlots, true);

   //////////////////////////////////////////////////     BUILD SIGNAL MASS SPECTRUM
   TRandom3* RNG = new TRandom3();
   for(unsigned int s=0;s<signals.size();s++){
      stPlots_Init(HistoFile,SignPlots[s],signals[s].Name       , CutPt.size());

      bool* HSCPTk          = new bool[CutPt.size()];
      bool* HSCPTk_SystP    = new bool[CutPt.size()];
      bool* HSCPTk_SystI    = new bool[CutPt.size()];
      bool* HSCPTk_SystT    = new bool[CutPt.size()];
      bool* HSCPTk_SystM    = new bool[CutPt.size()];
      bool* HSCPTk_SystPU   = new bool[CutPt.size()];
      double* MaxMass       = new double[CutPt.size()];
      double* MaxMass_SystP = new double[CutPt.size()];
      double* MaxMass_SystI = new double[CutPt.size()];
      double* MaxMass_SystT = new double[CutPt.size()];
      double* MaxMass_SystM = new double[CutPt.size()];
      double* MaxMass_SystPU= new double[CutPt.size()];


      printf("Progressing Bar                                    :0%%       20%%       40%%       60%%       80%%       100%%\n");
      //Do two loops through signal for samples with and without trigger change.  Period before has 325 1/pb and rest of luminosity is after
      for (int period=0; period<RunningPeriods; period++) {

      std::vector<string> SignFileName;
      GetInputFiles(SignFileName, signals[s].FileName, period);
      fwlite::ChainEvent treeS(SignFileName);

      if (period==0) printf("Building Mass for %10s for before RPC change :",signals[s].Name.c_str());
      if (period==1) printf("\nBuilding Mass for %10s for after RPC change  :",signals[s].Name.c_str());

      //get PU reweighted total # MC events.
      double NMCevents=0;
      for(Long64_t ientry=0;ientry<treeS.size();ientry++){
	treeS.to(ientry);
	if(MaxEntry>0 && ientry>MaxEntry)break;
        double PUSystFactor;
	double puwt= GetPUWeight(treeS, signals[s].IsS4PileUp, PUSystFactor);
	NMCevents+=puwt;
      }


      TreeStep = treeS.size()/50;if(TreeStep==0)TreeStep=1;

      double SampleWeight = GetSampleWeight(IntegratedLuminosity,IntegratedLuminosityBeforeTriggerChange,signals[s].XSec,NMCevents, period);
      for(Long64_t ientry=0;ientry<treeS.size();ientry++){
         treeS.to(ientry);
         if(MaxEntry>0 && ientry>MaxEntry)break;
         if(ientry%TreeStep==0){printf(".");fflush(stdout);}
	 double PUSystFactor;


         Event_Weight = SampleWeight * GetPUWeight(treeS, signals[s].IsS4PileUp, PUSystFactor);



         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         genCollHandle.getByLabel(treeS, "genParticles");
         if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
         std::vector<reco::GenParticle> genColl = *genCollHandle;
         int NChargedHSCP=HowManyChargedHSCP(genColl);
         //skip event with the wrong number of charged HSCP
         if(signals[s].NChargedHSCP>=0 && signals[s].NChargedHSCP!=NChargedHSCP)continue;

         double HSCPGenBeta1, HSCPGenBeta2;
         GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,false);
         if(HSCPGenBeta1>=0)SignPlots[s].Beta_Gen->Fill(HSCPGenBeta1, Event_Weight);        if(HSCPGenBeta2>=0)SignPlots[s].Beta_Gen->Fill(HSCPGenBeta2, Event_Weight);
         GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,true);
         if(HSCPGenBeta1>=0)SignPlots[s].Beta_GenCharged->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SignPlots[s].Beta_GenCharged->Fill(HSCPGenBeta2, Event_Weight);



         SignPlots[s]               .TotalE   ->Fill(0.0,Event_Weight);
         SignPlots[s]               .TotalEPU ->Fill(0.0,Event_Weight*PUSystFactor);
         if(!PassTrigger(treeS) )continue;
         SignPlots[s]               .TotalTE->Fill(0.0,Event_Weight);

         if(HSCPGenBeta1>=0)SignPlots[s].Beta_Triggered->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SignPlots[s].Beta_Triggered->Fill(HSCPGenBeta2, Event_Weight);



         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeS,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;



         fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
         dEdxSCollH.getByLabel(treeS, dEdxS_Label.c_str());
         if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

         fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
         dEdxMCollH.getByLabel(treeS, dEdxM_Label.c_str());
         if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFCollH;
         TOFCollH.getByLabel(treeS, "muontiming",TOF_Label.c_str());
         if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
         TOFDTCollH.getByLabel(treeS, "muontiming",TOFdt_Label.c_str());
         if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
         TOFCSCCollH.getByLabel(treeS, "muontiming",TOFcsc_Label.c_str());
         if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");continue;}




         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk        [CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystP  [CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystI  [CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystT  [CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystM  [CutIndex] = false;   }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystPU [CutIndex] = false; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass       [CutIndex] = -1; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystP [CutIndex] = -1; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystI [CutIndex] = -1; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystT [CutIndex] = -1; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystM [CutIndex] = -1; }
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystPU[CutIndex] = -1; }
         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

            int ClosestGen;
            if(DistToHSCP(hscp, genColl, ClosestGen)>0.03)continue;

	    const DeDxData* dedxSObj = NULL;
	    const DeDxData* dedxMObj = NULL;
	    if(TypeMode!=3 && !track.isNull()) {
	      dedxSObj  = &dEdxSCollH->get(track.key());
	      dedxMObj  = &dEdxMCollH->get(track.key());
	    }

            const reco::MuonTimeExtra* tof = NULL;
            const reco::MuonTimeExtra* dttof = NULL;
            const reco::MuonTimeExtra* csctof = NULL;
            if(TypeMode==2 && !hscp.muonRef().isNull()){ tof  = &TOFCollH->get(hscp.muonRef().key()); dttof  = &TOFDTCollH->get(hscp.muonRef().key()); csctof  = &TOFCSCCollH->get(hscp.muonRef().key()); }

            ///////////// START COMPUTATION OF THE SYSTEMATIC //////////

            bool PRescale = true;
            double IRescale = RNG->Gaus(0, 0.083)+0.015; // added to the Ias value
            double MRescale = 1.036;
            double TRescale = -0.02; // added to the 1/beta value
	    if(tof) if(csctof->nDof()==0) TRescale = -0.003;

            // Systematic on P
            if(PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,  NULL, -1,   PRescale, 0, 0)){
               double Mass     = GetMass(track->p()*PRescale,dedxMObj->dEdx());
               double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p()*PRescale,tof->inverseBeta());
               double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p()*PRescale, (GetIBeta(dedxMObj->dEdx()) + (1/tof->inverseBeta()))*0.5 ) ;

               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, NULL, -1,   PRescale, 0, 0)){
                     HSCPTk_SystP[CutIndex] = true;
		     if(Mass>MaxMass_SystP[CutIndex]) MaxMass_SystP[CutIndex]=Mass;
                     SignPlots[s].Mass_SystP->Fill(CutIndex, Mass,Event_Weight);
                     if(tof){
                        SignPlots[s].MassTOF_SystP ->Fill(CutIndex, MassTOF , Event_Weight);
                     }
                     SignPlots[s].MassComb_SystP->Fill(CutIndex, MassComb, Event_Weight);
                  }
               }
            }

            // Systematic on I (both Ias and Ih)
            if(PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,  NULL, -1,   0, IRescale, 0)){
               double Mass     = GetMass(track->p(),dedxMObj->dEdx()*MRescale);
               double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
               double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()*MRescale) + (1/tof->inverseBeta()))*0.5 ) ;

               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, NULL, -1,   0, IRescale, 0)){
                     HSCPTk_SystI[CutIndex] = true;
                     if(Mass>MaxMass_SystI[CutIndex]) MaxMass_SystI[CutIndex]=Mass;
                     SignPlots[s].Mass_SystI->Fill(CutIndex, Mass,Event_Weight);
                     if(tof){
                        SignPlots[s].MassTOF_SystI ->Fill(CutIndex, MassTOF , Event_Weight);
                     }
                     SignPlots[s].MassComb_SystI->Fill(CutIndex, MassComb, Event_Weight);
                  }
               }
            }


            // Systematic on M
            if(PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,  NULL, -1,   0, 0, 0)){
               double Mass     = GetMass(track->p(),dedxMObj->dEdx()*MRescale);
               double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
               double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()*MRescale) + (1/tof->inverseBeta()))*0.5 ) ;

               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, NULL, -1,   0, 0, 0)){
                     HSCPTk_SystM[CutIndex] = true;
                     if(Mass>MaxMass_SystM[CutIndex]) MaxMass_SystM[CutIndex]=Mass;
                     SignPlots[s].Mass_SystM->Fill(CutIndex, Mass,Event_Weight);
                     if(tof){
                        SignPlots[s].MassTOF_SystM ->Fill(CutIndex, MassTOF , Event_Weight);
                     }
                     SignPlots[s].MassComb_SystM->Fill(CutIndex, MassComb, Event_Weight);
                  }
               }
            }


            // Systematic on T
            if(PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,  NULL, -1,   0, 0, TRescale)){
               double Mass     = GetMass(track->p(),dedxMObj->dEdx());
               double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta()*TRescale);
               double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()) + ((1/tof->inverseBeta())*TRescale ))*0.5 ) ;

               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, NULL, -1,   0, 0, TRescale)){
                     HSCPTk_SystT[CutIndex] = true;
                     if(Mass>MaxMass_SystT[CutIndex]) MaxMass_SystT[CutIndex]=Mass;
                     SignPlots[s].Mass_SystT->Fill(CutIndex, Mass,Event_Weight);
                     if(tof){
                        SignPlots[s].MassTOF_SystT ->Fill(CutIndex, MassTOF , Event_Weight);
                     }
                     SignPlots[s].MassComb_SystT->Fill(CutIndex, MassComb, Event_Weight);
                  }
               }
            }

            // Systematic on PU
            if(PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,  NULL, -1,   PRescale, 0, 0)){
               double Mass     = GetMass(track->p()*PRescale,dedxMObj->dEdx());
               double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p()*PRescale,tof->inverseBeta());
               double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p()*PRescale, (GetIBeta(dedxMObj->dEdx()) + (1/tof->inverseBeta()))*0.5 ) ;

               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, NULL, -1,   PRescale, 0, 0)){
                     HSCPTk_SystPU[CutIndex] = true;
		     if(Mass>MaxMass_SystPU[CutIndex]) MaxMass_SystPU[CutIndex]=Mass;
                     SignPlots[s].Mass_SystPU->Fill(CutIndex, Mass,Event_Weight*PUSystFactor);
                     if(tof){
                        SignPlots[s].MassTOF_SystPU ->Fill(CutIndex, MassTOF , Event_Weight*PUSystFactor);
                     }
                     SignPlots[s].MassComb_SystPU->Fill(CutIndex, MassComb, Event_Weight*PUSystFactor);
                  }
               }
            }

            ///////////// END   COMPUTATION OF THE SYSTEMATIC //////////



            if(!PassPreselection(hscp,  dedxSObj, dedxMObj, tof, dttof, csctof, treeS,           &SignPlots[s], genColl[ClosestGen].p()/genColl[ClosestGen].energy()))continue;         

            double Mass     = GetMass(track->p(),dedxMObj->dEdx());
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
            double MassComb = Mass;if(tof)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()) + (1/tof->inverseBeta()))*0.5 ) ;


            for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
               if(!PassSelection   (hscp,  dedxSObj, dedxMObj, tof, treeS, CutIndex, &SignPlots[s], genColl[ClosestGen].p()/genColl[ClosestGen].energy()))continue;    

               HSCPTk[CutIndex] = true;
	       if(Mass>MaxMass[CutIndex]) MaxMass[CutIndex]=Mass;

               SignPlots[s].Mass->Fill(CutIndex, Mass,Event_Weight);
               if(tof){
                  SignPlots[s].MassTOF ->Fill(CutIndex, MassTOF , Event_Weight);
               }
               SignPlots[s].MassComb->Fill(CutIndex, MassComb, Event_Weight);
            } //end of Cut loop
            if(track->pt()>35 && Mass>35)stPlots_FillTree(SignPlots[s] , treeS.eventAuxiliary().run(),treeS.eventAuxiliary().event(), c, track->pt(), dedxSObj->dEdx(), tof ? tof->inverseBeta() : -1, Mass);
         } // end of Track Loop 
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
	   if(HSCPTk[CutIndex]){
	     SignPlots[s].HSCPE             ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass      ->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);
           }
           if(HSCPTk_SystP[CutIndex]){
             SignPlots[s].HSCPE_SystP       ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass_SystP->Fill(CutIndex,MaxMass_SystP[CutIndex], Event_Weight);
           }
           if(HSCPTk_SystI[CutIndex]){
             SignPlots[s].HSCPE_SystI       ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass_SystI->Fill(CutIndex,MaxMass_SystI[CutIndex], Event_Weight);
           }
           if(HSCPTk_SystM[CutIndex]){
             SignPlots[s].HSCPE_SystM       ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass_SystM->Fill(CutIndex,MaxMass_SystM[CutIndex], Event_Weight);
           }
           if(HSCPTk_SystT[CutIndex]){
             SignPlots[s].HSCPE_SystT       ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass_SystT->Fill(CutIndex,MaxMass_SystT[CutIndex], Event_Weight);
           }
           if(HSCPTk_SystPU[CutIndex]){
             SignPlots[s].HSCPE_SystPU       ->Fill(CutIndex,Event_Weight);
             SignPlots[s].MaxEventMass_SystPU->Fill(CutIndex,MaxMass_SystPU[CutIndex], Event_Weight);
           }
        }
      }// end of Event Loop
      }
      printf("\n");
      delete [] HSCPTk;
      delete [] HSCPTk_SystP;
      delete [] HSCPTk_SystI;
      delete [] HSCPTk_SystT;
      delete [] HSCPTk_SystM;
      delete [] HSCPTk_SystPU;
      delete [] MaxMass;
      delete [] MaxMass_SystP;
      delete [] MaxMass_SystI;
      delete [] MaxMass_SystT;
      delete [] MaxMass_SystM;
      delete [] MaxMass_SystPU;

      stPlots_Clear(SignPlots[s], true);
   }// end of signal Type loop
   delete RNG;
}

TH1D* GetPDF(TH1D* pdf){
   char NewName[2048];
   sprintf(NewName,"%s_PDF", pdf->GetName());

   TH1D* PDF = new TH1D(NewName,NewName,pdf->GetNbinsX(),pdf->GetXaxis()->GetXmin(),pdf->GetXaxis()->GetXmax());
   for(int i=0;i<=pdf->GetNbinsX();i++){
      if(i==0){
         PDF->SetBinContent(i, pdf->GetBinContent(i) );
      }else{
         PDF->SetBinContent(i, pdf->GetBinContent(i)+PDF->GetBinContent(i-1) );
      }
   }
   PDF->Scale(1.0/PDF->GetBinContent(PDF->GetNbinsX()));
   return PDF;
}

double GetRandValue(TH1D* PDF){
   int randNumber = rand();
   double uniform = randNumber / (double)RAND_MAX;
   for(int i=1;i<=PDF->GetNbinsX();i++){
      if(PDF->GetBinContent(i)>uniform){
         return PDF->GetXaxis()->GetBinUpEdge(i);
//         return PDF->GetXaxis()->GetBinUpEdge(i-1)+(rand()/(double)RAND_MAX)*PDF->GetXaxis()->GetBinWidth(i-1);
      }
   }
   return PDF->GetXaxis()->GetBinLowEdge(PDF->GetNbinsX());
}


void InitHistos(){
   for(unsigned int m=0;m<MCsample.size();m++){
      stPlots tmp;
      MCPlots.push_back(tmp);
   }

   for(unsigned int s=0;s<signals.size();s++){
      stPlots tmp;
      SignPlots.push_back(tmp);
   }
   HistoFile->cd();

   HCuts_Pt  = new TH1D("HCuts_Pt" ,"HCuts_Pt" ,CutPt.size(),0,CutPt.size());
   HCuts_I   = new TH1D("HCuts_I"  ,"HCuts_I"  ,CutPt.size(),0,CutPt.size());
   HCuts_TOF = new TH1D("HCuts_TOF","HCuts_TOF",CutPt.size(),0,CutPt.size());
   for(unsigned int i=0;i<CutPt.size();i++){  HCuts_Pt->Fill(i,CutPt[i]);     HCuts_I->Fill(i,CutI[i]);    HCuts_TOF->Fill(i,CutTOF[i]);   }
   /*
   if(DataFileName.size() || MCsample.size()){
      H_A = new TH1D("H_A" ,"H_A" ,CutPt.size(),0,CutPt.size());
      H_B = new TH1D("H_B" ,"H_B" ,CutPt.size(),0,CutPt.size());
      H_C = new TH1D("H_C" ,"H_C" ,CutPt.size(),0,CutPt.size());
      H_D = new TH1D("H_D" ,"H_D" ,CutPt.size(),0,CutPt.size());
      H_E = new TH1D("H_E" ,"H_E" ,CutPt.size(),0,CutPt.size());
      H_F = new TH1D("H_F" ,"H_F" ,CutPt.size(),0,CutPt.size());
      H_G = new TH1D("H_G" ,"H_G" ,CutPt.size(),0,CutPt.size());
      H_H = new TH1D("H_H" ,"H_H" ,CutPt.size(),0,CutPt.size());
      H_P = new TH1D("H_P" ,"H_P" ,CutPt.size(),0,CutPt.size());

      char Name   [1024];
      sprintf(Name,"Is");
      Hist_Is         = new TH1D(Name,Name, 200,0,dEdxS_UpLim);
      Hist_Is->Sumw2(); 

      sprintf(Name,"Pt");
      Hist_Pt       = new TH1D(Name,Name,200,0,PtHistoUpperBound);
      Hist_Pt->Sumw2();

      sprintf(Name,"TOF");
      Hist_TOF       = new TH1D(Name,Name,200,-10,20);
      Hist_TOF->Sumw2();

      sprintf(Name,"Pred_Mass");
      Pred_Mass = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),MassNBins,0,MassHistoUpperBound);
      Pred_Mass->Sumw2();

      sprintf(Name,"Pred_MassTOF");
      Pred_MassTOF = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(), MassNBins,0,MassHistoUpperBound);
      Pred_MassTOF->Sumw2();

      sprintf(Name,"Pred_MassComb");
      Pred_MassComb = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),MassNBins,0,MassHistoUpperBound);
      Pred_MassComb->Sumw2();

      sprintf(Name,"Pred_I");
      Pred_I  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   200,GlobalMinIm,dEdxM_UpLim);
      Pred_I->Sumw2();

      sprintf(Name,"Pred_EtaB");
      Pred_EtaB  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   50,-3,3);
      Pred_EtaB->Sumw2();

      sprintf(Name,"Pred_EtaS");
      Pred_EtaS  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   50,-3,3);
      Pred_EtaS->Sumw2();

      sprintf(Name,"Pred_EtaS2");
      Pred_EtaS2  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   50,-3,3);
      Pred_EtaS2->Sumw2();


      sprintf(Name,"Pred_EtaP");
      Pred_EtaP  = new TH3D(Name,Name,CutPt.size(),0,CutPt.size(),   50, -3, 3, 200,GlobalMinPt,PtHistoUpperBound);
      Pred_EtaP->Sumw2();

      sprintf(Name,"Pred_TOF");
      Pred_TOF  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   200,GlobalMinTOF,5);
      Pred_TOF->Sumw2();


      sprintf(Name,"RegionD_I");
      RegionD_I  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   200,GlobalMinIm,dEdxM_UpLim);
      RegionD_I->Sumw2();

      sprintf(Name,"RegionD_P");
      RegionD_P  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   200,GlobalMinPt,PtHistoUpperBound);
      RegionD_P->Sumw2();

      sprintf(Name,"RegionD_TOF");
      RegionD_TOF  = new TH2D(Name,Name,CutPt.size(),0,CutPt.size(),   200,GlobalMinTOF,5);
      RegionD_TOF->Sumw2();
   } 
   */
}


double DistToHSCP (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest){
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

   double RMin = 9999; IndexOfClosest=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000)continue;    

      double dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
      if(dR<RMin){RMin=dR;IndexOfClosest=g;}
   }
   return RMin;
}

double GetSampleWeight(const double& IntegratedLuminosityInPb, const double& IntegratedLuminosityInPbBeforeTriggerChange, const double& CrossSection, const double& MCEvents, int period){
  double Weight = 1.0;
  if(IntegratedLuminosityInPb>=IntegratedLuminosityInPbBeforeTriggerChange && IntegratedLuminosityInPb>0){
    double NMCEvents = MCEvents;
    //if(MaxEntry>0)NMCEvents=std::min(MCEvents,(double)MaxEntry);
    if      (period==0)Weight = (CrossSection * IntegratedLuminosityInPbBeforeTriggerChange) / NMCEvents;
    else if (period==1)Weight = (CrossSection * (IntegratedLuminosityInPb-IntegratedLuminosityInPbBeforeTriggerChange)) / NMCEvents;
  }
  return Weight;
}


double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent){
  double Weight = 1.0;
   unsigned long InitNumberOfEvents = GetInitialNumberOfMCEvent(fileNames); 
   double SampleEquivalentLumi = InitNumberOfEvents / XSection;
   if(MaxEvent<0)MaxEvent=SampleSize;
   printf("GetSampleWeight MC: IntLumi = %6.2E  SampleLumi = %6.2E --> EventWeight = %6.2E --> ",IntegratedLuminosityInPb,SampleEquivalentLumi, IntegratedLuminosityInPb/SampleEquivalentLumi);
//   printf("Sample NEvent = %6.2E   SampleEventUsed = %6.2E --> Weight Rescale = %6.2E\n",SampleSize, MaxEvent, SampleSize/MaxEvent);
   Weight = (IntegratedLuminosityInPb/SampleEquivalentLumi) * (SampleSize/MaxEvent);
   printf("FinalWeight = %6.2f\n",Weight);
   return Weight;
}

double GetPUWeight(const fwlite::ChainEvent& ev, const bool& Iss4pileup, double &PUSystFactor){
   //get pile up weight for this event
   fwlite::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   PupInfo.getByLabel(ev, "addPileupInfo");
   if(!PupInfo.isValid()){printf("PileupSummaryInfo Collection NotFound\n");return 1.0;}
   double PUWeight_thisevent=1;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int npv = -1;
   if(Iss4pileup){
      float sum_nvtx = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         npv = PVI->getPU_NumInteractions();
         sum_nvtx += float(npv);
      }
      float ave_nvtx = sum_nvtx/3.;
      PUWeight_thisevent = LumiWeightsMC_.weight3BX( ave_nvtx );
      PUSystFactor = PShift_.ShiftWeight( ave_nvtx );
   }else{
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      PUWeight_thisevent = LumiWeightsMC_.weight( npv );
      PUSystFactor = PShift_.ShiftWeight( npv );
   }


   return PUWeight_thisevent;
}


int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl){
   int toReturn = 0;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000)continue;
      if(AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324)continue; //Skip neutral gluino RHadrons
      if(AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333)continue;  //skip neutral stop RHadrons
      toReturn++;
   }
   return toReturn;
}


void  GetGenHSCPBeta (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, bool onlyCharged){
   beta1=-1; beta2=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000)continue;
      if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))continue; //Skip neutral gluino RHadrons
      if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))continue;  //skip neutral stop RHadrons
      if(beta1<0){beta1=genColl[g].p()/genColl[g].energy();}else if(beta2<0){beta2=genColl[g].p()/genColl[g].energy();return;}
   }
}

double RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge)
{
   double newInvPt = 1/pt+0.000236-0.000135*pow(eta,2)+charge*0.000282*TMath::Sin(phi-1.337);
   return 1/newInvPt;
}

unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames)
{
   unsigned long Total = 0;
   fwlite::ChainEvent tree(fileNames);

   for(unsigned int f=0;f<fileNames.size();f++){
     TFile *file;
     size_t place=fileNames[f].find("dcache");
     if(place!=string::npos) {
       string name=fileNames[f];
       name.replace(place, 7, "dcap://cmsgridftp.fnal.gov:24125");
       file = new TDCacheFile (name.c_str());
     }
     else file = new TFile (fileNames[f].c_str());
      fwlite::LuminosityBlock ls( file );
      for(ls.toBegin(); !ls.atEnd(); ++ls){
         fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
         nEventsTotalCounter.getByLabel(ls,"nEventsBefSkim");
         if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
         Total+= nEventsTotalCounter->value;
      }
   }
   return Total;
}

double SegSep(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double& minPhi, double& minEta) {
  reco::MuonRef muon = hscp.muonRef(); if(muon.isNull()) {cout << "No muon" << endl; return 10;}
  reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull()) {cout << "No track " << endl; return 10;}

  fwlite::Handle<MuonSegmentCollection> SegCollHandle;
  SegCollHandle.getByLabel(ev, "MuonSegmentProducer");
  if(!SegCollHandle.isValid()){printf("Segment Collection Not Found\n");return 0;}
  MuonSegmentCollection SegCollection = *SegCollHandle;

  double minDr=10;

  //Look for segment on opposite side of detector from track
  for (MuonSegmentCollection::const_iterator segment = SegCollection.begin(); segment!=SegCollection.end();++segment) {  
    GlobalPoint gp = segment->getGP();

    double eta_seg = gp.eta();
    double phi_seg = gp.phi();
    double eta_hscp = -1*track->eta();
    double phi_hscp= track->phi();
    //Flip phi to opposite side of detector
    if(phi_hscp<0) phi_hscp+=pi;
    else phi_hscp-=pi;

    double dEta=eta_seg-eta_hscp;
    double dPhi=phi_seg-phi_hscp;
    if(dPhi>pi) dPhi=-2*pi+dPhi;
    else if(dPhi<-1*pi) dPhi=2*pi+dPhi;

    if(fabs(dEta)<fabs(minEta) && fabs(dPhi)<(pi-0.5)) {
      minEta=dEta;
    }
    if(fabs(dPhi)<fabs(minPhi)) {
      minPhi=dPhi;
    }

    double dR=sqrt(dEta*dEta+dPhi*dPhi);
    if(dR<minDr) minDr=dR;
  }
  return minDr;
}
