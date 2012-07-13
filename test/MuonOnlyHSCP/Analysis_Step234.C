#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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

namespace reco    { class Vertex; class GenParticle; class MuonTimeExtra; class PFMET;}
namespace susybsm { class HSCParticle; class MuonSegment; class HSCPIsolation;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag; class LumiReWeighting;}
namespace reweight{class PoissonMeanShifter;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "DataFormats/FWLite/interface/Record.h"
#include "DataFormats/FWLite/interface/EventSetup.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

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
void  GetGenHSCPBeta   (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, double& pt1, double& pt2, bool onlyCharged=true);
bool   PassPreselection(const susybsm::HSCParticle& hscp, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, int& DzType, bool Control=false, const double& GenBeta=-1, const double& GenPt=-1, const double& GenCharge=-1);
bool   PassTkPreselection(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev);
bool PassTrigger      (const fwlite::ChainEvent& ev);
double GetPUWeight(const fwlite::ChainEvent& ev, const bool& Iss4pileup, double &PUSystFactor);
double GetSampleWeight(const double& IntegratedLuminosityInPb=-1, const double& CrossSection=0, const double& MCEvents=0);
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent);
unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames);
double SegSep(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double& minPhi, double& minEta);
double DistToTrigger (const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev);
double Zed(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev);
double RescalePt(double pt);
/////////////////////////// VARIABLE DECLARATION /////////////////////////////

TFile* HistoFile;

TH1D*  HCuts_Pt;
TH1D*  HCuts_TOF;
TH1D*  HCuts_TOFSyst;

std::vector<double>  CutPt ;
std::vector<double>  CutTOF;
std::vector<double>  CutTOFSyst;

std::vector<stSignal> signals;
std::vector<string>   DataFileName;

stPlots              DataPlots;  
stPlots              DataPlotsNoTrack;
stPlots              DataPlotsTrack;
stPlots              DataPlotsControl;
std::vector<stPlots> SignPlots; 

//for initializing PileUpReweighting utility.
const   float TrueDist2011_f[35] = {0.00285942, 0.0125603, 0.0299631, 0.051313, 0.0709713, 0.0847864, 0.0914627, 0.0919255, 0.0879994, 0.0814127, 0.0733995, 0.0647191, 0.0558327, 0.0470663, 0.0386988, 0.0309811, 0.0241175, 0.018241, 0.0133997, 0.00956071, 0.00662814, 0.00446735, 0.00292946, 0.00187057, 0.00116414, 0.000706805, 0.000419059, 0.000242856, 0.0001377, 7.64582e-05, 4.16101e-05, 2.22135e-05, 1.16416e-05, 5.9937e-06, 5.95542e-06};//from 2011 Full dataset

const   float Pileup_MC[35]= {1.45346E-01, 6.42802E-02, 6.95255E-02, 6.96747E-02, 6.92955E-02, 6.84997E-02, 6.69528E-02, 6.45515E-02, 6.09865E-02, 5.63323E-02, 5.07322E-02, 4.44681E-02, 3.79205E-02, 3.15131E-02, 2.54220E-02, 2.00184E-02, 1.53776E-02, 1.15387E-02, 8.47608E-03, 6.08715E-03, 4.28255E-03, 2.97185E-03, 2.01918E-03, 1.34490E-03, 8.81587E-04, 5.69954E-04, 3.61493E-04, 2.28692E-04, 1.40791E-04, 8.44606E-05, 5.10204E-05, 3.07802E-05, 1.81401E-05, 1.00201E-05, 5.80004E-06};

edm::LumiReWeighting LumiWeightsMC_;
std::vector< float > BgLumiMC; //MC                                           
std::vector< float > TrueDist2011;                                    
reweight::PoissonMeanShifter PShift_(0.6);//0.6 for upshift, -0.6 for downshift

TRandom3* RNG;

/////////////////////////// CODE PARAMETERS /////////////////////////////

void Analysis_Step234(string MODE_="COMPILE", string File="", float MinPt_=GlobalMinPt, float MaxEta_=GlobalMaxEta)
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

   GetSignalDefinition(signals, File);

   char Buffer[2048];   
   char Command[2048];
   DataFileName.clear();

   GlobalMaxEta = MaxEta_;
   GlobalMinPt    = MinPt_;

   for(double Pt =GlobalMinPt+30 ; Pt <=500;Pt+=30){
     for(double TOF=1.00; TOF<1.35;TOF+=0.01){
       CutPt .push_back(Pt); CutTOF.push_back(TOF);
     }
     for(double TOF=1.00; TOF>0.65;TOF-=0.01){
       CutTOFSyst.push_back(TOF);
     }
   }

   printf("%i Different Final Selection will be tested\n",(int)CutPt.size());

   //initialize LumiReWeighting
   for(int i=0; i<35; ++i)   BgLumiMC.push_back(Pileup_MC[i]);
   for(int i=0; i<35; ++i)    TrueDist2011.push_back(TrueDist2011_f[i]);
   LumiWeightsMC_ = edm::LumiReWeighting(BgLumiMC, TrueDist2011);

   sprintf(Buffer,"Results/"       );                                          sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sEta%02.0f/"  ,Buffer,10.0*GlobalMaxEta);                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sPtMin%02.0f/",Buffer,GlobalMinPt);                        sprintf(Command,"mkdir %s",Buffer); system(Command);

   time_t start = time(NULL);
   if(MODE=="ANALYSE_DATA"){
      signals.clear();  //Remove all signal samples
      GetInputFiles(DataFileName, "Data", File);
      HistoFile = new TFile((string(Buffer) + "/Histos_Data_" + File + ".root").c_str(),"RECREATE");
   }else if(MODE=="ANALYSE_SIGNAL"){
      DataFileName.clear();  //Remove all data files
      HistoFile = new TFile((string(Buffer) + "/Histos_" + File + ".root").c_str(),"RECREATE");
   }else if(MODE=="ANALYSE_COSMIC"){
      DataFileName.clear();  //Remove all data files
      signals.clear();  //Remove all signal samples
      GetInputFiles(DataFileName, "Cosmic", File);
      HistoFile = new TFile((string(Buffer) + "/Histos_Cosmic_" + File + ".root").c_str(),"RECREATE");
   }else{
      printf("You must select a MODE:\n");
      printf("MODE='ANALYSE_DATA'   : Will run the analysis on Data\n"); 
      printf("MODE='ANALYSE_SIGNAL' : Will run the analysis on Signal MC\n");
      printf("MODE='ANALYSE_COSMIC' : Will run the analysis on Cosmics\n");
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

bool PassTrigger(const fwlite::ChainEvent& ev)
{
  bool toReturn=false;
      edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");

      if(MODE=="ANALYSE_SIGNAL") {
	fwlite::Handle<reco::PFMETCollection> pfMETCollection;
	pfMETCollection.getByLabel(ev,"pfMet");

	if(tr.accept(tr.triggerIndex("HSCPPathSAMU")) && pfMETCollection->begin()->et()>60) toReturn = true;}
      else if(MODE!="ANALYSE_COSMIC") {if(tr.accept(tr.triggerIndex("HSCPHLTTriggerL2MuFilter")))toReturn = true;}
      else if(tr.accept(tr.triggerIndex("HSCPHLTTriggerCosmicFilter")))toReturn = true;

      return toReturn;
}

bool PassTkPreselection(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev) {
  edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");
  if(!tr.isValid())return false;

  if(!tr.accept("HSCPHLTTriggerMetDeDxFilter") && !tr.accept("HSCPHLTTriggerMuFilter"))return false;

  reco::TrackRef track = hscp.trackRef(); if(track.isNull())return false;
  if(fabs(track->eta())>GlobalMaxTkEta) return false;
  if(track->found()<GlobalMinTkNOH)return false;
  if(track->validFraction()<0.80)return false;
  if(track->hitPattern().numberOfValidPixelHits()<2)return false;
  if(track->qualityMask()<GlobalMinTkQual )return false;
  if(track->chi2()/track->ndof()>GlobalMaxTkChi2 )return false;
  if(track->pt()<GlobalMinTkPt)return false;

  fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
  dEdxSCollH.getByLabel(ev, dEdxS_Label.c_str());
  if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");return false;}
  fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
  dEdxMCollH.getByLabel(ev, dEdxM_Label.c_str());
  if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");return false;}
  const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
  const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());
  if(dedxSObj.numberOfMeasurements()<GlobalMinTkNOM)return false;
  if(dedxSObj.dEdx()<GlobalMinTkIs)return false;
  if(dedxMObj.dEdx()<GlobalMinTkIm)return false;

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
  if(v3d>GlobalMaxTkV3D )return false;

  fwlite::Handle<HSCPIsolationValueMap> IsolationH;
  IsolationH.getByLabel(ev, "HSCPIsolation03");
  if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
  const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();
  HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
  if(hscpIso.Get_TK_SumEt()>GlobalMaxTkTIsol)return false;
  double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
  if(EoP>GlobalMaxTkEIsol)return false;
  if((track->ptError()/track->pt())>GlobalMaxTkPterr)return false;

  return true;
}

bool PassPreselection(const susybsm::HSCParticle& hscp, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, int& DzType, bool Control, const double& GenBeta, const double& GenPt, const double &GenCharge)
{   

  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return false;  

   reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull())return false;
   reco::TrackRef innertrack = hscp.trackRef();
   bool isGlobal=(muon->isGlobalMuon() && muon->isTrackerMuon() && !innertrack.isNull());

   if(!tof) return false;

   //Make distributions without any cuts
   if(st) {
     st->BS_Pt_All->Fill(track->pt(), Event_Weight);
     if(tof->nDof()>0)st->BS_TOF_All->Fill(tof->inverseBeta(), Event_Weight);
     st->Total->Fill(0.0,Event_Weight);
     if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);
   }

   //Require track to match trigger object
   st->DistTrigger->Fill(DistToTrigger(hscp, ev),Event_Weight);

   //if(DistToTrigger(hscp, ev)>MaxDistTrigger) return false;
   st->TriggerMatch->Fill(0.0, Event_Weight);

   //Match to a SA track without vertex constraint for IP cuts
   fwlite::Handle< std::vector<reco::Track> > noVertexTrackCollHandle;
   noVertexTrackCollHandle.getByLabel(ev,"RefitSAMuons", "");
   //noVertexTrackCollHandle.getByLabel(ev,"RefitMTSAMuons", "");
   if(!noVertexTrackCollHandle.isValid()){
     //noVertexTrackCollHandle.getByLabel(ev,"refittedStandAloneMuons", "");
     //if(!noVertexTrackCollHandle.isValid()){printf("No Vertex Track Collection Not Found\n");return false;}
     noVertexTrackCollHandle.getByLabel(ev,"refittedStandAloneMuons", "");
     if(!noVertexTrackCollHandle.isValid()){printf("No Vertex Track Collection Not Found\n");return false;}
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
   if(fabs(track->eta())>GlobalMaxEta) return false;
   if(st) st->Eta->Fill(0.0, Event_Weight);

   //Cut on min pt
   if(track->pt()<GlobalMinPt)return false;
   if(st) st->MinPt->Fill(0.0, Event_Weight);

   //Cut on number of matched muon stations
   int count=track->hitPattern().dtStationsWithValidHits()+track->hitPattern().cscStationsWithValidHits();
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
   int ZedSegs=6;
   ZedSegs=Zed(hscp, ev);

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
     if(fabs(dz)>GlobalMaxDz) {
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
   if(fabs(dxy)>GlobalMaxDxy && !isGlobal) return false;
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

   //Cut on dz for SA only tracks but not if this for the control region
   if(fabs(dz)>GlobalMaxDz && !Control) return false;

   //Split into different dz regions, each different region used to predict cosmic background and find systematic
   if(Control) {
     if(fabs(dz)<GlobalMaxDz) DzType=0;
     else if(fabs(dz)<30) DzType=1;
     else if(fabs(dz)<50) DzType=2;
     else if(fabs(dz)<70) DzType=3;
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) DzType=4;
     if(fabs(dz)>CosmicMaxDz) DzType=5;
   }

   //Require control region cuts
   if(Control && (fabs(dz)<CosmicMinDz || fabs(dz)>CosmicMaxDz)) return false;
   
   //Fill a bunch of plots for tracks passing preselection
   if(st){
     st->Dz  ->Fill(0.0,Event_Weight);
     if(ZedSegs==0) st->BS_Dz_NoZed->Fill(dz,Event_Weight);
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

void Analysis_FillControlAndPredictionHist(const susybsm::HSCParticle& hscp, const reco::MuonTimeExtra* tof, stPlots* st=NULL, int DzType=-1){
	 reco::MuonRef muon = hscp.muonRef();
	 reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull()){ cout << "No track" << endl; return;}

	 st->H_DzCounts->Fill(DzType, Event_Weight);
	 if(fabs(track->eta())<DTRegion) st->H_DzCounts_DT->Fill(DzType, Event_Weight);
	 else st->H_DzCounts_CSC->Fill(DzType, Event_Weight);

         double MuonTOF = tof->inverseBeta();

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
	 }

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
	   //This loop is to determine the number of tracks in collision control region
	   //defined as having TOF less than 0.9 or between 0.9 and 1
	   //Background prediction is done only for tracks > 1, these are used to find systematic
	   if(MuonTOF>1) continue;
            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassTOFCut = MuonTOF<=0.9;

	    //DzType=-1 for tracks passing all the preselection cuts, other values are for tracks passing all cuts except for dz cut
	    //Different values of DzType are for different ranges of Dz values.  DzType==4 is the region used to predict the cosmic 
	    //background, the others are used to find a systematic on this prediction.  The others are only used later
	    if(DzType==-1 || DzType==4) {
            if(       PassTOFCut &&  PassPtCut){   //Region D
	      st->H_D_Low->Fill(CutIndex,                Event_Weight);
	      if(fabs(track->eta())<DTRegion) st->H_D_Cen_Low->Fill(CutIndex, Event_Weight);
	      else st->H_D_For_Low->Fill(CutIndex, Event_Weight);

            }else if( !PassTOFCut &&  PassPtCut){   //Region C
	      st->H_C_Low->Fill(CutIndex,                 Event_Weight);
	      if(fabs(track->eta())<DTRegion) st->H_C_Cen_Low->Fill(CutIndex, Event_Weight);
	      else st->H_C_For_Low->Fill(CutIndex, Event_Weight);                                                                      

            }else if( PassTOFCut && !PassPtCut){   //Region B
	      st->H_B_Low->Fill(CutIndex,                 Event_Weight);
	      if(fabs(track->eta())<DTRegion) st->H_B_Cen_Low->Fill(CutIndex, Event_Weight);
	      else st->H_B_For_Low->Fill(CutIndex, Event_Weight);

            }else if( !PassTOFCut && !PassPtCut){   //Region A
	      st->H_A_Low->Fill(CutIndex,                 Event_Weight);
	      if(fabs(track->eta())<DTRegion) st->H_A_Cen_Low->Fill(CutIndex, Event_Weight);
	      else st->H_A_For_Low->Fill(CutIndex, Event_Weight);
	    }
	    }
	 }

	 //Now loop again for main prediction region
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
           if(MuonTOF<1) continue;
            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassTOFCut = MuonTOF>=CutTOF[CutIndex];

	    //Again only for two main types of tracks
	    if(DzType==-1 || DzType==4) {
	      if(PassPtCut) st->Pt->Fill(CutIndex, Event_Weight);
	      if(PassTOFCut) st->TOFExclusive->Fill(CutIndex, Event_Weight);
	      if(PassPtCut && PassTOFCut) st->TOF->Fill(CutIndex, Event_Weight);

	      if(       PassTOFCut &&  PassPtCut){   //Region D
		st->Eta_RegionD->Fill(CutIndex,track->eta(), Event_Weight);
		st->H_D      ->Fill(CutIndex,                Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_D_Cen->Fill(CutIndex, Event_Weight);
		else st->H_D_For->Fill(CutIndex, Event_Weight);

	      }else if( !PassTOFCut &&  PassPtCut){   //Region C
		st->Eta_RegionC->Fill(CutIndex,track->eta(), Event_Weight);
		st->H_C->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_C_Cen->Fill(CutIndex, Event_Weight);
		else st->H_C_For->Fill(CutIndex, Event_Weight);                                         

	      }else if( PassTOFCut && !PassPtCut){   //Region B
		st->Eta_RegionB->Fill(CutIndex,track->eta(), Event_Weight);
		st->H_B     ->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_B_Cen->Fill(CutIndex, Event_Weight);
		else st->H_B_For->Fill(CutIndex, Event_Weight);

	      }else if( !PassTOFCut && !PassPtCut){   //Region A
		st->Eta_RegionA->Fill(CutIndex,track->eta(), Event_Weight);
		st->H_A     ->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_A_Cen->Fill(CutIndex, Event_Weight);
		else st->H_A_For->Fill(CutIndex, Event_Weight);
	      }
	    }

	    //Now fill plots for cosmic systematic allowing for background to be predicted again
	    if(DzType>-1) {
	      if(       PassTOFCut &&  PassPtCut){   //Region D
		st->H_D_Syst[DzType]->Fill(CutIndex,                Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_D_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
		else st->H_D_For_Syst[DzType]->Fill(CutIndex, Event_Weight);

	      }else if( !PassTOFCut &&  PassPtCut){   //Region C
		st->H_C_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_C_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
		else st->H_C_For_Syst[DzType]->Fill(CutIndex, Event_Weight);                                         

	      }else if( PassTOFCut && !PassPtCut){   //Region B
		st->H_B_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_B_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
		else st->H_B_For_Syst[DzType]->Fill(CutIndex, Event_Weight);

	      }else if( !PassTOFCut && !PassPtCut){   //Region A
		st->H_A_Syst[DzType]->Fill(CutIndex,                 Event_Weight);
		if(fabs(track->eta())<DTRegion) st->H_A_Cen_Syst[DzType]->Fill(CutIndex, Event_Weight);
		else st->H_A_For_Syst[DzType]->Fill(CutIndex, Event_Weight);
	      }
	    }
	 }
}

void Analysis_Step3(char* SavePath)
{
   printf("Step3: Building Mass Spectrum for B and S\n");

   int TreeStep;
   //////////////////////////////////////////////////     BUILD BACKGROUND MASS SPECTRUM

   if(MODE=="ANALYSE_COSMIC") {
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlots,"Cosmic", CutPt.size(), false);
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsTrack,"Cosmic_Track", CutPt.size());
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsNoTrack,"Cosmic_NoTrack", CutPt.size());
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsControl,"Cosmic_Control", CutPt.size());
   }
   else if(MODE=="ANALYSE_DATA") {
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsTrack,"Data_Track", CutPt.size());
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsNoTrack,"Data_NoTrack", CutPt.size());
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlots,"Data", CutPt.size(), false);
     if(DataFileName.size())stPlots_Init(HistoFile, DataPlotsControl,"Data_Control", CutPt.size());
   }

   HistoFile->cd();
   fwlite::ChainEvent treeD(DataFileName);

   double SampleWeight = GetSampleWeight(-1);
   Event_Weight = SampleWeight;
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Building Mass Spectrum for D :");
   TreeStep = treeD.size()/50;if(TreeStep==0)TreeStep=1;


   for(Long64_t ientry=0;ientry<treeD.size();ientry++){
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
      treeD.to(ientry);

      DataPlots.TotalE->Fill(0.0,Event_Weight);  
      DataPlotsNoTrack.TotalE->Fill(0.0,Event_Weight); DataPlotsTrack.TotalE->Fill(0.0,Event_Weight);

      if(!PassTrigger(treeD) )continue;

      DataPlots.TotalTE->Fill(0.0,Event_Weight);
      DataPlotsNoTrack.TotalTE->Fill(0.0,Event_Weight); DataPlotsTrack.TotalTE->Fill(0.0,Event_Weight);

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeD,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(treeD, "muontiming","combined");
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

      fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
      TOFDTCollH.getByLabel(treeD, "muontiming","dt");
      if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");return;}

      fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
      TOFCSCCollH.getByLabel(treeD, "muontiming","csc");
      if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");return;}

      for(unsigned int c=0;c<hscpColl.size();c++){

         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::TrackRef innertrack = hscp.trackRef();
	 reco::MuonRef muon = hscp.muonRef();

         if(muon.isNull())continue;
	 reco::TrackRef   SAtrack = muon->standAloneMuon(); if(SAtrack.isNull()) continue;

         if(PassTkPreselection(hscp, treeD)) continue;

	 bool isGlobal=(muon->isGlobalMuon() && muon->isTrackerMuon() && !innertrack.isNull());

         if(isGlobal) DataPlotsTrack.Reconstructed->Fill(0.0, Event_Weight);
         if(!isGlobal) DataPlotsNoTrack.Reconstructed->Fill(0.0, Event_Weight);
	 DataPlots.Reconstructed->Fill(0.0, Event_Weight);

         const reco::MuonTimeExtra* tof = NULL;
         const reco::MuonTimeExtra* dttof = NULL;
         const reco::MuonTimeExtra* csctof = NULL;
         tof  = &TOFCollH->get(hscp.muonRef().key()); dttof = &TOFDTCollH->get(hscp.muonRef().key());  csctof = &TOFCSCCollH->get(hscp.muonRef().key());
	 if(!tof) continue;

         if(isGlobal) DataPlotsTrack.tofFound->Fill(0.0, Event_Weight);
         if(!isGlobal) DataPlotsNoTrack.tofFound->Fill(0.0, Event_Weight);
         DataPlots.tofFound->Fill(0.0, Event_Weight);

	 //Dz type is used to determine in which control region in the dz distribution a track falls
	 //The different regions are used to predict the cosmic background and its uncertainty
	 int DzType=-1;
         if(isGlobal) PassPreselection(hscp, tof, dttof, csctof, treeD, &DataPlotsTrack, DzType);
	 else PassPreselection(hscp, tof, dttof, csctof, treeD, &DataPlotsNoTrack, DzType);

         if(!isGlobal) PassPreselection(hscp, tof, dttof, csctof, treeD, &DataPlotsControl, DzType, true);

         if(DzType>-1) Analysis_FillControlAndPredictionHist(hscp, tof, &DataPlotsControl, DzType);

	 if(!PassPreselection(hscp, tof, dttof, csctof, treeD, &DataPlots, DzType)) continue;

         if(isGlobal) Analysis_FillControlAndPredictionHist(hscp, tof, &DataPlotsTrack);
	 else Analysis_FillControlAndPredictionHist(hscp, tof, &DataPlotsNoTrack);
	 Analysis_FillControlAndPredictionHist(hscp, tof, &DataPlots);

	 stPlots_FillTree(DataPlots, treeD.eventAuxiliary().run(),treeD.eventAuxiliary().event(),c,  SAtrack->pt(), tof->inverseBeta(), -1);
      } // end of Track Loop
   }// end of Event Loop

   printf("\n");
   if(DataFileName.size())stPlots_Clear(DataPlots, true, false);
   if(DataFileName.size())stPlots_Clear(DataPlotsNoTrack, true);
   if(DataFileName.size())stPlots_Clear(DataPlotsTrack, true);
   if(DataFileName.size())stPlots_Clear(DataPlotsControl, true);


   //////////////////////////////////////////////////     BUILD SIGNAL MASS SPECTRUM
   RNG = new TRandom3();
   for(unsigned int s=0;s<signals.size();s++){
      stPlots_Init(HistoFile,SignPlots[4*s+0],signals[s].Name       , CutPt.size());
      stPlots_Init(HistoFile,SignPlots[4*s+1],signals[s].Name+"_NC0", CutPt.size());//, true);
      stPlots_Init(HistoFile,SignPlots[4*s+2],signals[s].Name+"_NC1", CutPt.size());//, true);
      stPlots_Init(HistoFile,SignPlots[4*s+3],signals[s].Name+"_NC2", CutPt.size());//, true);

      bool* HSCPTk         = new bool[CutPt.size()];
      bool* HSCPTk_SystTOF = new bool[CutPt.size()];
      bool* HSCPTk_SystPt  = new bool[CutPt.size()];

      printf("Progressing Bar                                    :0%%       20%%       40%%       60%%       80%%       100%%\n");

      std::vector<string> SignFileName;
      GetInputFiles(SignFileName, signals[s].Name);

      fwlite::ChainEvent treeS(SignFileName);

      printf("Running over sample %10s :",signals[s].Name.c_str());

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

      double SampleWeight = GetSampleWeight(IntegratedLuminosity,signals[s].XSec,NMCevents);
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

	 //Fill plots on gen beta distribution
         double HSCPGenBeta1, HSCPGenBeta2;
         double HSCPGenPt1, HSCPGenPt2;
         GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,HSCPGenPt1,HSCPGenPt2,true);
         if(HSCPGenBeta1>=0)SignPlots[4*s].Beta_GenCharged->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SignPlots[4*s].Beta_GenCharged->Fill(HSCPGenBeta2, Event_Weight);
         GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,HSCPGenPt1,HSCPGenPt2,false);
         if(HSCPGenBeta1>=0) {SignPlots[4*s].Beta_Gen->Fill(HSCPGenBeta1, Event_Weight); SignPlots[4*s].Pt_Gen->Fill(HSCPGenPt1, Event_Weight);}
	 if(HSCPGenBeta2>=0) {SignPlots[4*s].Beta_Gen->Fill(HSCPGenBeta2, Event_Weight); SignPlots[4*s].Pt_Gen->Fill(HSCPGenPt2, Event_Weight);}

         SignPlots[4*s]               .TotalE   ->Fill(0.0,Event_Weight);
         SignPlots[4*s+NChargedHSCP+1].TotalE   ->Fill(0.0,Event_Weight);
         SignPlots[4*s]               .TotalEPU ->Fill(0.0,Event_Weight*PUSystFactor);
         SignPlots[4*s+NChargedHSCP+1].TotalEPU ->Fill(0.0,Event_Weight*PUSystFactor);

         if(!PassTrigger(treeS) )continue;
         SignPlots[4*s]               .TotalTE->Fill(0.0,Event_Weight);
         SignPlots[4*s+NChargedHSCP+1].TotalTE->Fill(0.0,Event_Weight);

         if(HSCPGenBeta1>=0)SignPlots[4*s].Beta_Triggered->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SignPlots[4*s].Beta_Triggered->Fill(HSCPGenBeta2, Event_Weight);

         fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
         hscpCollHandle.getByLabel(treeS,"HSCParticleProducer");
         if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
         const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

         fwlite::Handle<MuonTimeExtraMap> TOFCollH;
         TOFCollH.getByLabel(treeS, "muontiming","combined");
         if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
         TOFDTCollH.getByLabel(treeS, "muontiming","dt");
         if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

         fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
         TOFCSCCollH.getByLabel(treeS, "muontiming","csc");
         if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");continue;}

         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  
	   HSCPTk[CutIndex] = false;
           HSCPTk_SystTOF[CutIndex] = false;
           HSCPTk_SystPt[CutIndex] = false;
	 }

         for(unsigned int c=0;c<hscpColl.size();c++){
	   susybsm::HSCParticle hscp  = hscpColl[c];
           int ClosestGen;

	   reco::TrackRef track = hscp.trackRef();
	   reco::MuonRef muon = hscp.muonRef();

	   if(muon.isNull())continue;
	   reco::TrackRef   SAtrack = muon->standAloneMuon(); if(SAtrack.isNull())continue;
	   if(PassTkPreselection(hscp, treeS)) continue;

	   SignPlots[4*s+NChargedHSCP+1].DistToGen->Fill(DistToHSCP(hscp, genColl, ClosestGen),Event_Weight);
	   SignPlots[4*s].DistToGen->Fill(DistToHSCP(hscp, genColl, ClosestGen),Event_Weight);

	   //Require it matches a generator level HSCP
           if(DistToHSCP(hscp, genColl, ClosestGen)>0.04)continue;
           SignPlots[4*s]               .Reconstructed->Fill(0.0,Event_Weight);
           SignPlots[4*s+NChargedHSCP+1].Reconstructed->Fill(0.0,Event_Weight);


	   const reco::MuonTimeExtra* tof = NULL;
	   const reco::MuonTimeExtra* dttof = NULL;
	   const reco::MuonTimeExtra* csctof = NULL;
	   tof  = &TOFCollH->get(hscp.muonRef().key()); dttof = &TOFDTCollH->get(hscp.muonRef().key());  csctof = &TOFCSCCollH->get(hscp.muonRef().key());
	   if(!tof) continue;
           SignPlots[4*s]               .tofFound->Fill(0.0,Event_Weight);
           SignPlots[4*s+NChargedHSCP+1].tofFound->Fill(0.0,Event_Weight);

	   double MuonTOF = GlobalMinTOF;
	   MuonTOF = tof->inverseBeta();

	   double TRescale = -0.02; // added to the 1/beta value
	   if(tof) if(csctof->nDof()==0) TRescale = -0.003; //If any CSC hits use the larger error coming from CSC

	   //Require it pass preselection.  DzType here is just a dummy need to be passed to the function
	   int DzType=0;
	   PassPreselection(    hscp, tof, dttof, csctof, treeS,&SignPlots[4*s+NChargedHSCP+1],DzType,false,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), genColl[ClosestGen].pt(), genColl[ClosestGen].charge());
	   if(!PassPreselection(hscp, tof, dttof, csctof, treeS,&SignPlots[4*s               ],DzType,false,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), genColl[ClosestGen].pt(), genColl[ClosestGen].charge()))continue;         

	   //Loop over all the different selections
            for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
	      //Calculate TOF Systematic
              if(SAtrack->pt()>=CutPt[CutIndex] && MuonTOF+TRescale>=CutTOF[CutIndex]) {
		HSCPTk_SystTOF[CutIndex] = true;
	      }

	      //Calculate Pt Systematic
              if(RescalePt(SAtrack->pt())>=CutPt[CutIndex] && MuonTOF>=CutTOF[CutIndex]) {
	        HSCPTk_SystPt[CutIndex] = true;
	      }

	      //Now find the actual number passing
	      if(SAtrack->pt()<CutPt[CutIndex]) continue;

	      SignPlots[4*s               ].Beta_SelectedP->Fill(CutIndex,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), Event_Weight);
              SignPlots[4*s               ].Pt  ->Fill(CutIndex,Event_Weight);
              SignPlots[4*s+NChargedHSCP+1].Beta_SelectedP->Fill(CutIndex,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), Event_Weight);
              SignPlots[4*s+NChargedHSCP+1].Pt  ->Fill(CutIndex,Event_Weight);

	      if(MuonTOF<CutTOF[CutIndex])continue;
              SignPlots[4*s               ].Beta_SelectedT->Fill(CutIndex,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), Event_Weight);
              SignPlots[4*s               ].TOF  ->Fill(CutIndex,Event_Weight);
              SignPlots[4*s+NChargedHSCP+1].Beta_SelectedT->Fill(CutIndex,genColl[ClosestGen].p()/genColl[ClosestGen].energy(), Event_Weight);
              SignPlots[4*s+NChargedHSCP+1].TOF  ->Fill(CutIndex,Event_Weight);

	      HSCPTk[CutIndex] = true;
            } //end of Cut loop

         } // end of Track Loop 

	 //Loop again over the cuts to determine number of events with at least one track passing
         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
	   if(HSCPTk[CutIndex]){
	     SignPlots[4*s               ].HSCPE->Fill(CutIndex,Event_Weight);
	     SignPlots[4*s+NChargedHSCP+1].HSCPE->Fill(CutIndex,Event_Weight);}
           if(HSCPTk_SystTOF[CutIndex]){
	     SignPlots[4*s               ].HSCPE_SystT->Fill(CutIndex,Event_Weight);
             SignPlots[4*s+NChargedHSCP+1].HSCPE_SystT->Fill(CutIndex,Event_Weight);}
           if(HSCPTk_SystPt[CutIndex]){
	     SignPlots[4*s               ].HSCPE_SystP->Fill(CutIndex,Event_Weight);
             SignPlots[4*s+NChargedHSCP+1].HSCPE_SystP->Fill(CutIndex,Event_Weight);}
	 }
      }// end of Event Loop
      
      printf("\n");
      delete [] HSCPTk;
      delete [] HSCPTk_SystTOF;
      delete [] HSCPTk_SystPt;

      stPlots_Clear(SignPlots[4*s+0], true);
      stPlots_Clear(SignPlots[4*s+1], true);
      stPlots_Clear(SignPlots[4*s+2], true);
      stPlots_Clear(SignPlots[4*s+3], true);
   }// end of signal Type loop
   delete RNG;
}

void InitHistos(){
   for(unsigned int s=0;s<signals.size();s++){
   for(int NC=0;NC<4;NC++){
      stPlots tmp;
      if(NC==0){
      }else{
         char buffer[256];sprintf(buffer,"_NC%i",NC-1);
      }
      SignPlots.push_back(tmp);
   }}
   HistoFile->cd();

   HCuts_Pt  = new TH1D("HCuts_Pt" ,"HCuts_Pt" ,CutPt.size(),0,CutPt.size());
   HCuts_TOF = new TH1D("HCuts_TOF","HCuts_TOF",CutPt.size(),0,CutPt.size());
   HCuts_TOFSyst = new TH1D("HCuts_TOFSyst","HCuts_TOFSyst",CutPt.size(),0,CutPt.size());
   for(unsigned int i=0;i<CutPt.size();i++){  HCuts_Pt->Fill(i,CutPt[i]);   HCuts_TOF->Fill(i,CutTOF[i]); HCuts_TOFSyst->Fill(i,CutTOFSyst[i]);}
}

double DistToHSCP (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest){
  reco::MuonRef   muon;
  muon = hscp.muonRef();
  if(muon.isNull())return false;

  reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull()) return false;

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

double GetSampleWeight(const double& IntegratedLuminosityInPb, const double& CrossSection, const double& MCEvents){
  double Weight = 1.0;
  if(IntegratedLuminosityInPb>0){
    double NMCEvents = MCEvents;
    if(MaxEntry>0)NMCEvents=std::min(MCEvents,(double)MaxEntry);
    Weight = (CrossSection * IntegratedLuminosityInPb) / NMCEvents;
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


void  GetGenHSCPBeta (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, double& pt1, double& pt2, bool onlyCharged){
   beta1=-1; beta2=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000)continue;
      if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))continue; //Skip neutral gluino RHadrons
      if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))continue;  //skip neutral stop RHadrons
      if(beta1<0){beta1=genColl[g].p()/genColl[g].energy(); pt1=genColl[g].pt();}
      else if(beta2<0){beta2=genColl[g].p()/genColl[g].energy();pt2=genColl[g].pt();return;}
   }
}

unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames)
{
   unsigned long Total = 0;
   fwlite::ChainEvent tree(fileNames);

   for(unsigned int f=0;f<fileNames.size();f++){
      TFile file(fileNames[f].c_str() );
      fwlite::LuminosityBlock ls( &file );
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

//Function finds the number of DT segments with Zed projection
//If track has any CSC hits it returns six
//Code mostly lifted from RecoMuon/TrackingTools/src/MuonSegmentMatcher.cc
double Zed(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev) {
   reco::MuonRef muon = hscp.muonRef();
   reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull())return 0;

   fwlite::Handle<DTRecSegment4DCollection> dtRecHits;
   //dtRecHits.getByLabel(ev, "dt4DSegmentsMT");
   dtRecHits.getByLabel(ev, "dt4DSegments");

   if(!dtRecHits.isValid()){printf("DT Segment Collection NotFound\n");return 0;}

   double matchRatioZ=0;
   double matchRatioPhi=0;

   double totalDTSegs=0;
   double WithZed=0;

   bool segments=true;

   for(trackingRecHit_iterator hit = track->recHitsBegin(); hit != track->recHitsEnd(); ++hit) {
     if ( !(*hit)->isValid()) continue;
     if ( (*hit)->geographicalId().det() != DetId::Muon ) continue;
     if ( (*hit)->geographicalId().subdetId() != MuonSubdetId::CSC ) continue;
     return 6;
   }


   TrackingRecHitRefVector dtHits;

   for(trackingRecHit_iterator hit = track->recHitsBegin(); hit != track->recHitsEnd(); ++hit) {
     if ( !(*hit)->isValid()) continue; 
     if ( (*hit)->geographicalId().det() != DetId::Muon ) continue; 
     if ( (*hit)->geographicalId().subdetId() != MuonSubdetId::DT ) continue; 
     if (!(*hit)->isValid()) continue; 
     if ((*hit)->recHits().size()>1) segments = true;
     dtHits.push_back(*hit);
   }

   if(MODE=="ANALYSE_SIGNAL") segments=false;
   for (DTRecSegment4DCollection::const_iterator rechit = dtRecHits->begin(); rechit!=dtRecHits->end();++rechit) {
     if ( !rechit->isValid()) continue; 
     LocalPoint pointLocal = rechit->localPosition();

     if (segments) {
       // Loop over muon recHits
       for(trackingRecHit_iterator hit = dtHits.begin(); hit != dtHits.end(); ++hit) {
	 if ( !(*hit)->isValid()) continue; 
	 // Pick the one in the same DT Chamber as the muon
	 DetId idT = (*hit)->geographicalId();
	 if(!(rechit->geographicalId().rawId()==idT.rawId())) continue; 

	 // and compare the local positions
	 LocalPoint segLocal = (*hit)->localPosition();
	 if ((fabs(pointLocal.x()-segLocal.x())<1) && 
	     (fabs(pointLocal.y()-segLocal.y())<1)) 
	   if(rechit->hasZed()) WithZed++;
       }
     }

     else{
     double nhitsPhi = 0;
     double nhitsZ = 0;

     if(rechit->hasZed()) {
       double countMuonDTHits = 0;
       double countAgreeingHits=0;

       const DTRecSegment2D* segmZ;
       segmZ = dynamic_cast<const DTRecSegment2D*>(rechit->zSegment());
       nhitsZ = segmZ->recHits().size(); 

       const vector<DTRecHit1D> hits1d = segmZ->specificRecHits();
       DTChamberId chamberSegIdT((segmZ->geographicalId()).rawId());

       // Loop over muon recHits
       for(trackingRecHit_iterator hit = track->recHitsBegin(); hit != track->recHitsEnd(); ++hit) {
	 if ( !(*hit)->isValid()) continue; 
	 if ( (*hit)->geographicalId().det() != DetId::Muon ) continue; 
	 if ( (*hit)->geographicalId().subdetId() != MuonSubdetId::DT ) continue; 
	 if (!(*hit)->isValid()) continue;

	 if ( !(*hit)->isValid()) continue; 
	 
	 DetId idT = (*hit)->geographicalId();
	 DTChamberId dtDetIdHitT(idT.rawId());
	 DTSuperLayerId dtDetLayerIdHitT(idT.rawId());

	 LocalPoint  pointLocal = (*hit)->localPosition();
	 
	 if ((chamberSegIdT==dtDetIdHitT) && (dtDetLayerIdHitT.superlayer()==2)) countMuonDTHits++;

	 for (vector<DTRecHit1D>::const_iterator hiti=hits1d.begin(); hiti!=hits1d.end(); hiti++) {

	   if ( !hiti->isValid()) continue; 

	   // Pick the one in the same DT Layer as the 1D hit
	   if(!(hiti->geographicalId().rawId()==idT.rawId())) continue; 

	   // and compare the local positions
	   LocalPoint segLocal = hiti->localPosition();
	   //  cout << "Zed Segment Point = "<<pointLocal<<"    Muon Point = "<<segLocal<<"  Dist:  "
	   //       << (fabs(pointLocal.x()-segLocal.x()))+(fabs(pointLocal.y()-segLocal.y()))<< endl;
	   if ((fabs(pointLocal.x()-segLocal.x())<1) && 
	       (fabs(pointLocal.y()-segLocal.y())<1)) 
	     countAgreeingHits++;
	 } //End Segment Hit Iteration
       } //End Muon Hit Iteration
       
       matchRatioZ = countMuonDTHits == 0 ? 0 : countAgreeingHits/countMuonDTHits;
       if (nhitsZ)
	 if (countAgreeingHits/nhitsZ>matchRatioZ) matchRatioZ=countAgreeingHits/nhitsZ;
     } //End HasZed Check
     
     if(rechit->hasPhi()) {
       double countMuonDTHits = 0;
       double countAgreeingHits=0;

       //PREPARE PARAMETERS FOR SEGMENT DETECTOR GEOMETRY
       const DTRecSegment2D* segmPhi;
       segmPhi = dynamic_cast<const DTRecSegment2D*>(rechit->phiSegment());
       nhitsPhi = segmPhi->recHits().size();

       const vector<DTRecHit1D> hits1d = segmPhi->specificRecHits();
       DTChamberId chamberSegIdT((segmPhi->geographicalId()).rawId());

       // Loop over muon recHits
       for(trackingRecHit_iterator hit = track->recHitsBegin(); hit != track->recHitsEnd(); ++hit) {
	 if ( !(*hit)->isValid()) continue; 
         if ( (*hit)->geographicalId().det() != DetId::Muon ) continue;
         if ( (*hit)->geographicalId().subdetId() != MuonSubdetId::DT ) continue;

	 DetId idT = (*hit)->geographicalId();

	 DTChamberId dtDetIdHitT(idT.rawId());

	 DTSuperLayerId dtDetLayerIdHitT(idT.rawId());

	 LocalPoint pointLocal = (*hit)->localPosition(); //Localposition is in DTLayer http://cmslxr.fnal.gov/lxr/source/DataFormats/DTRecHit/interface/DTRecHit1D.h

	 if ((chamberSegIdT==dtDetIdHitT)&&((dtDetLayerIdHitT.superlayer()==1)||(dtDetLayerIdHitT.superlayer()==3))) 
	   countMuonDTHits++;

	 for (vector<DTRecHit1D>::const_iterator hiti=hits1d.begin(); hiti!=hits1d.end(); hiti++) {

	   if ( !hiti->isValid()) continue; 

	   // Pick the one in the same DT Layer as the 1D hit
	   if(!(hiti->geographicalId().rawId()==idT.rawId())) continue; 

	   // and compare the local positions
	   LocalPoint segLocal = hiti->localPosition();
	   //  cout << "     Phi Segment Point = "<<pointLocal<<"    Muon Point = "<<segLocal<<"  Dist:   " 
	   //       << (fabs(pointLocal.x()-segLocal.x()))+(fabs(pointLocal.y()-segLocal.y()))<< endl;

	   if ((fabs(pointLocal.x()-segLocal.x())<1) && 
	       (fabs(pointLocal.y()-segLocal.y())<1))
	     countAgreeingHits++; 
	 } // End Segment Hit Iteration
       } // End Muon Hit Iteration

       matchRatioPhi = countMuonDTHits != 0 ? countAgreeingHits/countMuonDTHits : 0;
       if (nhitsPhi)
	 if (countAgreeingHits/nhitsPhi>matchRatioPhi) matchRatioPhi=countAgreeingHits/nhitsPhi;
     } // End HasPhi Check
     //    DTChamberId chamberSegId2((rechit->geographicalId()).rawId());
        if((matchRatioPhi>0.9 && nhitsPhi)||(matchRatioZ>0.9 && nhitsZ)) {
	 //cout<<"Making a loose match in Chamber "<<chamberSegId2<<endl;
	  totalDTSegs++;
	  if(matchRatioZ>0.9 && nhitsZ) WithZed++;
	}
     }
   }
   return WithZed;
}

double DistToTrigger (const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev){
  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return 9999;
  reco::TrackRef   track = muon->standAloneMuon(); if(track.isNull())return 9999;

  fwlite::Handle< trigger::TriggerEvent > trEvHandle;
  trEvHandle.getByLabel(ev, "hltTriggerSummaryAOD");
  trigger::TriggerEvent trEv = *trEvHandle;

  unsigned int filterIndex = trEv.filterIndex(InputTag("hltL2fL1sMu70Eta2p1L1f0L2Filtered70Q","","HLT"));
  if (filterIndex>=trEv.sizeFilters()) filterIndex = trEv.filterIndex(InputTag("hltL2Mu20L2Filtered20","","HLT"));
  //if (filterIndex>=trEv.sizeFilters()) filterIndex = trEv.filterIndex(InputTag("hltL2Mu60Eta2p1L2Filtered60","","HLT"));
  if (MODE=="ANALYSE_COSMIC") filterIndex = trEv.filterIndex(InputTag("hltL2fL1sMu6NoBPTXL1f0L2Filtered20","","HLT"));
  if (MODE=="ANALYSE_COSMIC" && filterIndex>=trEv.sizeFilters()) filterIndex = trEv.filterIndex(InputTag("hltL2fL1sMu6NoBPTXL1f0L2Filtered10","","HLT"));

  double RMin=999;

  if (filterIndex<trEv.sizeFilters()){

    const trigger::Vids& VIDS(trEv.filterIds(filterIndex));
    const trigger::Keys& KEYS(trEv.filterKeys(filterIndex));
    const int nI(VIDS.size());
    const int nK(KEYS.size());
    assert(nI==nK);
    const int n(std::max(nI,nK));
    const trigger::TriggerObjectCollection& TOC(trEv.getObjects());

    for (int i=0; i!=n; ++i) {
      double eta=TOC[KEYS[i]].eta();
      double phi=TOC[KEYS[i]].phi();
      double dR = deltaR(track->eta(), track->phi(), eta, phi);
      if(dR<RMin)RMin=dR;
    }
  }
  else {
    for ( size_t ia = 0; ia < trEv.sizeFilters(); ++ ia) {
      std::string fullname = trEv.filterTag(ia).encode();
      std::string name;
      size_t p = fullname.find_first_of(':');
      if ( p != std::string::npos) {
	name = fullname.substr(0, p);
      }
      else {
	name = fullname;
      }
      std::cout << "Path name " << name << std::endl;
    }
  }
  return RMin;
}

double RescalePt(double pt) {
  double invpt = 1./pt;
  //invpt=invpt*1.2;
  //invpt=invpt*0.992;
  invpt+=RNG->Gaus(0, 0.19);
  return 1./invpt;
}

