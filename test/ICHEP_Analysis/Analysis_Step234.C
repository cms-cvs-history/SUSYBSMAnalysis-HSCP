
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


namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData;}
namespace susybsm { class HSCParticle; class HSCPIsolation;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}


#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
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

void Analysis_Step2a();
void Analysis_Step2b();
void Analysis_Step3(char* SavePath);
void Analysis_Step4(char* SavePath);

void InitHistos();

double DistToHSCP     (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest);
int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl);
bool   isGoodCandidate(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const fwlite::ChainEvent& ev, const double& PtCut=0, const double& ICut=0,  const double& TOFCut=0, stPlots* st=NULL, const double& PtRescale=1.0, const double& IRescale=1.0);
void DumpCandidateInfo(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, FILE* pFile);
bool PassTrigger      (const fwlite::ChainEvent& ev);
bool hasGoodPtHat     (const fwlite::ChainEvent& ev, const double& PtMax);

void SetWeight(const double& IntegratedLuminosityInPb=-1, const double& CrossSection=0, const double& MCEvents=0);
void SetWeightMC(const double& IntegratedLuminosityInPb, const double& SampleEquivalentLumi, const double& SampleSize, double MaxEvent);

bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold=false);

/////////////////////////// VARIABLE DECLARATION /////////////////////////////

TH1D*  Data_Pt [NSUBSAMPLE];   double Data_Pt_MaxValue;
TH1D*  Data_I  [NSUBSAMPLE];   double Data_I_MaxValue;
TH1D*  Data_TOF[NSUBSAMPLE];   double Data_TOF_MaxValue;

TH1D** Sign_Mass_Syst_PtLow;
TH1D** Sign_Mass_Syst_ILow;
TH1D*  Pred_Mass;
TH1D*  Pred_Mass2;
TH1D*  Pred_Mass3;
TH1D*  Pred_Mass4;

double N_A[NSUBSAMPLE];	double N_Aerr[NSUBSAMPLE];
double N_B[NSUBSAMPLE];	double N_Berr[NSUBSAMPLE];
double N_C[NSUBSAMPLE];	double N_Cerr[NSUBSAMPLE];
double N_D[NSUBSAMPLE];	double N_Derr[NSUBSAMPLE];
double N_E[NSUBSAMPLE]; double N_Eerr[NSUBSAMPLE];
double N_F[NSUBSAMPLE]; double N_Ferr[NSUBSAMPLE];
double N_G[NSUBSAMPLE]; double N_Gerr[NSUBSAMPLE];
double N_H[NSUBSAMPLE]; double N_Herr[NSUBSAMPLE];
double N_P1[NSUBSAMPLE]; double N_P1err[NSUBSAMPLE];
double N_P2[NSUBSAMPLE]; double N_P2err[NSUBSAMPLE];
double N_P3[NSUBSAMPLE]; double N_P3err[NSUBSAMPLE];
double N_P4[NSUBSAMPLE]; double N_P4err[NSUBSAMPLE];
double N_P5[NSUBSAMPLE]; double N_P5err[NSUBSAMPLE];

TH1D*  Pred_P    [NSUBSAMPLE];
TH1D*  Pred_I    [NSUBSAMPLE];
//TH2D*  Data_PI_A [NSUBSAMPLE];
//TH2D*  Data_PI_B [NSUBSAMPLE];
//TH2D*  Data_PI_C [NSUBSAMPLE];
//TH2D*  Data_PI_D [NSUBSAMPLE];
//TH2D*  Pred_PI   [NSUBSAMPLE];
TH1D*  CtrlPt_BckgIs;
TH1D*  CtrlPt_BckgIm;
TH1D*  CtrlPt_BckgTOF;
TH1D*  CtrlPt_SignIs;
TH1D*  CtrlPt_SignIm;
TH1D*  CtrlPt_SignTOF;

TH1D*  CtrlIs_BckgPt;
TH1D*  CtrlIs_BckgTOF;
TH1D*  CtrlIs_SignPt;
TH1D*  CtrlIs_SignTOF;

TH1D*  CtrlTOF_BckgPt;
TH1D*  CtrlTOF_BckgIs;
TH1D*  CtrlTOF_SignPt;
TH1D*  CtrlTOF_SignIs;


TH1D* Pred_Expected_Entries;
TH1D* Pred_Observed_Entries;
//TH1D* Pred_Correlation_A;
//TH1D* Pred_Correlation_B;
//TH1D* Pred_Correlation_C;
//TH1D* Pred_Correlation_D;


stPlots              DataPlots;  
std::vector<stPlots> SignPlots; 
std::vector<stPlots> MCPlots;  
stPlots              MCTrPlots;

std::vector<stSignal> signals;
std::vector<stMC>     MCsample;

/////////////////////////// CODE PARAMETERS /////////////////////////////

std::vector<string> DataFileName;

float Event_Weight = 1;
int MaxEntry = -1;


//void Analysis_Step234(string MODE="COMPILE", double WP_Pt=-1.0, double WP_I=-1, double WP_TOF=-1, int SplitMode_=2, string dEdxSel_="dedxASmi", string dEdxMass_="dedxCNPHarm2", int TypeMode_=0, float MaxEta_=2.5, float MaxPtErr_=0.25)
void Analysis_Step234(string MODE="COMPILE", int TypeMode_=0, int SplitMode_=2, string dEdxSel_="dedxASmi", string dEdxMass_="dedxHarm2", string TOF_Label_="combined", double WP_Pt=-1.0, double WP_I=-1, double WP_TOF=-1, float MinPt_=15, float MaxEta_=2.5, float MaxPtErr_=0.25)
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
   TH1::AddDirectory(kTRUE);

   GetSignalDefinition(signals);
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
   SplitMode = SplitMode_;
   GlobalMaxEta = MaxEta_;
   GlobalMaxPterr = MaxPtErr_;
   GlobalMinPt    = MinPt_;

   if(WP_Pt <=0){   SelectionCutPt  = pow(10,WP_Pt);   DefaultCutPt    =    0;}else{SelectionCutPt  = -1;  DefaultCutPt  = WP_Pt; }
   if(WP_I  <=0){   SelectionCutI   = pow(10,WP_I);    DefaultCutI     =    0;}else{SelectionCutI   = -1;  DefaultCutI   = WP_I;  }
   if(WP_TOF<=0){   SelectionCutTOF = pow(10,WP_TOF);  DefaultCutTOF   = 9999;}else{SelectionCutTOF = -1;  DefaultCutTOF = WP_TOF;}

   if(SplitMode==0){
      GlobalMinNOH = 8;
   }

   if(TypeMode!=2){
      SelectionCutTOF = 1;
      GlobalMinNDOF   = 0; 
      GlobalMinTOF    = 0;
      //GlobalMaxEta    = std::min(GlobalMaxEta,1.2f);
   }



   sprintf(Buffer,"Results/"       );                                          sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%s%s/"         ,Buffer,dEdxS_Label.c_str());                sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%s%s/"         ,Buffer,TOF_Label.c_str());                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sEta%02.0f/"  ,Buffer,10.0*GlobalMaxEta);                  sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sPtMin%02.0f/",Buffer,GlobalMinPt);                        sprintf(Command,"mkdir %s",Buffer); system(Command);
   sprintf(Buffer,"%sType%i/"     ,Buffer,TypeMode);                           sprintf(Command,"mkdir %s",Buffer); system(Command);

   if(MODE=="CUTFINDER"){
      TFile* OutputHisto = new TFile((string(Buffer) + "/CutHistos.root").c_str(),"RECREATE");
      Analysis_Step2a();
      OutputHisto->Write();
      OutputHisto->Close();
      return;
   }else{
      TFile* InputHisto = new TFile((string(Buffer) + "/CutHistos.root").c_str(),"READ");
      for(unsigned int i=0;i<NSUBSAMPLE;i++){
         Data_I  [i] = (TH1D*)GetObjectFromPath(InputHisto, string("CutFinder_I"  ) + GetNameFromIndex(i), true);
         Data_Pt [i] = (TH1D*)GetObjectFromPath(InputHisto, string("CutFinder_Pt" ) + GetNameFromIndex(i), true);
         Data_TOF[i] = (TH1D*)GetObjectFromPath(InputHisto, string("CutFinder_TOF") + GetNameFromIndex(i), true);
      }
      Analysis_Step2b();
      //InputHisto->Close(); //Can't be closed otherwise it releases all the cutfinder histos...
   }




   sprintf(Buffer,"%sSplitMode%i/",Buffer,SplitMode);                          sprintf(Command,"mkdir %s",Buffer); system(Command);

   if(SelectionCutPt >=0){sprintf(Buffer,"%sWPPt%02.0f/" ,Buffer,fabs((10*log10(SelectionCutPt ))));}else{sprintf(Buffer,"%sPt%02.0f/" ,Buffer,DefaultCutPt );} sprintf(Command,"mkdir %s",Buffer); system(Command);
   if(SelectionCutI  >=0){sprintf(Buffer,"%sWPI%02.0f/"  ,Buffer,fabs((10*log10(SelectionCutI  ))));}else{sprintf(Buffer,"%sI%01.2f/"  ,Buffer,DefaultCutI  );} sprintf(Command,"mkdir %s",Buffer); system(Command);
   if(SelectionCutTOF>=0){sprintf(Buffer,"%sWPTOF%02.0f/",Buffer,fabs((10*log10(SelectionCutTOF))));}else{sprintf(Buffer,"%sTOF%01.2f/",Buffer,DefaultCutTOF);} sprintf(Command,"mkdir %s",Buffer); system(Command);

   time_t start = time(NULL);
   TFile* FinalHisto = new TFile((string(Buffer) + "/DumpHistos.root").c_str(),"RECREATE");
   InitHistos();
   Analysis_Step3(Buffer);
   Analysis_Step4(Buffer);
   FinalHisto->Write();
   FinalHisto->Close();
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


bool isGoodCandidate(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const fwlite::ChainEvent& ev, const double& PtCut, const double& ICut, const double& TOFCut, stPlots* st, const double& PtRescale, const double& IRescale)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
   if(TypeMode==2 && hscp.type() != HSCParticleType::globalMuon)return false;
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;
   if(st){st->WN_Total+=Event_Weight;	st->UN_Total++;}

   if(fabs(track->eta())>GlobalMaxEta) return false;

   if(st){st->BS_Hits->Fill(track->found(),Event_Weight);}
   if(track->found()<GlobalMinNOH)return false;
   if(dedxSObj.numberOfMeasurements()<GlobalMinNOM)return false;
   if(st){st->AS_Hits->Fill(track->found(),Event_Weight);}
   if(st){st->WN_Hits  +=Event_Weight;   st->UN_Hits++;}

   double MuonTOF = GlobalMinTOF;
   double NDOF     = 9999;
   if(TypeMode==2 && !hscp.muonRef().isNull()){
      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(ev, "muontiming",TOF_Label.c_str());
      const MuonTimeExtra& tof      = TOFCollH->get(hscp.muonRef().key());
      MuonTOF = tof.inverseBeta();
      NDOF = tof.nDof();
   }

   if(st){st->BS_nDof->Fill(NDOF,Event_Weight);}
   if(NDOF<GlobalMinNDOF)return false;
   if(st){st->AS_nDof->Fill(NDOF,Event_Weight);}
   if(st){st->WN_nDof  +=Event_Weight;   st->UN_nDof++;}

   if(st){st->BS_Qual->Fill(track->qualityMask(),Event_Weight);}
   if(track->qualityMask()<GlobalMinQual )return false;
   if(st){st->AS_Qual->Fill(track->qualityMask(),Event_Weight);}
   if(st){st->WN_Qual  +=Event_Weight;   st->UN_Qual++;}

   if(st){st->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);}
   if(track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(st){st->AS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);}
   if(st){st->WN_Chi2  +=Event_Weight;   st->UN_Chi2++;}

   if(st){st->BS_MPt ->Fill(track->pt(),Event_Weight);}
   if(track->pt()*PtRescale<GlobalMinPt)return false;
   if(st){st->AS_MPt ->Fill(track->pt(),Event_Weight);}
   if(st){st->WN_MPt   +=Event_Weight;   st->UN_MPt ++;}

   if(st){st->BS_MIs->Fill(dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_MIm->Fill(dedxMObj.dEdx(),Event_Weight);}
   if(dedxSObj.dEdx()*IRescale<GlobalMinI)return false;
   if(st){st->AS_MIs->Fill(dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_MIm->Fill(dedxMObj.dEdx(),Event_Weight);}
   if(st){st->WN_MI   +=Event_Weight;   st->UN_MI++;}

   if(st){st->BS_MTOF ->Fill(MuonTOF,Event_Weight);}
   if(MuonTOF<GlobalMinTOF)return false;
   if(st){st->AS_MTOF ->Fill(MuonTOF,Event_Weight);}
   if(st){st->WN_MTOF +=Event_Weight;   st->UN_MTOF++;}

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
   if(st){st->AS_V3D->Fill(v3d,Event_Weight);}
   if(st){st->WN_V3D  +=Event_Weight;   st->UN_V3D++;}

   if(st){st->BS_DZ->Fill(fabs(dz),Event_Weight);}
   if(fabs(dz)>GlobalMaxDZ )return false; 
   if(st){st->AS_DZ->Fill(fabs(dz),Event_Weight);} 
   if(st){st->WN_DZ   +=Event_Weight;   st->UN_DZ++;}

   if(st){st->BS_DXY->Fill(fabs(dxy),Event_Weight);}  
   if(fabs(dxy)>GlobalMaxDXY )return false; 
   if(st){st->AS_DXY->Fill(fabs(dxy),Event_Weight);}
   if(st){st->WN_DXY  +=Event_Weight;   st->UN_DXY++;}

   fwlite::Handle<HSCPIsolationValueMap> IsolationH;
   IsolationH.getByLabel(ev, "HSCPIsolation03");
   if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
   const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();

   HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());

   if(st){st->BS_CIsol ->Fill(hscpIso.Get_TK_Count(),Event_Weight);}
//   if(hscpIso.Get_TK_Count()>0)return false;
   if(st){st->AS_CIsol ->Fill(hscpIso.Get_TK_Count(),Event_Weight);}
   if(st){st->WN_CIsol   +=Event_Weight;   st->UN_CIsol++;}

   if(st){st->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
    if(hscpIso.Get_TK_SumEt()>3)return false;
   if(st){st->AS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
   if(st){st->WN_TIsol   +=Event_Weight;   st->UN_TIsol++;}

   double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
   if(st){st->BS_EIsol ->Fill(EoP,Event_Weight);}
   if(EoP>0.25)return false;
   if(st){st->AS_EIsol ->Fill(EoP,Event_Weight);}
   if(st){st->WN_EIsol   +=Event_Weight;   st->UN_EIsol++;}

   if(st){st->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if((track->ptError()/track->pt())>GlobalMaxPterr)return false;
   if(st){st->AS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if(st){st->WN_Pterr   +=Event_Weight;   st->UN_Pterr ++;}

   if(st){st->BS_P  ->Fill(track->p(),Event_Weight);}
   if(st){st->BS_Pt ->Fill(track->pt(),Event_Weight);}
   if(st){st->BS_Is ->Fill(dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_Im ->Fill(dedxMObj.dEdx(),Event_Weight);}
   if(st){st->BS_TOF->Fill(MuonTOF,Event_Weight);}

   if(st){st->BS_EtaIs->Fill(track->eta(),dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_EtaIm->Fill(track->eta(),dedxMObj.dEdx(),Event_Weight);}
   if(st){st->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);}
   if(st){st->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);}
   if(st){st->BS_PIs  ->Fill(track->p()  ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_PIm  ->Fill(track->p()  ,dedxMObj.dEdx(),Event_Weight);}
   if(st){st->BS_PtIs ->Fill(track->pt() ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_PtIm ->Fill(track->pt() ,dedxMObj.dEdx(),Event_Weight);}
   if(st){st->BS_TOFIs->Fill(MuonTOF     ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->BS_TOFIm->Fill(MuonTOF     ,dedxMObj.dEdx(),Event_Weight);}

   if(track->pt()*PtRescale<PtCut)return false;
   if(st){st->WN_Pt    +=Event_Weight;   st->UN_Pt ++;}
   if(dedxSObj.dEdx()*IRescale<ICut)return false;
   if(st){st->WN_I    +=Event_Weight;   st->UN_I++;}
   if(MuonTOF<TOFCut)return false;
   if(st){st->WN_TOF  +=Event_Weight;   st->UN_TOF++;}

   if(st){st->AS_P  ->Fill(track->p(),Event_Weight);}
   if(st){st->AS_Pt ->Fill(track->pt(),Event_Weight);}
   if(st){st->AS_Is ->Fill(dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_Im ->Fill(dedxMObj.dEdx(),Event_Weight);}
   if(st){st->AS_TOF->Fill(MuonTOF,Event_Weight);}

   if(st){st->AS_EtaIs->Fill(track->eta(),dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_EtaIm->Fill(track->eta(),dedxMObj.dEdx(),Event_Weight);}
   if(st){st->AS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);}
   if(st){st->AS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);}
   if(st){st->AS_PIs  ->Fill(track->p()  ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_PIm  ->Fill(track->p()  ,dedxMObj.dEdx(),Event_Weight);}
   if(st){st->AS_PtIs ->Fill(track->pt() ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_PtIm ->Fill(track->pt() ,dedxMObj.dEdx(),Event_Weight);}
   if(st){st->AS_TOFIs->Fill(MuonTOF     ,dedxSObj.dEdx(),Event_Weight);}
   if(st){st->AS_TOFIm->Fill(MuonTOF     ,dedxMObj.dEdx(),Event_Weight);}

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

   fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
   dEdxSCollH.getByLabel(ev, dEdxS_Label.c_str());
   if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");return;}
   DeDxData dedxSObj  = dEdxSCollH->get(track.key());

   fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
   dEdxMCollH.getByLabel(ev, dEdxM_Label.c_str());
   if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");return;}
   DeDxData dedxMObj  = dEdxMCollH->get(track.key());


   double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
   double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()));
   double Mass = GetMass(PBinned,IBinned);   
   double MassExact = GetMass(track->p(),dedxMObj.dEdx(), true);
   double dz  = track->dz (vertex.position());
   double dxy = track->dxy(vertex.position());

   int HitIndex,EtaIndex;
   GetIndices(dedxSObj.numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
   int CutIndex = GetCutIndex(HitIndex,EtaIndex);

   fprintf(pFile,"\n");
   fprintf(pFile,"---------------------------------------------------------------------------------------------------\n");
   fprintf(pFile,"Candidate Type = %i --> Mass (Binned): %7.2f GeV  Mass (UnBinned): %7.2f\n",hscp.type(),Mass, MassExact);
   fprintf(pFile,"------------------------------------------ EVENT INFO ---------------------------------------------\n");
   fprintf(pFile,"Run=%i Lumi=%i Event=%i BX=%i  Orbit=%i Store=%i\n",ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().event(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().orbitNumber(),ev.eventAuxiliary().storeNumber());
   fprintf(pFile,"------------------------------------------ INNER TRACKER ------------------------------------------\n");
   fprintf(pFile,"Quality = %i Chi2/NDF=%6.2f dz=+%6.2f dxy=%+6.2f charge:%+i\n",track->qualityMask(), track->chi2()/track->ndof(), dz, dxy, track->charge());
   fprintf(pFile,"P=%7.2f  Pt=%7.2f+-%6.2f (Cut=%6.2f) Eta=%+6.2f  Phi=%+6.2f  NOH=%2i\n",track->p(),track->pt(), track->ptError(), CutPt[CutIndex], track->eta(), track->phi(), track->found() );

   fprintf(pFile,"------------------------------------------ DEDX INFO ----------------------------------------------\n");
   fprintf(pFile,"dEdx for selection:%6.2f (Cut=%6.2f) NOM %2i NOS %2i\n",dedxSObj.dEdx(),CutI[CutIndex],dedxSObj.numberOfMeasurements(),dedxSObj.numberOfSaturatedMeasurements());
   fprintf(pFile,"dEdx for mass reco:%6.2f             NOM %2i NOS %2i\n",dedxMObj.dEdx(),dedxMObj.numberOfMeasurements(),dedxMObj.numberOfSaturatedMeasurements());
   if(!muon.isNull()){
      fprintf(pFile,"------------------------------------------ MUON INFO ----------------------------------------------\n");
      fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
      TOFDTCollH.getByLabel(ev, "muontiming","dt");
      if(!TOFDTCollH.isValid()){printf("Invalid TOF DT collection\n");return;}
      MuonTimeExtra tofDT      = TOFDTCollH->get(hscp.muonRef().key());

      fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
      TOFCSCCollH.getByLabel(ev, "muontiming","csc");
      if(!TOFDTCollH.isValid()){printf("Invalid TOF CSC collection\n");return;}
      MuonTimeExtra tofCSC      = TOFCSCCollH->get(hscp.muonRef().key());

      fwlite::Handle<MuonTimeExtraMap> TOFCombCollH;
      TOFCombCollH.getByLabel(ev, "muontiming","combined");
      if(!TOFCombCollH.isValid()){printf("Invalid TOF Combined collection\n");return;}
      MuonTimeExtra tofComb      = TOFCombCollH->get(hscp.muonRef().key());

      fprintf(pFile,"Quality=%i type=%i P=%7.2f  Pt=%7.2f Eta=%+6.2f Phi=%+6.2f #Chambers=%i\n" ,muon->isQualityValid(),muon->type(),muon->p(),muon->pt(),muon->eta(),muon->phi(),muon->numberOfChambers());
      fprintf(pFile,"muonTimeDT      : NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f\n",tofDT  .nDof(),tofDT  .inverseBeta(), tofDT  .inverseBetaErr(), CutTOF[CutIndex], (1.0/tofDT  .inverseBeta()), tofDT  .freeInverseBeta(),tofDT  .freeInverseBetaErr());
      fprintf(pFile,"muonTimeCSC     : NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f\n",tofCSC .nDof(),tofCSC .inverseBeta(), tofCSC .inverseBetaErr(), CutTOF[CutIndex], (1.0/tofCSC .inverseBeta()), tofCSC .freeInverseBeta(),tofCSC .freeInverseBetaErr());
      fprintf(pFile,"muonTimeCombined: NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f\n",tofComb.nDof(),tofComb.inverseBeta(), tofComb.inverseBetaErr(), CutTOF[CutIndex], (1.0/tofComb.inverseBeta()), tofComb.freeInverseBeta(),tofComb.freeInverseBetaErr());
   }
   if(hscp.hasRpcInfo()){
      fprintf(pFile,"------------------------------------------ RPC INFO -----------------------------------------------\n");
      fprintf(pFile,"isCandidate %i Beta=%6.2f\n",hscp.rpc().isCandidate,hscp.rpc().beta);
   }
   if(hscp.hasCaloInfo() && hscp.caloInfoRef()->ecalTime!=-9999){
      fprintf(pFile,"------------------------------------------ CALO INFO ----------------------------------------------\n");
      fprintf(pFile,"HCAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f HO E=%6.2f\n",hscp.caloInfoRef()->hcalCrossedEnergy,hscp.caloInfoRef()->hcal3by3dir, hscp.caloInfoRef()->hcal5by5dir, hscp.caloInfoRef()->hoCrossedEnergy);
      fprintf(pFile,"ECAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f\n"           ,hscp.caloInfoRef()->ecalCrossedEnergy,hscp.caloInfoRef()->ecal3by3dir, hscp.caloInfoRef()->ecal5by5dir);
      fprintf(pFile,"ECAL: time=%6.2f beta=%6.2f trkisodr=%6.2f\n"    ,hscp.caloInfoRef()->ecalTime  ,hscp.caloInfoRef()->ecalBeta   , hscp.caloInfoRef()->trkIsoDr);
   }
   fprintf(pFile,"------------------------------------------ ISOL INFO ----------------------------------------------\n");
   fwlite::Handle<HSCPIsolationValueMap> IsolationH05;
   IsolationH05.getByLabel(ev, "HSCPIsolation05");
   if(!IsolationH05.isValid()){printf("Invalid IsolationH\n");return;}
   const ValueMap<HSCPIsolation>& IsolationMap05 = *IsolationH05.product();

   fwlite::Handle<HSCPIsolationValueMap> IsolationH03;
   IsolationH03.getByLabel(ev, "HSCPIsolation03");
   if(!IsolationH03.isValid()){printf("Invalid IsolationH\n");return;}
   const ValueMap<HSCPIsolation>& IsolationMap03 = *IsolationH03.product();

   fwlite::Handle<HSCPIsolationValueMap> IsolationH01;
   IsolationH01.getByLabel(ev, "HSCPIsolation01");
   if(!IsolationH01.isValid()){printf("Invalid IsolationH\n");return;}
   const ValueMap<HSCPIsolation>& IsolationMap01 = *IsolationH01.product();

   HSCPIsolation hscpIso05 = IsolationMap05.get((size_t)track.key());
   HSCPIsolation hscpIso03 = IsolationMap03.get((size_t)track.key());
   HSCPIsolation hscpIso01 = IsolationMap01.get((size_t)track.key());
   fprintf(pFile,"Isolation05 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",hscpIso05.Get_TK_Count(), hscpIso05.Get_TK_SumEt(), hscpIso05.Get_ECAL_Energy()/track->p(), hscpIso05.Get_HCAL_Energy()/track->p(), (hscpIso05.Get_ECAL_Energy()+hscpIso05.Get_HCAL_Energy())/track->p());
   fprintf(pFile,"Isolation03 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",hscpIso03.Get_TK_Count(), hscpIso03.Get_TK_SumEt(), hscpIso03.Get_ECAL_Energy()/track->p(), hscpIso03.Get_HCAL_Energy()/track->p(), (hscpIso03.Get_ECAL_Energy()+hscpIso03.Get_HCAL_Energy())/track->p());
   fprintf(pFile,"Isolation01 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",hscpIso01.Get_TK_Count(), hscpIso01.Get_TK_SumEt(), hscpIso01.Get_ECAL_Energy()/track->p(), hscpIso01.Get_HCAL_Energy()/track->p(), (hscpIso01.Get_ECAL_Energy()+hscpIso01.Get_HCAL_Energy())/track->p());
   fprintf(pFile,"\n");
}

bool PassTrigger(const fwlite::ChainEvent& ev)
{
      edm::TriggerResultsByName tr = ev.triggerResultsByName("Merge");
      if(!tr.isValid())return false;
      if(!tr.accept(tr.triggerIndex("HscpPathMu")))return true;
      if(!tr.accept(tr.triggerIndex("HscpPathMet")))return true;
      return false;

/*
      //THIS IS ONLY USE FOR SIGNAL SAMPLE (TRIGGER IS ALREADY APPLIED IN BOTH DATA AND MC!)

      edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
      if(!tr.isValid()){
         tr = ev.triggerResultsByName("REDIGI36X");
         if(!tr.isValid())return false;
      }

      fwlite::Handle< trigger::TriggerEvent > trEvHandle;
      trEvHandle.getByLabel(ev,"hltTriggerSummaryAOD");
      trigger::TriggerEvent trEv = *trEvHandle;

      if(tr.accept("HLT_DoubleMu3")                                                    ) return true;
      if(IncreasedTreshold(trEv, InputTag("hltSingleMu11L3Filtered11","","HLT"),15,1)  ) return true;
      if(tr.accept("HLT_MET100")                                                       ) return true;
      //if(IncreasedTreshold(trEv, InputTag("hlt1jet100U"       ,"","HLT"), 140, 1      )) return true;
      //if(IncreasedTreshold(trEv, InputTag("hltDiJetAve70U"    ,"","HLT"), 140, 2, true)) return true;
      //if(tr.accept("HLT_QuadJet25U")                                                   ) return true;
      return false;
*/
}

void Analysis_Step2a()
{
   printf("Step2a: Scanning for Cuts\n");

   char Name   [1024];
   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      sprintf(Name,"CutFinder_I%s",GetNameFromIndex(i).c_str());
      Data_I[i]         = new TH1D(Name,Name, 10000,0,dEdxS_UpLim);
      Data_I[i]->Sumw2();

      sprintf(Name,"CutFinder_Pt%s",GetNameFromIndex(i).c_str());
      Data_Pt[i]       = new TH1D(Name,Name,10000,0,PtHistoUpperBound);
      Data_Pt[i]->Sumw2();

      sprintf(Name,"CutFinder_TOF%s",GetNameFromIndex(i).c_str());
      Data_TOF[i]       = new TH1D(Name,Name,10000,0,50);
      Data_TOF[i]->Sumw2();
   }



   fwlite::ChainEvent tree(DataFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Finding Cuts                 :");
   int TreeStep = tree.size()/50;if(TreeStep==0)TreeStep=1;

   for(Long64_t ientry=0;ientry<tree.size();ientry++){
      tree.to(ientry);
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}
      if(!PassTrigger(tree) )continue;

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(tree,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

      fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
      dEdxSCollH.getByLabel(tree, dEdxS_Label.c_str());
      if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

      fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
      dEdxMCollH.getByLabel(tree, dEdxM_Label.c_str());
      if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(tree, "muontiming",TOF_Label.c_str());
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();

         const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
         const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());

         if(!isGoodCandidate(hscp, dedxSObj, dedxMObj, tree))continue;

         int HitIndex, EtaIndex;
         GetIndices(dedxSObj.numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         
         FillHisto(HitIndex, EtaIndex, Data_Pt, track->pt()                     ,Event_Weight);
         FillHisto(HitIndex, EtaIndex, Data_I , dedxSObj.dEdx() ,Event_Weight);

         if(!hscp.muonRef().isNull()){
            const MuonTimeExtra& tof      = TOFCollH->get(hscp.muonRef().key());
            FillHisto(HitIndex, EtaIndex, Data_TOF , tof.inverseBeta() ,Event_Weight);
         }

      } // end of Track Loop
   }// end of Event Loop
   printf("\n");
}

void Analysis_Step2b()
{
   printf("Step2b: Optimizing Cuts\n");
   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      if(!isSubSampleExist(i))continue;

      if(SelectionCutPt>=0){
         CutPt[i] = CutFromEfficiency(Data_Pt[i],SelectionCutPt);
         if(CutPt[i]<GlobalMinPt)CutPt[i]=GlobalMinPt;
         if(CutPt[i]>Data_Pt[i]->GetXaxis()->GetXmax())CutPt[i]=99999;
         if(CutPt[i]!=99999 && Data_Pt[i]->Integral()*SelectionCutPt<1){printf("Bug --> Set to 99999\n"); CutPt[i]=99999;        }
         printf("SubSample %20s --> %6.2E Entries Eff=%6.2E -->Pt Cut=%5.2f\n",GetNameFromIndex(i).c_str(),Data_Pt[i]->Integral(),SelectionCutPt,CutPt[i]);
      }else{
         CutPt[i]=DefaultCutPt;
         printf("SubSample %20s -->Pt  Cut=%5.2f\n",GetNameFromIndex(i).c_str(),CutPt[i]);  
      }

      if(SelectionCutI>=0){
         CutI[i] = CutFromEfficiency(Data_I[i],SelectionCutI);
         if(CutI[i]<GlobalMinI)CutI[i]=GlobalMinI;
         if(CutI[i]>Data_I[i]->GetXaxis()->GetXmax())CutI[i]=99999;   
         if(CutI[i]!=99999 && Data_I[i]->Integral()*SelectionCutI<1){printf("Bug --> Set to 99999\n"); CutI[i]=99999;        }
         printf("SubSample %20s --> %6.2E Entries Eff=%6.2E -->I Cut=%5.2f\n",GetNameFromIndex(i).c_str(),Data_I[i]->Integral(),SelectionCutI,CutI[i]);
      }else{
         CutI[i]=DefaultCutI;
         printf("SubSample %20s -->I   Cut=%5.2f\n",GetNameFromIndex(i).c_str(),CutI[i]);
      }

      if(SelectionCutTOF>=0){
         CutTOF[i] = CutFromEfficiency(Data_TOF[i],SelectionCutTOF);
         if(CutTOF[i]<GlobalMinTOF)CutTOF[i]=GlobalMinTOF;
         if(CutTOF[i]>Data_TOF[i]->GetXaxis()->GetXmax())CutTOF[i]=99999;
         if(CutTOF[i]!=99999 && Data_TOF[i]->Integral()*SelectionCutTOF<1){printf("Bug --> Set to 99999\n"); CutTOF[i]=99999;        }
         printf("SubSample %20s --> %6.2E Entries Eff=%6.2E -->TOF Cut=%5.2f\n",GetNameFromIndex(i).c_str(),Data_TOF[i]->Integral(),SelectionCutTOF,CutTOF[i]);
      }else{
         CutTOF[i]=DefaultCutTOF;
         printf("SubSample %20s -->TOF Cut=%5.2f\n",GetNameFromIndex(i).c_str(),CutTOF[i]); 
      }
   }
}



void Analysis_Step3(char* SavePath)
{
   printf("Step4: Building Mass Spectrum for B and S\n");

   int TreeStep;
   int HitIndex, EtaIndex;
   FILE* pFile = NULL;
   //FILE* pFileTrg = NULL;

   //////////////////////////////////////////////////     BUILD BACKGROUND MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_D_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");

      //sprintf(Buffer,"%s/Candidate_D_Trigger.txt",SavePath);
      //pFileTrg = fopen(Buffer,"w");
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
      if(!PassTrigger(treeD) )continue;
      DataPlots.WN_TotalTE+=Event_Weight;      DataPlots.UN_TotalTE++;

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
      
      
      bool HSCPTk = false;
      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;

         const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
         const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());

         double MuonTOF = GlobalMinTOF;
         if(TypeMode==2 && !hscp.muonRef().isNull()){
            const MuonTimeExtra& tof      = TOFCollH->get(hscp.muonRef().key());
            MuonTOF = tof.inverseBeta();
         }
 
         GetIndices(dedxSObj.numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
         int CutIndex = GetCutIndex(HitIndex,EtaIndex);
   
         ///////////////////////////////  PREDICTION BEGINS ////////////////////////////////
         if(isGoodCandidate(hscp, dedxSObj, dedxMObj, treeD)){

          //printf("Pt %6.2f|%6.2f  I %6.2f|%6.2f\n", track->pt(),CutPt[CutIndex], dedxSObj.dEdx(), CutI[CutIndex]);

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

            if(track->pt()<25){
               CtrlPt_BckgIs->Fill(dedxSObj.dEdx(), Event_Weight);
               CtrlPt_BckgIm->Fill(dedxMObj.dEdx(), Event_Weight);
               if(MuonTOF!=GlobalMinTOF)CtrlPt_BckgTOF->Fill(MuonTOF, Event_Weight);
            }else{
               CtrlPt_SignIs->Fill(dedxSObj.dEdx(), Event_Weight);
               CtrlPt_SignIm->Fill(dedxMObj.dEdx(), Event_Weight);
               if(MuonTOF!=GlobalMinTOF)CtrlPt_SignTOF->Fill(MuonTOF, Event_Weight);
            }

            if(dedxSObj.dEdx()<0.2){
               CtrlIs_BckgPt->Fill(track->pt(), Event_Weight);
               if(MuonTOF!=GlobalMinTOF)CtrlIs_BckgTOF->Fill(MuonTOF, Event_Weight);
            }else{
               CtrlIs_SignPt->Fill(track->pt(), Event_Weight);
               if(MuonTOF!=GlobalMinTOF)CtrlIs_SignTOF->Fill(MuonTOF, Event_Weight);
            }

            if(MuonTOF!=GlobalMinTOF){
            if(MuonTOF>0.9){
               CtrlTOF_BckgPt->Fill(track->pt()    , Event_Weight);
               CtrlTOF_BckgIs->Fill(dedxSObj.dEdx(), Event_Weight);
            }else{
               CtrlTOF_SignPt->Fill(track->pt()    , Event_Weight);
               CtrlTOF_SignIs->Fill(dedxSObj.dEdx(), Event_Weight);
            }
            }

            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassICut   = (dedxSObj.dEdx()>=CutI[CutIndex]);
            bool PassTOFCut = MuonTOF>=CutTOF[CutIndex];
 
            if(PassTOFCut){
               if(PassPtCut){
                  if(PassICut){ //Region D
                      FillArray(HitIndex, EtaIndex, N_D             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Derr          , Event_Weight*Event_Weight);
                  }else{        //Region C
                      FillArray(HitIndex, EtaIndex, N_C             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Cerr          , Event_Weight*Event_Weight);
                      FillHisto(HitIndex, EtaIndex, Pred_P, track->p(), Event_Weight);
                  }
               }else{
                  if(PassICut){ //Region B
                      FillArray(HitIndex, EtaIndex, N_B             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Berr          , Event_Weight*Event_Weight);
                      FillHisto(HitIndex, EtaIndex, Pred_I, dedxMObj.dEdx(), Event_Weight);
                  }else{        //Region A
                      FillArray(HitIndex, EtaIndex, N_A             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Aerr          , Event_Weight*Event_Weight);
                  }
               }
            }else{
               if(PassPtCut){
                  if(PassICut){ //Region H
                      FillArray(HitIndex, EtaIndex, N_H             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Herr          , Event_Weight*Event_Weight);
                  }else{        //Region G
                      FillArray(HitIndex, EtaIndex, N_G             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Gerr          , Event_Weight*Event_Weight);
                      FillHisto(HitIndex, EtaIndex, Pred_P, track->p(), Event_Weight);
                   }
               }else{
                  if(PassICut){ //Region F
                      FillArray(HitIndex, EtaIndex, N_F             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Ferr          , Event_Weight*Event_Weight);
                      FillHisto(HitIndex, EtaIndex, Pred_I, dedxMObj.dEdx(), Event_Weight);
                  }else{        //Region E
                      FillArray(HitIndex, EtaIndex, N_E             , Event_Weight);
                      FillArray(HitIndex, EtaIndex, N_Eerr          , Event_Weight*Event_Weight);
                  }
               }
            }

         }
         ///////////////////////////////  PREDICTION ENDS   ////////////////////////////////


         //Full Selection
         if(!isGoodCandidate(hscp, dedxSObj, dedxMObj, treeD,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], &DataPlots))continue;



         //DEBUG
         double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
         double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()));
         double Mass = GetMass(PBinned,IBinned);
         if(!isnan(Mass)){
 
            HSCPTk = true;
	    DataPlots.Mass->Fill(Mass,Event_Weight);
            if(SavePath && (Mass>=MinCandidateMass || COUNT_PRINTED_DATA<1000) ){DumpCandidateInfo(hscp, treeD, pFile);   COUNT_PRINTED_DATA++;}
         }
      } // end of Track Loop
      if(HSCPTk){DataPlots.WN_HSCPE+=Event_Weight;  DataPlots.UN_HSCPE++;          }
   }// end of Event Loop
   printf("\n");
   if(pFile){fclose(pFile);pFile=NULL;};
   //if(pFileTrg){fclose(pFileTrg);pFileTrg=NULL;};


      printf("A=%6.2E B=%6.E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E\n",N_A[0], N_B[0], N_C[0], N_D[0], N_E[0], N_F[0], N_G[0], N_H[0] );



   //////////////////////////////////////////////////     BUILD MCTRUTH MASS SPECTRUM

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_M_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");
   }

   for(unsigned int m=0;m<MCsample.size();m++){
      unsigned int COUNT_PRINTED_SIGNAL=0;
      fprintf(pFile,"\nXXXXXXXXXXXXXXXXXXXXXXX %10s XXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n",MCsample[m].Name.c_str());

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


         MCTrPlots .WN_TotalE+=Event_Weight;       MCTrPlots .UN_TotalE++;
         MCPlots[m].WN_TotalE+=Event_Weight;       MCPlots[m].UN_TotalE++;
         if(!PassTrigger(treeM) )continue;
         MCTrPlots .WN_TotalTE+=Event_Weight;      MCTrPlots .UN_TotalTE++;
         MCPlots[m].WN_TotalTE+=Event_Weight;      MCPlots[m].UN_TotalTE++;

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

         
         bool HSCPTk = false;
         for(unsigned int c=0;c<hscpColl.size();c++){
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())continue;

            const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
            const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());

            GetIndices(dedxSObj.numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
            int CutIndex = GetCutIndex(HitIndex,EtaIndex);
                isGoodCandidate(hscp, dedxSObj, dedxMObj, treeM,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], &MCPlots[m]);
            if(!isGoodCandidate(hscp, dedxSObj, dedxMObj, treeM,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], &MCTrPlots))continue;

            double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
            double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()));
            double Mass = GetMass(PBinned,IBinned);
            if(!isnan(Mass)){
               HSCPTk = true;
               MCTrPlots.Mass->Fill(Mass,Event_Weight);
               MCPlots[m].Mass->Fill(Mass,Event_Weight);

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
   }
   if(pFile){fclose(pFile);pFile=NULL;};


   //////////////////////////////////////////////////     BUILD SIGNAL MASS SPECTRUM

   MaxEntry = 25000; //Only look at the first 10K events

   if(SavePath){
      char Buffer[2048];
      sprintf(Buffer,"%s/Candidate_S_Dump.txt",SavePath);
      pFile = fopen(Buffer,"w");

         //sprintf(Buffer,"%s/Candidate_%s_Trigger.txt",SavePath,signals[s].Name.c_str());
         //pFileTrg = fopen(Buffer,"w");
   }

   for(unsigned int s=0;s<signals.size();s++){
      unsigned int COUNT_PRINTED_SIGNAL=0;
      if(SavePath){
         fprintf(pFile,"\n\n\n\n");
         fprintf(pFile,"#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$\n");
         fprintf(pFile,"#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$\n");
         fprintf(pFile,"#$#SIGNAL IS %55s$#$\n",signals[s].Name.c_str());
         fprintf(pFile,"#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$\n");
         fprintf(pFile,"#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$\n\n");
      }

      std::vector<string> SignFileName;
      GetInputFiles(SignFileName, signals[s].Name);

      fwlite::ChainEvent treeS(SignFileName);
      SetWeight(IntegratedLuminosity,signals[s].XSec,(double)std::min((double)treeS.size(),(double)MaxEntry));
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
         int NChargedHSCP=HowManyChargedHSCP(genColl);

         SignPlots[4*s]               .WN_TotalE +=Event_Weight;       SignPlots[4*s]               .UN_TotalE++;
         SignPlots[4*s+NChargedHSCP+1].WN_TotalE +=Event_Weight;       SignPlots[4*s+NChargedHSCP+1].UN_TotalE++;
         if(!PassTrigger(treeS) )continue;
         SignPlots[4*s]               .WN_TotalTE+=Event_Weight;       SignPlots[4*s]               .UN_TotalTE++;
         SignPlots[4*s+NChargedHSCP+1].WN_TotalTE+=Event_Weight;       SignPlots[4*s+NChargedHSCP+1].UN_TotalTE++;

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

            const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
            const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());

            GetIndices(dedxSObj.numberOfMeasurements(), track->eta(),HitIndex,EtaIndex);
            int CutIndex = GetCutIndex(HitIndex,EtaIndex);

            //FOR SYSTEMATIC COMPUTATION (START)
            //A Signal Pt(&P) -->0.95*Pt(&P)
            if(isGoodCandidate(hscp,  dedxSObj, dedxMObj, treeS,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], NULL, 0.95, 1.0)){
               SignPlots[4*s               ].WN_TOF_SYSTA+=Event_Weight; SignPlots[4*s               ].UN_TOF_SYSTA++;
               SignPlots[4*s+NChargedHSCP+1].WN_TOF_SYSTA+=Event_Weight; SignPlots[4*s+NChargedHSCP+1].UN_TOF_SYSTA++;
               double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()*0.95));
               double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()));
               double Mass    = GetMass(PBinned,IBinned);
               if(!isnan(Mass)){
                  HSCPTkSystA    = true;
                  Sign_Mass_Syst_PtLow[4*s               ]->Fill(Mass, Event_Weight);
                  Sign_Mass_Syst_PtLow[4*s+NChargedHSCP+1]->Fill(Mass, Event_Weight);
               }
            }

            //B Signal I -->0.95*I
            if(isGoodCandidate(hscp,  dedxSObj, dedxMObj, treeS,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], NULL, 1.0, 0.95)){
               SignPlots[4*s               ].WN_TOF_SYSTB+=Event_Weight; SignPlots[4*s               ].UN_TOF_SYSTB++;
               SignPlots[4*s+NChargedHSCP+1].WN_TOF_SYSTB+=Event_Weight; SignPlots[4*s+NChargedHSCP+1].UN_TOF_SYSTB++;
               double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
               double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()*0.95));
               double Mass    = GetMass(PBinned,IBinned);
               if(!isnan(Mass)){
                  HSCPTkSystB    = true;
                  Sign_Mass_Syst_ILow[4*s               ]->Fill(Mass, Event_Weight);
                  Sign_Mass_Syst_ILow[4*s+NChargedHSCP+1]->Fill(Mass, Event_Weight);
               }
            }
            //FOR SYSTEMATIC COMPUTATION (END)

                isGoodCandidate(hscp,  dedxSObj, dedxMObj, treeS,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], &SignPlots[4*s+NChargedHSCP+1]);
            if(!isGoodCandidate(hscp,  dedxSObj, dedxMObj, treeS,CutPt[CutIndex], CutI[CutIndex], CutTOF[CutIndex], &SignPlots[4*s               ]))continue;         

            //DEBUG
            double PBinned = Pred_P[0]->GetXaxis()->GetBinCenter(Pred_P[0]->GetXaxis()->FindBin(track->p()));
            double IBinned = Pred_I[0]->GetXaxis()->GetBinCenter(Pred_I[0]->GetXaxis()->FindBin(dedxMObj.dEdx()));
            double Mass = GetMass(PBinned,IBinned);
            if(!isnan(Mass)){

               HSCPTk = true;

               SignPlots[4*s               ].Mass->Fill(Mass,Event_Weight);
               SignPlots[4*s+NChargedHSCP+1].Mass->Fill(Mass,Event_Weight);

               if(SavePath && (Mass>=MinCandidateMass && COUNT_PRINTED_SIGNAL<25)){
                  DumpCandidateInfo(hscp, treeS, pFile);
                  COUNT_PRINTED_SIGNAL++;
               }
            }
         } // end of Track Loop 
         if(HSCPTk     ){SignPlots[4*s               ].WN_HSCPE      +=Event_Weight;  SignPlots[4*s               ].UN_HSCPE++;          }
         if(HSCPTkSystA){SignPlots[4*s               ].WN_HSCPE_SYSTA+=Event_Weight;  SignPlots[4*s               ].UN_HSCPE_SYSTA++;    }
         if(HSCPTkSystB){SignPlots[4*s               ].WN_HSCPE_SYSTB+=Event_Weight;  SignPlots[4*s               ].UN_HSCPE_SYSTB++;    }
         if(HSCPTk     ){SignPlots[4*s+NChargedHSCP+1].WN_HSCPE      +=Event_Weight;  SignPlots[4*s+NChargedHSCP+1].UN_HSCPE++;          }
         if(HSCPTkSystA){SignPlots[4*s+NChargedHSCP+1].WN_HSCPE_SYSTA+=Event_Weight;  SignPlots[4*s+NChargedHSCP+1].UN_HSCPE_SYSTA++;    }
         if(HSCPTkSystB){SignPlots[4*s+NChargedHSCP+1].WN_HSCPE_SYSTB+=Event_Weight;  SignPlots[4*s+NChargedHSCP+1].UN_HSCPE_SYSTB++;    }
       }// end of Event Loop
      printf("\n");
   }// end of signal Type loop
   if(pFile){fclose(pFile);pFile=NULL;};
   //if(pFileTrg){fclose(pFileTrg);pFileTrg=NULL;};

}


void Analysis_Step4(char* SavePath)
{
   //////////////////////////////////////////////////      MAKING THE PREDICTION
   printf("Predicting (Finding Prob)    :");
   int TreeStep = (NSUBSAMPLE)/50;if(TreeStep==0)TreeStep=1;
   int CountStep = 0;

std::cout << "A\n";

   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      if(i%TreeStep==0 && CountStep<=50){printf(".");fflush(stdout);CountStep++;}
      if(!isSubSampleExist(i))continue;
std::cout << "B\n";

      TH1D* tmp_Mass = (TH1D*)Pred_Mass->Clone("PredMassInSubSample");
      tmp_Mass->Reset();
std::cout << "C\n";

      TH1D* tmp_Mass2 = (TH1D*)Pred_Mass2->Clone("PredMassInSubSample2");
      tmp_Mass2->Reset();
std::cout << "D\n";

      TH1D* tmp_Mass3 = (TH1D*)Pred_Mass3->Clone("PredMassInSubSample3");
      tmp_Mass3->Reset();
std::cout << "E\n";

      TH1D* tmp_Mass4 = (TH1D*)Pred_Mass4->Clone("PredMassInSubSample4");
      tmp_Mass4->Reset();
std::cout << "F\n";


      if(N_A[i]>0){
         Pred_Expected_Entries->SetBinContent(i, ((N_C[i]*N_B[i])/N_A[i])  );
         Pred_Expected_Entries->SetBinError  (i, sqrt((pow(N_B[i]/N_A[i],2)*N_Cerr[i]) + (pow(N_C[i]/N_A[i],2)*N_Berr[i]) + (pow((N_B[i]*(N_C[i])/(N_A[i]*N_A[i])),2)*N_Aerr[i])) );
	 Pred_Observed_Entries->SetBinContent(i, N_D[i]  );
         Pred_Observed_Entries->SetBinError  (i, sqrt(N_Derr[i])  );
      }

std::cout << "G\n";


      double IntegralP = Pred_P[i]->Integral(0, Pred_P[i]->GetNbinsX()+1);
      double IntegralI = Pred_I[i]->Integral(0, Pred_I[i]->GetNbinsX()+1);
      if(IntegralP>0)Pred_P[i]->Scale(1.0/IntegralP);
      if(IntegralI>0)Pred_I[i]->Scale(1.0/IntegralI);

std::cout << "H\n";


      double NExpectedBckgEntriesC       = 0;
      double NExpectedBckgEntriesCSqErr  = 0;

      double NExpectedBckgEntriesC2      = 0;
      double NExpectedBckgEntriesC2SqErr = 0;

      double NExpectedBckgEntriesC3      = 0;
      double NExpectedBckgEntriesC3SqErr = 0;

      double NExpectedBckgEntriesC4      = 0;
      double NExpectedBckgEntriesC4SqErr = 0;


std::cout << "I\n";


      N_P1   [i] = N_A[i]>0 ? ((N_C[i]*N_B[i])/N_A[i])                                                                                                             : 0;
      N_P1err[i] = N_A[i]>0 ? sqrt( (pow(N_B[i]/N_A[i],2)*N_Cerr[i]) + (pow(N_C[i]/N_A[i],2)*N_Berr[i]) + (pow((N_B[i]*(N_C[i])/(N_A[i]*N_A[i])),2)*N_Aerr[i])  )  : 0;

      N_P2   [i] = N_E[i]>0 ? N_F[i]*N_C[i]/N_E[i]                                                                                                                 : 0;
      N_P2err[i] = N_E[i]>0 ? sqrt( (pow(N_F[i]/N_E[i],2)*N_Cerr[i]) + (pow(N_C[i]/N_E[i],2)*N_Ferr[i]) + (pow((N_F[i]*(N_C[i])/(N_E[i]*N_E[i])),2)*N_Eerr[i])  )  : 0;

      N_P3   [i] = N_E[i]>0 ? N_G[i]*N_B[i]/N_E[i]                                                                                                                 : 0;
      N_P3err[i] = N_E[i]>0 ? sqrt( (pow(N_B[i]/N_E[i],2)*N_Gerr[i]) + (pow(N_G[i]/N_E[i],2)*N_Berr[i]) + (pow((N_B[i]*(N_G[i])/(N_E[i]*N_E[i])),2)*N_Eerr[i])  )  : 0;

      N_P4   [i] = N_E[i]>0 ? N_A[i]*N_H[i]/N_E[i]                                                                                                                 : 0;
      N_P4err[i] = N_E[i]>0 ? sqrt( (pow(N_H[i]/N_E[i],2)*N_Aerr[i]) + (pow(N_A[i]/N_E[i],2)*N_Herr[i]) + (pow((N_H[i]*(N_A[i])/(N_E[i]*N_E[i])),2)*N_Eerr[i])  )  : 0;

std::cout << "J\n";


      N_P5   [i] = N_E[i]>0 ? (N_A[i]*N_F[i]*N_G[i])/(N_E[i]*N_E[i])                                                                                               : 0;
      N_P5err[i] = N_E[i]>0 ? sqrt( ((pow(N_F[i]*N_G[i],2)*N_Aerr[i] + pow(N_A[i]*N_G[i],2)*N_Ferr[i] + pow(N_A[i]*N_F[i],2)*N_Gerr[i])/pow(N_E[i],4)) + (pow((2*N_A[i]*N_F[i]*N_G[i])/pow(N_E[i],3),2)*N_Eerr[i])) : 0;

std::cout << "K\n";


      if(N_E[i]>0){
         NExpectedBckgEntriesC          = N_P2   [i];
         NExpectedBckgEntriesCSqErr     = N_P2err[i];

         NExpectedBckgEntriesC2         = N_P3   [i];
         NExpectedBckgEntriesC2SqErr    = N_P3err[i];

         NExpectedBckgEntriesC3         = N_P4   [i];
         NExpectedBckgEntriesC3SqErr    = N_P4err[i];

         NExpectedBckgEntriesC4         = N_P5   [i];
         NExpectedBckgEntriesC4SqErr    = N_P5err[i];
      }else if(N_A[i]>0){
         NExpectedBckgEntriesC          = N_P1   [i];
         NExpectedBckgEntriesCSqErr     = N_P1err[i];

         NExpectedBckgEntriesC2         = N_P1   [i];
         NExpectedBckgEntriesC2SqErr    = N_P1err[i];

         NExpectedBckgEntriesC3         = N_P1   [i];
         NExpectedBckgEntriesC3SqErr    = N_P1err[i];

         NExpectedBckgEntriesC4         = N_P1   [i];
         NExpectedBckgEntriesC4SqErr    = N_P1err[i];
      }

std::cout << "L\n";


      printf("A=%6.2E B=%6.E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E\n",N_A[i], N_B[i], N_C[i], N_D[i], N_E[i], N_F[i], N_G[i], N_H[i] );
      printf("B*C/A  = %6.2E +- %6.2E (%6.2E%%)\n", N_P1[i],  N_P1err[i], 100.0*N_P1err[i]/N_P1[i] );
      printf("F*C/E  = %6.2E +- %6.2E (%6.2E%%)\n", N_P2[i],  N_P2err[i], 100.0*N_P2err[i]/N_P2[i] );
      printf("G*B/E  = %6.2E +- %6.2E (%6.2E%%)\n", N_P3[i],  N_P3err[i], 100.0*N_P3err[i]/N_P3[i] );
      printf("A*H/E  = %6.2E +- %6.2E (%6.2E%%)\n", N_P4[i],  N_P4err[i], 100.0*N_P4err[i]/N_P4[i] );
      printf("AFG/EE = %6.2E +- %6.2E (%6.2E%%)\n", N_P5[i],  N_P5err[i], 100.0*N_P5err[i]/N_P5[i] );
      printf("D=%6.2E\n",N_D[i]);

std::cout << "M\n";


      //Loop on Mass Line
      for(int m=0;m<tmp_Mass->GetNbinsX()+1;m++){
         //Find which bins contributes to this particular mass bin
         std::vector<std::pair<int,int> > BinThatGivesThisMass;
         for(int x=1;x<Pred_P[i]->GetNbinsX()+1;x++){
         for(int y=1;y<Pred_I[i]->GetNbinsX()+1;y++){
            double Mass = GetMass( Pred_P[i]->GetXaxis()->GetBinCenter(x) , Pred_I[i]->GetXaxis()->GetBinCenter(y) );
            if(Mass>tmp_Mass->GetXaxis()->GetBinLowEdge(m) && Mass<tmp_Mass->GetXaxis()->GetBinUpEdge(m)){
               BinThatGivesThisMass.push_back(std::make_pair(x,y));
            }
         }}


std::cout << "N\n";


         /// bx1 -->i1; by1-->j1; vx1 --> Ci1 ; vy1 --> Bj1 ; ****MISTAKE - MUST USE MBinContent INSTEAD ***NExpectedBckgEntriesC --> N ********** ; N_A[i] --> A ;         
         double MBinContent =0;
	 double Err_Numer_ijSum=0.;  
	 double Err_Denom_ijSum=0.;
	 double Err_Numer_CorrelSum=0.;
	 double ErrSquared=0.;

         //Loops on the bins that contribute to this mass bin.
         for(unsigned int b1=0;b1<BinThatGivesThisMass.size();b1++){
std::cout << "O\n";

            double bx1 = BinThatGivesThisMass[b1].first;
            double by1 = BinThatGivesThisMass[b1].second;
            double vx1 = Pred_P[i]->GetBinContent(bx1);
            double vy1 = Pred_I[i]->GetBinContent(by1);
            //double ex1 = Pred_P[i]->GetBinError(bx1);
            //double ey1 = Pred_I[i]->GetBinError(by1);
            double vz1 = vx1*vy1;

	    double vxN1=vx1 *IntegralP; 
	    double vyN1=vy1 *IntegralI; 

	    Err_Numer_ijSum += (vxN1*vyN1*(vxN1+vyN1)); 
	    Err_Denom_ijSum += (vxN1*vyN1); // will square at the end of the loop 

            MBinContent  += vz1;

std::cout << "P\n";


            //Compute the errors with a covariance matrix (on the fly) --> Only vertical and horizontal term contributes.
	    ///bx2 -->i2; by2-->j2; vxN2 --> Ci2 ; vyN2 --> Bj2 ;
            for(unsigned int b2=0;b2<BinThatGivesThisMass.size();b2++){
               double bx2 = BinThatGivesThisMass[b2].first;
               double by2 = BinThatGivesThisMass[b2].second;
               double vx2 = Pred_P[i]->GetBinContent(bx2);
               double vy2 = Pred_I[i]->GetBinContent(by2);
	       double vxN2=vx2*IntegralP; 
	       double vyN2=vy2*IntegralI; 

               if(bx1==bx2 && by1==by2){
                  //Correlation with itself!
               }else if(by1==by2){
                  //Vertical term
		  Err_Numer_CorrelSum += vyN2*vxN1*vxN2;
               }else if(bx1==bx2){
                  //Horizontal term
		  Err_Numer_CorrelSum += vxN2*vyN1*vyN2;
               }else{
                  //Diagonal term... do nothing
               }
            }

std::cout << "Q\n";

//            printf("Interval %i --> M = %i --> %f +- %f, %f+-%f\n",i,m,vx1,ex1,vy1,ey1);
//            printf("Interval %i --> M = %i --> %i Bins concerned --> %fEntries -->  %f +- %f\n",i,m,BinThatGivesThisMass.size(),NExpectedBckgEntriesC,MBinContent,NExpectedBckgEntriesC*sqrt(MBinError));
         }

std::cout << "R\n";


	 // squared error on predicted background in considered  mass bin
         if ( (N_A[i]+N_E[i]) != 0 && Err_Denom_ijSum != 0 ) ErrSquared= (1./(N_A[i]+N_E[i]) + (Err_Numer_ijSum + Err_Numer_CorrelSum)/(Err_Denom_ijSum*Err_Denom_ijSum) );  //EXTENDED ABCD
	 //final statistical error
         tmp_Mass ->SetBinContent(m, NExpectedBckgEntriesC  * MBinContent);
         tmp_Mass->SetBinError   (m, MBinContent * sqrt((NExpectedBckgEntriesC *NExpectedBckgEntriesC *ErrSquared) + (NExpectedBckgEntriesCSqErr *NExpectedBckgEntriesCSqErr ) ) );

         tmp_Mass2->SetBinContent(m, NExpectedBckgEntriesC2 * MBinContent);
         tmp_Mass2->SetBinError  (m, MBinContent * sqrt((NExpectedBckgEntriesC2*NExpectedBckgEntriesC2*ErrSquared) + (NExpectedBckgEntriesC2SqErr*NExpectedBckgEntriesC2SqErr) ) );

         tmp_Mass3->SetBinContent(m, NExpectedBckgEntriesC3 * MBinContent);
         tmp_Mass3->SetBinError  (m, MBinContent * sqrt((NExpectedBckgEntriesC3*NExpectedBckgEntriesC3*ErrSquared) + (NExpectedBckgEntriesC3SqErr*NExpectedBckgEntriesC3SqErr) ) );

         tmp_Mass4->SetBinContent(m, NExpectedBckgEntriesC4 * MBinContent);
         tmp_Mass4->SetBinError  (m, MBinContent * sqrt((NExpectedBckgEntriesC4*NExpectedBckgEntriesC4*ErrSquared) + (NExpectedBckgEntriesC4SqErr*NExpectedBckgEntriesC4SqErr) ) );

         BinThatGivesThisMass.clear();
std::cout << "S\n";

      }
      Pred_Mass->Add(tmp_Mass ,1);     
      Pred_Mass2->Add(tmp_Mass2,1);
      Pred_Mass3->Add(tmp_Mass3,1);
      Pred_Mass4->Add(tmp_Mass4,1);
std::cout << "T\n";

      delete tmp_Mass;
      delete tmp_Mass2;
      delete tmp_Mass3;
      delete tmp_Mass4;
std::cout << "U\n";

   }
   printf("\n");



   //////////////////////////////////////////////////     DUMP USEFUL INFORMATION

   char Buffer[2048];
   sprintf(Buffer,"%s/CUT_Dump.txt",SavePath);
   FILE* pFile = fopen(Buffer,"w");
   fprintf(pFile,"MODE          = %i\n",SplitMode);
   fprintf(pFile,"Selection     = %s\n",dEdxS_Label.c_str());
   fprintf(pFile,"Mass          = %s\n",dEdxM_Label.c_str());
   fprintf(pFile,"TOF           = %s\n",TOF_Label.c_str());
   fprintf(pFile,"WP PT         = %4.3E\n",SelectionCutPt);
   fprintf(pFile,"WP I          = %4.3E\n",SelectionCutI);
   fprintf(pFile,"WP TOF        = %4.3E\n",SelectionCutTOF);
   fprintf(pFile,"GlobalMaxEta  = %f\n",GlobalMaxEta);
   fprintf(pFile,"GlobalMaxPterr= %f\n",GlobalMaxPterr);
   fprintf(pFile,"GlobalMinNOH  = %02i\n",GlobalMinNOH);
   fprintf(pFile,"GlobalMinNOM  = %02i\n",GlobalMinNOM);
   fprintf(pFile,"GlobalMinNDOF = %02i\n",GlobalMinNOH);
   fprintf(pFile,"GlobalMaxChi2 = %6.2f\n",GlobalMaxChi2);
   fprintf(pFile,"--------------------\n");

   double CutMin_I   = 999999;   double CutMax_I  = 0;  double CutMean_I   = 0;
   double CutMin_Pt  = 999999;   double CutMax_Pt = 0;  double CutMean_Pt  = 0;
   double CutMin_TOF = 999999;   double CutMax_TOF= 0;  double CutMean_TOF = 0;

   int NCutsI=0;  int NCutsPt=0;   int NCutsTOF=0;
   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      if(!isSubSampleExist(i))continue;

      if(CutI [i]<CutMin_I                                                         )CutMin_I =CutI [i];
      if(CutI [i]>CutMax_I  && CutI [i]<Data_I [0]->GetXaxis()->GetXmax())CutMax_I =CutI [i];
      if(CutI [i]<Data_I [0]->GetXaxis()->GetXmax()){CutMean_I+=CutI [i];NCutsI++;}
      if(CutPt[i]<CutMin_Pt                                                        )CutMin_Pt=CutPt[i];
      if(CutPt[i]>CutMax_Pt && CutPt[i]<Data_Pt[0]->GetXaxis()->GetXmax())CutMax_Pt=CutPt[i];
      if(CutPt[i]<Data_Pt[0]->GetXaxis()->GetXmax()){CutMean_Pt+=CutPt[i];NCutsPt++;}
      if(CutTOF[i]<CutMin_TOF && CutTOF[i]>Data_TOF[0]->GetXaxis()->GetXmin())CutMin_TOF=CutTOF[i];
      if(CutTOF[i]>CutMax_TOF)CutMax_TOF=CutTOF[i];
      if(CutTOF[i]>Data_TOF[0]->GetXaxis()->GetXmin()){CutMean_TOF+=CutTOF[i];NCutsTOF++;}
      fprintf(pFile,"CutIndex=%03i  %20s PtCut=%14.5f   ICut=%14.5f  TOFCut=%4.6f\n",i,GetNameFromIndex(i).c_str(),CutPt[i],CutI[i], CutTOF[i]);
   }
   CutMean_I /= NCutsI;   CutMean_Pt /= NCutsPt;   CutMean_TOF /= NCutsTOF;

   DataPlots.MeanICut = CutMean_I;
   MCTrPlots.MeanICut = CutMean_I;
   for(unsigned int s=0;s<signals.size();s++){
      SignPlots[4*s].MeanICut = CutMean_I;
   }

   DataPlots.MeanPtCut = CutMean_Pt;
   MCTrPlots.MeanPtCut = CutMean_Pt;
   for(unsigned int s=0;s<signals.size();s++){
      SignPlots[4*s].MeanPtCut = CutMean_Pt;
   }

   fprintf(pFile,"--------------------\n");
   fprintf(pFile,"Pt Cut Range =[%10.5f,%10.5f]  Mean = %10.5f\n",CutMin_Pt ,CutMax_Pt,  CutMean_Pt);
   fprintf(pFile,"I  Cut Range =[%10.5f,%10.5f]  Mean = %10.5f\n",CutMin_I  ,CutMax_I,   CutMean_I);
   fprintf(pFile,"TOFCut Range =[%10.5f,%10.5f]  Mean = %10.5f\n",CutMin_TOF,CutMax_TOF, CutMean_TOF);
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
   for(unsigned int s=0;s<SignPlots.size();s++){   
      fprintf(pFile,"##### ##### %10s ##### #####\n",SignPlots[s].Name.c_str());
      stPlots_Dump(SignPlots[s], pFile);
   }

   fprintf(pFile,"\n\n--------------------\n");
   fprintf(pFile,"PREDICTION OF THE MASS DISTRIBUTION\n");
   fprintf(pFile,"--------------------\n");
   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      if(!isSubSampleExist(i))continue;
      fprintf(pFile,"CutIndex=%03i %20s --> N_A=%E, N_B=%E, N_C=%E, N_D=%E  N_E=%E, N_F=%E, N_G=%E, N_H=%E<--> %E+-%E or %E+-%E or %E+-%E\n",i,GetNameFromIndex(i).c_str(),N_A[i],N_B[i],N_C[i],N_D[i],N_E[i],N_F[i],N_G[i],N_H[i], N_P1[i],N_P1err[i],N_P2[i],N_P2err[i],N_P3[i],N_P3err[i] );      
   }
   fprintf(pFile,"--------------------\n");

   fprintf(pFile,"\nIntegral in range [0,2000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D ",GetEventInRange(0,2000,DataPlots.Mass));
   fprintf(pFile,"%15s = %5.3E\n","P ",GetEventInRange(0,2000,Pred_Mass));
   fprintf(pFile,"%15s = %5.3E\n","P2",GetEventInRange(0,2000,Pred_Mass2));
   fprintf(pFile,"%15s = %5.3E\n","P3",GetEventInRange(0,2000,Pred_Mass3));
   fprintf(pFile,"%15s = %5.3E\n","P4",GetEventInRange(0,2000,Pred_Mass4));
   fprintf(pFile,"%15s = %5.3E\n","M ",GetEventInRange(0,2000,MCTrPlots.Mass));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(0,2000,SignPlots[4*s].Mass));
   }
   fprintf(pFile,"\nIntegral in range [75,2000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D ",GetEventInRange(75,2000,DataPlots.Mass));
   fprintf(pFile,"%15s = %5.3E\n","P ",GetEventInRange(75,2000,Pred_Mass));
   fprintf(pFile,"%15s = %5.3E\n","P2",GetEventInRange(75,2000,Pred_Mass2));
   fprintf(pFile,"%15s = %5.3E\n","P3",GetEventInRange(75,2000,Pred_Mass3));
   fprintf(pFile,"%15s = %5.3E\n","P4",GetEventInRange(75,2000,Pred_Mass4));
   fprintf(pFile,"%15s = %5.3E\n","M ",GetEventInRange(75,2000,MCTrPlots.Mass));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(75,2000,SignPlots[4*s].Mass));
   }
   fprintf(pFile,"\nIntegral in range [100,2000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D ",GetEventInRange(100,2000,DataPlots.Mass));
   fprintf(pFile,"%15s = %5.3E\n","P ",GetEventInRange(100,2000,Pred_Mass));
   fprintf(pFile,"%15s = %5.3E\n","P2",GetEventInRange(100,2000,Pred_Mass2));
   fprintf(pFile,"%15s = %5.3E\n","P3",GetEventInRange(100,2000,Pred_Mass3));
   fprintf(pFile,"%15s = %5.3E\n","P4",GetEventInRange(100,2000,Pred_Mass4));
   fprintf(pFile,"%15s = %5.3E\n","M ",GetEventInRange(100,2000,MCTrPlots.Mass));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(100,2000,SignPlots[4*s].Mass));
   }
   fprintf(pFile,"\nIntegral in range [200,2000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D ",GetEventInRange(125,2000,DataPlots.Mass));
   fprintf(pFile,"%15s = %5.3E\n","P ",GetEventInRange(125,2000,Pred_Mass));
   fprintf(pFile,"%15s = %5.3E\n","P2",GetEventInRange(125,2000,Pred_Mass2));
   fprintf(pFile,"%15s = %5.3E\n","P3",GetEventInRange(125,2000,Pred_Mass3));
   fprintf(pFile,"%15s = %5.3E\n","P4",GetEventInRange(125,2000,Pred_Mass4));
   fprintf(pFile,"%15s = %5.3E\n","M ",GetEventInRange(125,2000,MCTrPlots.Mass));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(125,2000,SignPlots[4*s].Mass));
   }
   fprintf(pFile,"\nIntegral in range [300,2000]GeV:\n");
   fprintf(pFile,"%15s = %5.3E\n","D ",GetEventInRange(300,2000,DataPlots.Mass));
   fprintf(pFile,"%15s = %5.3E\n","P ",GetEventInRange(300,2000,Pred_Mass));
   fprintf(pFile,"%15s = %5.3E\n","P2",GetEventInRange(300,2000,Pred_Mass2));
   fprintf(pFile,"%15s = %5.3E\n","P3",GetEventInRange(300,2000,Pred_Mass3));
   fprintf(pFile,"%15s = %5.3E\n","P4",GetEventInRange(300,2000,Pred_Mass4));
   fprintf(pFile,"%15s = %5.3E\n","M ",GetEventInRange(300,2000,MCTrPlots.Mass));
   for(unsigned int s=0;s<signals.size();s++){
   fprintf(pFile,"%15s = %5.3E\n",signals[s].Name.c_str(),GetEventInRange(300,2000,SignPlots[4*s].Mass));
   }
   fprintf(pFile,"--------------------\n");
   fclose(pFile);

   //////////////////////////////////////////////////     CREATE EFFICIENCY FILE

   sprintf(Buffer,"%s/Aeff.tmp",SavePath);
   pFile = fopen(Buffer,"w");
   for(unsigned int s=0;s<signals.size();s++){
      fprintf(pFile,"%15s     Eff=%4.3E (%4.3E)",signals[s].Name.c_str(),  SignPlots[4*s  ].WN_TOF        /(  SignPlots[4*s  ].WN_HSCPE),      SignPlots[4*s  ].UN_TOF       /(  SignPlots[4*s  ].UN_HSCPE        ));
      fprintf(pFile,"SYSTA:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s  ].WN_TOF_SYSTA  /(  SignPlots[4*s  ].WN_HSCPE_SYSTA),SignPlots[4*s  ].UN_TOF_SYSTA /(  SignPlots[4*s  ].UN_HSCPE_SYSTA  ));
      fprintf(pFile,"SYSTB:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s  ].WN_TOF_SYSTB  /(  SignPlots[4*s  ].WN_HSCPE_SYSTB),SignPlots[4*s  ].UN_TOF_SYSTB /(  SignPlots[4*s  ].UN_HSCPE_SYSTB  ));
      fprintf(pFile,"\n");
      fprintf(pFile,"%15s NC0 Eff=%4.3E (%4.3E)",signals[s].Name.c_str(),  SignPlots[4*s+1].WN_TOF        /(  SignPlots[4*s+1].WN_HSCPE),      SignPlots[4*s+1].UN_TOF       /(  SignPlots[4*s+1].UN_HSCPE        ));
      fprintf(pFile,"SYSTA:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+1].WN_TOF_SYSTA  /(  SignPlots[4*s+1].WN_HSCPE_SYSTA),SignPlots[4*s+1].UN_TOF_SYSTA /(  SignPlots[4*s+1].UN_HSCPE_SYSTA  ));
      fprintf(pFile,"SYSTB:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+1].WN_TOF_SYSTB  /(  SignPlots[4*s+1].WN_HSCPE_SYSTB),SignPlots[4*s+1].UN_TOF_SYSTB /(  SignPlots[4*s+1].UN_HSCPE_SYSTB  ));
      fprintf(pFile,"\n");
      fprintf(pFile,"%15s NC1 Eff=%4.3E (%4.3E)",signals[s].Name.c_str(),  SignPlots[4*s+2].WN_TOF        /(  SignPlots[4*s+2].WN_HSCPE),      SignPlots[4*s+2].UN_TOF       /(  SignPlots[4*s+2].UN_HSCPE        ));
      fprintf(pFile,"SYSTA:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+2].WN_TOF_SYSTA  /(  SignPlots[4*s+2].WN_HSCPE_SYSTA),SignPlots[4*s+2].UN_TOF_SYSTA /(  SignPlots[4*s+2].UN_HSCPE_SYSTA  ));
      fprintf(pFile,"SYSTB:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+2].WN_TOF_SYSTB  /(  SignPlots[4*s+2].WN_HSCPE_SYSTB),SignPlots[4*s+2].UN_TOF_SYSTB /(  SignPlots[4*s+2].UN_HSCPE_SYSTB  ));
      fprintf(pFile,"\n");
      fprintf(pFile,"%15s NC2 Eff=%4.3E (%4.3E)",signals[s].Name.c_str(),  SignPlots[4*s+3].WN_TOF        /(  SignPlots[4*s+3].WN_HSCPE),      SignPlots[4*s+3].UN_TOF       /(  SignPlots[4*s+3].UN_HSCPE        ));
      fprintf(pFile,"SYSTA:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+3].WN_TOF_SYSTA  /(  SignPlots[4*s+3].WN_HSCPE_SYSTA),SignPlots[4*s+3].UN_TOF_SYSTA /(  SignPlots[4*s+3].UN_HSCPE_SYSTA  ));
      fprintf(pFile,"SYSTB:   Eff=%4.3E (%4.3E)",                          SignPlots[4*s+3].WN_TOF_SYSTB  /(  SignPlots[4*s+3].WN_HSCPE_SYSTB),SignPlots[4*s+3].UN_TOF_SYSTB /(  SignPlots[4*s+3].UN_HSCPE_SYSTB  ));
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
   for(int NC=0;NC<4;NC++){
      stPlots tmp;
      if(NC==0){
         stPlots_Init(tmp,signals[s].Name);
      }else{
         char buffer[256];sprintf(buffer,"_NC%i",NC-1);
         stPlots_Init(tmp,signals[s].Name + buffer,true);
      }
      SignPlots.push_back(tmp);
   }}

   Pred_Expected_Entries = new TH1D("Pred_Expected_Entries","Pred_Expected_Entries",NSUBSAMPLE,0,NSUBSAMPLE);
   Pred_Observed_Entries = new TH1D("Pred_Observed_Entries","Pred_Observed_Entries",NSUBSAMPLE,0,NSUBSAMPLE);

   //Pred_Correlation_A = new TH1D("Pred_Correlation_A","Pred_Correlation_A",NSUBSAMPLE,0,NSUBSAMPLE);
   //Pred_Correlation_B = new TH1D("Pred_Correlation_B","Pred_Correlation_B",NSUBSAMPLE,0,NSUBSAMPLE);
   //Pred_Correlation_C = new TH1D("Pred_Correlation_C","Pred_Correlation_C",NSUBSAMPLE,0,NSUBSAMPLE);
   //Pred_Correlation_D = new TH1D("Pred_Correlation_D","Pred_Correlation_D",NSUBSAMPLE,0,NSUBSAMPLE);

   CtrlPt_BckgIs   = new TH1D("CtrlPt_BckgIs" ,"CtrlPt_BckgIs" ,200,0,dEdxS_UpLim);  CtrlPt_BckgIs ->Sumw2();
   CtrlPt_BckgIm   = new TH1D("CtrlPt_BckgIm" ,"CtrlPt_BckgIm" ,200,0,dEdxM_UpLim);  CtrlPt_BckgIm ->Sumw2();
   CtrlPt_BckgTOF  = new TH1D("CtrlPt_BckgTOF","CtrlPt_BckgTOF",200,0,20);           CtrlPt_BckgTOF->Sumw2();
   CtrlPt_SignIs   = new TH1D("CtrlPt_SignIs" ,"CtrlPt_SignIs" ,200,0,dEdxS_UpLim);  CtrlPt_SignIs ->Sumw2();
   CtrlPt_SignIm   = new TH1D("CtrlPt_SignIm" ,"CtrlPt_SignIm" ,200,0,dEdxM_UpLim);  CtrlPt_SignIm ->Sumw2();
   CtrlPt_SignTOF  = new TH1D("CtrlPt_SignTOF","CtrlPt_SignTOF",200,0,20);           CtrlPt_SignTOF->Sumw2();

   CtrlIs_BckgPt   = new TH1D("CtrlIs_BckgPt" ,"CtrlIs_BckgPt" ,200,0,1500);         CtrlIs_BckgPt ->Sumw2();
   CtrlIs_BckgTOF  = new TH1D("CtrlIs_BckgTOF","CtrlIs_BckgTOF",200,0,20);           CtrlIs_BckgTOF->Sumw2();
   CtrlIs_SignPt   = new TH1D("CtrlIs_SignPt" ,"CtrlIs_SignPt" ,200,0,1500);         CtrlIs_SignPt ->Sumw2();
   CtrlIs_SignTOF  = new TH1D("CtrlIs_SignTOF","CtrlIs_SignTOF",200,0,20);           CtrlIs_SignTOF->Sumw2();

   CtrlTOF_BckgPt  = new TH1D("CtrlTOF_BckgPt","CtrlTOF_BckgPt",200,0,1500);         CtrlTOF_BckgPt->Sumw2();
   CtrlTOF_BckgIs  = new TH1D("CtrlTOF_BckgIs","CtrlTOF_BckgIs",200,0,dEdxS_UpLim);  CtrlTOF_BckgIs->Sumw2();
   CtrlTOF_SignPt  = new TH1D("CtrlTOF_SignPt","CtrlTOF_SignPt",200,0,1500);         CtrlTOF_SignPt->Sumw2();
   CtrlTOF_SignIs  = new TH1D("CtrlTOF_SignIs","CtrlTOF_SignIs",200,0,dEdxS_UpLim);  CtrlTOF_SignIs->Sumw2();


   char Name   [1024];
   sprintf(Name,"Pred_Mass");
   Pred_Mass = new TH1D(Name,Name,400,0,MassHistoUpperBound);
   Pred_Mass->Sumw2();

   sprintf(Name,"Pred_Mass2");
   Pred_Mass2 = new TH1D(Name,Name,400,0,MassHistoUpperBound);
   Pred_Mass2->Sumw2();

   sprintf(Name,"Pred_Mass3");
   Pred_Mass3 = new TH1D(Name,Name,400,0,MassHistoUpperBound);
   Pred_Mass3->Sumw2();

   sprintf(Name,"Pred_Mass4");
   Pred_Mass4 = new TH1D(Name,Name,400,0,MassHistoUpperBound);
   Pred_Mass4->Sumw2();



   for(unsigned int i=0;i<NSUBSAMPLE;i++){
      sprintf(Name,"Pred_I%s",GetNameFromIndex(i).c_str());
      Pred_I[i]  = new TH1D(Name,Name,200,0,dEdxM_UpLim);
      Pred_I[i]->Sumw2();

      sprintf(Name,"Pred_P%s",GetNameFromIndex(i).c_str());
      Pred_P[i]  = new TH1D(Name,Name,200,0,PtHistoUpperBound);
      Pred_P[i]->Sumw2();

      //sprintf(Name,"Pred_PI%s",GetNameFromIndex(i).c_str());
      //Pred_PI[i] = new TH2D(Name,Name,400,0,PtHistoUpperBound, 400, 0, dEdxM_UpLim);
      //Pred_PI[i]->Sumw2();

      //sprintf(Name,"Data_PI_A%s",GetNameFromIndex(i).c_str());
      //Data_PI_A[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxM_UpLim);
      //Data_PI_A[i]->Sumw2();

      //sprintf(Name,"Data_PI_B%s",GetNameFromIndex(i).c_str());
      //Data_PI_B[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxM_UpLim);
      //Data_PI_B[i]->Sumw2();

      //sprintf(Name,"Data_PI_C%s",GetNameFromIndex(i).c_str());
      //Data_PI_C[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxM_UpLim);
      //Data_PI_C[i]->Sumw2();

      //sprintf(Name,"Data_PI_D%s",GetNameFromIndex(i).c_str());
      //Data_PI_D[i] = new TH2D(Name,Name,100,0,PtHistoUpperBound, 100, 0, dEdxM_UpLim);
      //Data_PI_D[i]->Sumw2();

      N_A[i] = 0;    N_Aerr[i] = 0;
      N_B[i] = 0;    N_Berr[i] = 0;
      N_C[i] = 0;    N_Cerr[i] = 0;
      N_D[i] = 0;    N_Derr[i] = 0;
      N_E[i] = 0;    N_Eerr[i] = 0;
      N_F[i] = 0;    N_Ferr[i] = 0;
      N_G[i] = 0;    N_Gerr[i] = 0;
      N_H[i] = 0;    N_Herr[i] = 0;
   }

   Sign_Mass_Syst_PtLow = new TH1D*[4*signals.size()];
   Sign_Mass_Syst_ILow  = new TH1D*[4*signals.size()];
   for(unsigned int s=0;s<signals.size();s++){
      for(unsigned int n=0;n<4;n++){
         if(n==0){sprintf(Name,"%s_Mass_Syst_PtLow", signals[s].Name.c_str());
         }else{   sprintf(Name,"%s_NC%i_Mass_Syst_PtLow", signals[s].Name.c_str(),n-1); }
         Sign_Mass_Syst_PtLow[4*s+n] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
         Sign_Mass_Syst_PtLow[4*s+n]->Sumw2();

         if(n==0){sprintf(Name,"%s_Mass_Syst_ILow", signals[s].Name.c_str());
         }else{   sprintf(Name,"%s_NC%i_Mass_Syst_ILow", signals[s].Name.c_str(),n-1); }
         Sign_Mass_Syst_ILow[4*s+n] = new TH1D(Name,Name,200,0,MassHistoUpperBound);
         Sign_Mass_Syst_ILow[4*s+n]->Sumw2();
      }
   }
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

void SetWeight(const double& IntegratedLuminosityInPb, const double& CrossSection, const double& MCEvents){
   if(IntegratedLuminosityInPb>0){
      Event_Weight = (CrossSection * IntegratedLuminosityInPb) / MCEvents;
   }else{
      Event_Weight=1;
   }
}

void SetWeightMC(const double& IntegratedLuminosityInPb, const double& SampleEquivalentLumi, const double& SampleSize, double MaxEvent){
   if(MaxEvent<0)MaxEvent=SampleSize;
   printf("SetWeight MC: IntLumi = %6.2E  SampleLumi = %6.2E --> EventWeight = %6.2E\n",IntegratedLuminosityInPb,SampleEquivalentLumi, IntegratedLuminosityInPb/SampleEquivalentLumi);
   printf("Sample NEvent = %6.2E   SampleEventUsed = %6.2E --> Weight Rescale = %6.2E\n",SampleSize, MaxEvent, SampleSize/MaxEvent);
   Event_Weight = (IntegratedLuminosityInPb/SampleEquivalentLumi) * (SampleSize/MaxEvent);
   printf("FinalWeight = %6.2f\n",Event_Weight);
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



bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold)
{
   unsigned int filterIndex = trEv.filterIndex(InputPath);
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
         }
         if(NObjectAboveThresholdObserved>=NObjectAboveThreshold)return true;

      }else{
         std::vector<double> ObjPt;

         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            ObjPt.push_back(TO.pt());
         }
         if((int)(ObjPt.size())<NObjectAboveThreshold)return false;
         std::sort(ObjPt.begin(), ObjPt.end());

         double Average = 0;
         for(int i=0; i<NObjectAboveThreshold;i++){
            Average+= ObjPt[ObjPt.size()-1-i];
         }Average/=NObjectAboveThreshold;
         if(Average>NewThreshold)return true;
      }
   }
   return false;
}

