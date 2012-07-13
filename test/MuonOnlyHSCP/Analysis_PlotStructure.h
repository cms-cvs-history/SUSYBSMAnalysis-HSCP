#include "Analysis_Samples.h"
#include "TGraphErrors.h"

struct stPlots {
   bool        SelPlot;
   std::string Name;
   TDirectory* Directory;
   TTree*      Tree;
   unsigned int NCuts;
   unsigned int Tree_Run;
   unsigned int Tree_Event;
   unsigned int Tree_Hscp;
   float        Tree_Pt;
   float        Tree_TOF;

  TH2F* Eta_RegionA;
  TH2F* Eta_RegionB;
  TH2F* Eta_RegionC;
  TH2F* Eta_RegionD;

   TH1F* TotalE;
   TH1F* TotalEPU; 
   TH1F* TotalTE;
   TH1F* Total;

  TH1F* Reconstructed;
  TH1F* TriggerMatch;
  TH1F* tofFound;
  TH1F* Preselected;
  TH1F* Preselected_DT;
  TH1F* Preselected_CSC;
  TH1F* Eta;
  TH1F* NVTrack;
  TH1F* MinPt;
  TH1F* Stations;
  TH1F* nDof;
  TH1F* tofError;
  TH1F* Dxy;
  TH1F* Dz;
  TH1F* Pt;	
  TH1F* TOF;
  TH1F* TOFExclusive;
  TH1F* SegSep;
  TH1F* FailDz;
  TH1F* FailDz_DT;
  TH1F* FailDz_CSC;

  TH1F* HSCPE;
  TH1F* HSCPE_SystT;
  TH1F* HSCPE_SystP;

   TH1F* Beta_Gen;
   TH1F* Beta_GenCharged;
   TH1F* Beta_Triggered;
   TH1F* Beta_Matched;
   TH1F* Beta_PreselectedC;
   TH2F* Beta_SelectedP;
   TH2F* Beta_SelectedT;

  TH1F* Pt_Gen;

  TH1F* DistToGen;
  TH1F* DistTrigger;

  TH1F*  BS_Dz_FailDxy;
  TH1F*  BS_Dz_PassDxy;
  TH1F*  BS_Eta_FailDxy;
  TH1F*  BS_Eta_PassDxy;
  TH1F*  BS_Eta_FailDz;
  TH1F*  BS_Eta_PassDz;

  TH2F*  BS_Dxy_Dz;
  TH1F*  BS_Dxy;
  TH1F*  BS_Dz;
  TH1F*  BS_Dz_CSC;
  TH1F*  BS_Dz_DT;
  TH1F*  BS_Dz_FailPhi;
  TH1F*  BS_Dz_PassPhi;
  TH1F*  BS_Dz_NoZed;
  TH1F*  BS_Dz_Def;
  TH1F*  BS_V3D;
   TH1F*  BS_Chi2;
   TH1F*  BS_TNOH;
   TH1F*  BS_TNOH_Barrel;
   TH1F*  BS_TNOH_Endcap;
   TH1F*  BS_Eta;
  TH1F*  BS_Phi;
   TH1F*  BS_nDof;
   TH1F*  BS_Pterr;
   TH1F*  BS_PterrSq;
  TH1F*  BS_Pt_PassDz;
  TH1F*  BS_Pt_FailDz;
  TH1F*  BS_Pt_PassDz_DT;
  TH1F*  BS_Pt_FailDz_DT;
  TH1F*  BS_Pt_PassDz_CSC;
  TH1F*  BS_Pt_FailDz_CSC;
  TH1F*  BS_Pt_PassDxy;
  TH1F*  BS_Pt_FailDxy;
  TH1F*  BS_GlobalPt_FailDxy;
  TH1F*  BS_GlobalPt_PassDxy;
  TH1F*  BS_Pt_FailPhi;
  TH1F*  BS_Pt_PassPhi;
  TH1F*  BS_TOF_FailDz;
  TH1F*  BS_TOF_PassDz;
  TH1F*  BS_TOF_FailDz_DT;
  TH1F*  BS_TOF_PassDz_DT;
  TH1F*  BS_TOF_FailDz_CSC;
  TH1F*  BS_TOF_PassDz_CSC;

  TH1F*  BS_Time_FailDz;
  TH1F*  BS_Time_PassDz;
  TH1F*  BS_Time_FailDz_DT;
  TH1F*  BS_Time_PassDz_DT;
  TH1F*  BS_Time_FailDz_CSC;
  TH1F*  BS_Time_PassDz_CSC;

  TH1F*  BS_TOF_FailDxy;
  TH1F*  BS_TOF_PassDxy;
  TH1F*  BS_TOF_FailPhi;
  TH1F*  BS_TOF_PassPhi;

  TH1F*  BS_BetaDiff;
  TH1F*  BS_Pt_All;
  TH1F*  BS_TOF_All;
  TH1F*  BS_VertexTime;
  TH1F*  BS_SegSep;
  TH1F*  BS_SegMinEtaSep_CSC;
  TH1F*  BS_SegMinEtaSep_DT;
  TH1F*  BS_SegMinPhiSep;
  TH1F*  BS_SegMinEtaSep;
  TH1F*  BS_SegMinEtaSep_FailDz;
  TH1F*  BS_SegMinEtaSep_PassDz;
  TH2F*  BS_SegEta_MinEtaSep;
  TH1F*  BS_IsTracker;
  TH1F*  BS_IsUpdated;
  TH1F*  BS_MatchedStations;
  TH1F*  BS_InvBetaErr;
  TH1F*  BS_PV;
  TH1F*  BS_dR_NVTrack;
  TH1F*  BS_Pt_Global;

  TH1F*  BS_P; 	 
  TH1F*  BS_Pt;	 
  TH1F*  BS_Pt_Bar;
  TH1F*  BS_Pt_For;
  TH1F*  BS_QoverPt;
  TH1F*  BS_TOF; 
  TH1F*  BS_TOF_Bar;
  TH1F*  BS_TOF_For;
  TH1F*  BS_TOF_DT;
  TH1F*  BS_TOF_CSC;

   TH2F*  BS_EtaP;
   TH2F*  BS_EtaPt;
   TH2F*  BS_EtaTOF;
  TH2F*  BS_EtaTime;
  TH2F*  BS_EtaDz;
  TH2F*  BS_DzTime;
  TH2F*  BS_DzTime_DT;
  TH2F*  BS_DzTime_CSC;
  TH2F*  BS_PhiTime;
  TH2F*  BS_DzPt;
  TH2F*  BS_PtTOF;  

  TH1D* H_A;
  TH1D* H_B;
  TH1D* H_C;
  TH1D* H_D;

  TH1D* H_A_Cen;
  TH1D* H_B_Cen;
  TH1D* H_C_Cen;
  TH1D* H_D_Cen;

  TH1D* H_A_For;
  TH1D* H_B_For;
  TH1D* H_C_For;
  TH1D* H_D_For;

  TH1D* H_A_Low;
  TH1D* H_B_Low;
  TH1D* H_C_Low;
  TH1D* H_D_Low;

  TH1D* H_A_Cen_Low;
  TH1D* H_B_Cen_Low;
  TH1D* H_C_Cen_Low;
  TH1D* H_D_Cen_Low;

  TH1D* H_A_For_Low;
  TH1D* H_B_For_Low;
  TH1D* H_C_For_Low;
  TH1D* H_D_For_Low;

  TH1D* H_A_Syst[DzRegions];
  TH1D* H_B_Syst[DzRegions];
  TH1D* H_C_Syst[DzRegions];
  TH1D* H_D_Syst[DzRegions];

  TH1D* H_A_Cen_Syst[DzRegions];
  TH1D* H_B_Cen_Syst[DzRegions];
  TH1D* H_C_Cen_Syst[DzRegions];
  TH1D* H_D_Cen_Syst[DzRegions];

  TH1D* H_A_For_Syst[DzRegions];
  TH1D* H_B_For_Syst[DzRegions];
  TH1D* H_C_For_Syst[DzRegions];
  TH1D* H_D_For_Syst[DzRegions];

  TH1D* H_DzCounts;
  TH1D* H_DzCounts_DT;
  TH1D* H_DzCounts_CSC;

  TH1D*  CtrlPt_S1_TOF;
  TH1D*  CtrlPt_S2_TOF;
  TH1D*  CtrlPt_S3_TOF;
  TH1D*  CtrlPt_S4_TOF;

  TH1D*  CtrlCen_Pt_S1_TOF;
  TH1D*  CtrlCen_Pt_S2_TOF;
  TH1D*  CtrlCen_Pt_S3_TOF;
  TH1D*  CtrlCen_Pt_S4_TOF;

  TH1D*  CtrlFor_Pt_S1_TOF;
  TH1D*  CtrlFor_Pt_S2_TOF;
  TH1D*  CtrlFor_Pt_S3_TOF;
  TH1D*  CtrlFor_Pt_S4_TOF;

  TH1D*  CtrlTOF_S1_Pt;
  TH1D*  CtrlTOF_S2_Pt;
  TH1D*  CtrlTOF_S3_Pt;
  TH1D*  CtrlTOF_S4_Pt;

  TH1D*  CtrlCen_TOF_S1_Pt;
  TH1D*  CtrlCen_TOF_S2_Pt;
  TH1D*  CtrlCen_TOF_S3_Pt;
  TH1D*  CtrlCen_TOF_S4_Pt;

  TH1D*  CtrlFor_TOF_S1_Pt;
  TH1D*  CtrlFor_TOF_S2_Pt;
  TH1D*  CtrlFor_TOF_S3_Pt;
  TH1D*  CtrlFor_TOF_S4_Pt;

  //TH1D*  HCuts_Pt;
  //TH1D*  HCuts_I;
  //TH1D*  HCuts_TOF;
};

void stPlots_Init(TFile* HistoFile, stPlots& st, std::string BaseName, unsigned int NCuts, bool SkipTree=true)
{
  st.SelPlot = true;
   st.Name = BaseName;
   st.NCuts = NCuts;

   std::string Name;
   Name = BaseName;               st.Directory = HistoFile->mkdir(Name.c_str(), Name.c_str()); 
   st.Directory->cd();

   Name = "TotalE";   st.TotalE  = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "TotalEPU"; st.TotalEPU= new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "TotalTE";  st.TotalTE = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Total";    st.Total   = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Reconstructed";      st.Reconstructed     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "TriggerMatch";      st.TriggerMatch     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "tofFound";      st.tofFound     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Preselected";      st.Preselected     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Preselected_DT";      st.Preselected_DT     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Preselected_CSC";      st.Preselected_CSC     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Eta";      st.Eta     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "NVTrack";      st.NVTrack     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "MinPt";      st.MinPt     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Stations";      st.Stations     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "tofError";      st.tofError     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Dxy";      st.Dxy     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "Dz";      st.Dz     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "nDof";     st.nDof    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Pt";       st.Pt      = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     
   Name = "TOF";      st.TOF     = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     
   Name = "TOFExclusive"; st.TOFExclusive     = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "SegSep";      st.SegSep     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz";      st.FailDz     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz_DT";      st.FailDz_DT     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz_CSC";      st.FailDz_CSC     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);

   Name = "HSCPE";    st.HSCPE   = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystT";    st.HSCPE_SystT   = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystP";    st.HSCPE_SystP   = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);

   Name = "Eta_RegionA";    st.Eta_RegionA    = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, 60, -2.5, 2.5);   st.Eta_RegionA    ->Sumw2();
   Name = "Eta_RegionB";    st.Eta_RegionB    = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, 60, -2.5, 2.5);   st.Eta_RegionB    ->Sumw2();
   Name = "Eta_RegionC";    st.Eta_RegionC    = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, 60, -2.5, 2.5);   st.Eta_RegionC    ->Sumw2();
   Name = "Eta_RegionD";    st.Eta_RegionD    = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, 60, -2.5, 2.5);   st.Eta_RegionD    ->Sumw2();

   Name = "Beta_Gen"         ; st.Beta_Gen         = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Gen         ->Sumw2();
   Name = "Beta_Triggered"   ; st.Beta_Triggered   = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Triggered   ->Sumw2();
   Name = "Beta_Matched"     ; st.Beta_Matched     = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Matched     ->Sumw2();
   Name = "Beta_PreselectedC"; st.Beta_PreselectedC= new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_PreselectedC->Sumw2();
   Name = "Beta_SelectedP"   ; st.Beta_SelectedP   = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts,   20, 0,  1);  st.Beta_SelectedP   ->Sumw2();
   Name = "Beta_SelectedT"   ; st.Beta_SelectedT   = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts,   20, 0,  1);  st.Beta_SelectedT   ->Sumw2();

   Name = "Pt_Gen"         ; st.Pt_Gen         = new TH1F(Name.c_str(), Name.c_str(),                 100, 0,  1000);  st.Pt_Gen         ->Sumw2();

   Name = "DistToGen"  ; st.DistToGen  = new TH1F(Name.c_str(), Name.c_str(),  200,  0,  2);                st.DistToGen->Sumw2();
   Name = "DistTrigger"; st.DistTrigger  = new TH1F(Name.c_str(), Name.c_str(),  1000,  0,  10);                st.DistTrigger->Sumw2();

   Name = "BS_Dz_FailDxy"  ; st.BS_Dz_FailDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -100,  100);                st.BS_Dz_FailDxy->Sumw2();
   Name = "BS_Dz_PassDxy"  ; st.BS_Dz_PassDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -100,  100);                st.BS_Dz_PassDxy->Sumw2();
   Name = "BS_Eta_FailDxy"  ; st.BS_Eta_FailDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_FailDxy->Sumw2();
   Name = "BS_Eta_PassDxy"  ; st.BS_Eta_PassDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_PassDxy->Sumw2();
   Name = "BS_Eta_FailDz"  ; st.BS_Eta_FailDz   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_FailDz->Sumw2();
   Name = "BS_Eta_PassDz"  ; st.BS_Eta_PassDz   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_PassDz->Sumw2();
   Name = "BS_V3D"; st.BS_V3D   = new TH1F(Name.c_str(), Name.c_str(),  IPHistoUpperBound,  0,  IPHistoUpperBound);                st.BS_V3D->Sumw2();
   Name = "BS_Dxy";st.BS_Dxy   = new TH1F(Name.c_str(), Name.c_str(),  2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dxy->Sumw2();
   Name = "BS_Dz"; st.BS_Dz   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz->Sumw2();
   Name = "BS_Dz_CSC"; st.BS_Dz_CSC = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_CSC->Sumw2();
   Name = "BS_Dz_DT"; st.BS_Dz_DT=new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_DT->Sumw2();
   Name = "BS_Dz_FailPhi"; st.BS_Dz_FailPhi   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_FailPhi->Sumw2();
   Name = "BS_Dz_PassPhi"; st.BS_Dz_PassPhi   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_PassPhi->Sumw2();
   Name = "BS_Dz_NoZed"; st.BS_Dz_NoZed   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_NoZed->Sumw2();
   Name = "BS_Dz_Def"; st.BS_Dz_Def   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_Def->Sumw2();
   Name = "BS_Dxy_Dz";st.BS_Dxy_Dz = new TH2F(Name.c_str(), Name.c_str(),2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound,2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound); st.BS_Dxy_Dz->Sumw2();
   Name = "BS_Chi2" ; st.BS_Chi2  = new TH1F(Name.c_str(), Name.c_str(),  200,  0,  50);                st.BS_Chi2->Sumw2();
   Name = "BS_TNOH" ; st.BS_TNOH  = new TH1F(Name.c_str(), Name.c_str(),  61,  -0.5,  60.5);                st.BS_TNOH->Sumw2();
   Name = "BS_TNOH_Barrel" ; st.BS_TNOH_Barrel  = new TH1F(Name.c_str(), Name.c_str(),  61,  -0.5,  60.5);                st.BS_TNOH_Barrel->Sumw2();
   Name = "BS_TNOH_Endcap" ; st.BS_TNOH_Endcap  = new TH1F(Name.c_str(), Name.c_str(),  61,  -0.5,  60.5);                st.BS_TNOH_Endcap->Sumw2();
   Name = "BS_Eta" ; st.BS_Eta  = new TH1F(Name.c_str(), Name.c_str(),  50,  -2.6,  2.6);                st.BS_Eta->Sumw2();
   Name = "BS_Phi" ; st.BS_Phi  = new TH1F(Name.c_str(), Name.c_str(),  50,  -3.14,  3.14);              st.BS_Phi->Sumw2();
   Name = "BS_nDof" ; st.BS_nDof  = new TH1F(Name.c_str(), Name.c_str(),  20,  0, 40);  st.BS_nDof->Sumw2();
   Name = "BS_PtErr"; st.BS_Pterr = new TH1F(Name.c_str(), Name.c_str(),  80,  0,  5);                st.BS_Pterr->Sumw2();
   Name = "BS_PtErrSq"; st.BS_PterrSq = new TH1F(Name.c_str(), Name.c_str(),  80,  0,  0.04);                st.BS_PterrSq->Sumw2();
   Name = "BS_Pt_All"; st.BS_Pt_All = new TH1F(Name.c_str(), Name.c_str(),  700,  0,  PtHistoUpperBound);                st.BS_Pt_All->Sumw2();
   Name = "BS_TOF_All"; st.BS_TOF_All = new TH1F(Name.c_str(), Name.c_str(),  80,  -4,  6);                st.BS_TOF_All->Sumw2();

   Name = "BS_Pt_FailDz"; st.BS_Pt_FailDz = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_FailDz->Sumw2();
   Name = "BS_Pt_PassDz"; st.BS_Pt_PassDz = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_PassDz->Sumw2();
   Name = "BS_Pt_FailDz_DT"; st.BS_Pt_FailDz_DT = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_FailDz_DT->Sumw2();
   Name = "BS_Pt_PassDz_DT"; st.BS_Pt_PassDz_DT = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_PassDz_DT->Sumw2();
   Name = "BS_Pt_FailDz_CSC"; st.BS_Pt_FailDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_FailDz_CSC->Sumw2();
   Name = "BS_Pt_PassDz_CSC"; st.BS_Pt_PassDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_PassDz_CSC->Sumw2();
   Name = "BS_Pt_FailDxy"; st.BS_Pt_FailDxy = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_FailDxy->Sumw2();
   Name = "BS_Pt_PassDxy"; st.BS_Pt_PassDxy = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_PassDxy->Sumw2();
   Name = "BS_GlobalPt_FailDxy"; st.BS_GlobalPt_FailDxy = new TH1F(Name.c_str(), Name.c_str(),  200, 0, 1000); st.BS_GlobalPt_FailDxy->Sumw2();
   Name = "BS_GlobalPt_PassDxy"; st.BS_GlobalPt_PassDxy = new TH1F(Name.c_str(), Name.c_str(),  200, 0, 1000); st.BS_GlobalPt_PassDxy->Sumw2();
   Name = "BS_Pt_FailPhi"; st.BS_Pt_FailPhi = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_FailPhi->Sumw2();
   Name = "BS_Pt_PassPhi"; st.BS_Pt_PassPhi = new TH1F(Name.c_str(), Name.c_str(),  600, 0, PtHistoUpperBound); st.BS_Pt_PassPhi->Sumw2();
   Name = "BS_TOF_FailDz"; st.BS_TOF_FailDz = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_FailDz->Sumw2();
   Name = "BS_TOF_PassDz"; st.BS_TOF_PassDz = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_PassDz->Sumw2();
   Name = "BS_TOF_FailDz_DT"; st.BS_TOF_FailDz_DT = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_FailDz_DT->Sumw2();
   Name = "BS_TOF_PassDz_DT"; st.BS_TOF_PassDz_DT = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_PassDz_DT->Sumw2();
   Name = "BS_TOF_FailDz_CSC"; st.BS_TOF_FailDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_FailDz_CSC->Sumw2();
   Name = "BS_TOF_PassDz_CSC"; st.BS_TOF_PassDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_PassDz_CSC->Sumw2();
   Name = "BS_Time_FailDz"; st.BS_Time_FailDz = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_FailDz->Sumw2();
   Name = "BS_Time_PassDz"; st.BS_Time_PassDz = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_PassDz->Sumw2();
   Name = "BS_Time_FailDz_DT"; st.BS_Time_FailDz_DT = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_FailDz_DT->Sumw2();
   Name = "BS_Time_PassDz_DT"; st.BS_Time_PassDz_DT = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_PassDz_DT->Sumw2();
   Name = "BS_Time_FailDz_CSC"; st.BS_Time_FailDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_FailDz_CSC->Sumw2();
   Name = "BS_Time_PassDz_CSC"; st.BS_Time_PassDz_CSC = new TH1F(Name.c_str(), Name.c_str(),  50, -100, 100); st.BS_Time_PassDz_CSC->Sumw2();
   Name = "BS_TOF_FailDxy"; st.BS_TOF_FailDxy = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_FailDxy->Sumw2();
   Name = "BS_TOF_PassDxy"; st.BS_TOF_PassDxy = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_PassDxy->Sumw2();
   Name = "BS_TOF_FailPhi"; st.BS_TOF_FailPhi = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_FailPhi->Sumw2();
   Name = "BS_TOF_PassPhi"; st.BS_TOF_PassPhi = new TH1F(Name.c_str(), Name.c_str(),  50, -4, 6); st.BS_TOF_PassPhi->Sumw2();

   Name = "BS_BetaDiff"; st.BS_BetaDiff = new TH1F(Name.c_str(), Name.c_str(),  60,  -0.5,  0.5); st.BS_BetaDiff->Sumw2();
   Name = "BS_Pt_Global"   ; st.BS_Pt_Global    = new TH1F(Name.c_str(), Name.c_str(),                   200, 0, 1500); st.BS_Pt_Global->Sumw2();

   Name = "BS_P"    ; st.BS_P     = new TH1F(Name.c_str(), Name.c_str(),                   160, 0, PtHistoUpperBound); st.BS_P->Sumw2();
   Name = "BS_Pt"   ; st.BS_Pt    = new TH1F(Name.c_str(), Name.c_str(),                   160, 70, PtHistoUpperBound); st.BS_Pt->Sumw2();
   Name = "BS_Pt_Bar"   ; st.BS_Pt_Bar    = new TH1F(Name.c_str(), Name.c_str(),                   160, 70, PtHistoUpperBound); st.BS_Pt_Bar->Sumw2();
   Name = "BS_Pt_For"   ; st.BS_Pt_For    = new TH1F(Name.c_str(), Name.c_str(),                   160, 70, PtHistoUpperBound); st.BS_Pt_For->Sumw2();
   Name = "BS_QoverPt"   ; st.BS_QoverPt    = new TH1F(Name.c_str(), Name.c_str(),                   300, -5, 3); st.BS_QoverPt->Sumw2();
   Name = "BS_TOF"  ; st.BS_TOF   = new TH1F(Name.c_str(), Name.c_str(),                   150, -3, 5);                 st.BS_TOF->Sumw2();
   Name = "BS_TOF_Bar"  ; st.BS_TOF_Bar   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_Bar->Sumw2();
   Name = "BS_TOF_For"  ; st.BS_TOF_For   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_For->Sumw2();
   Name = "BS_TOF_DT"  ; st.BS_TOF_DT   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_DT->Sumw2();
   Name = "BS_TOF_CSC"  ; st.BS_TOF_CSC   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_CSC->Sumw2();

   Name = "BS_VertexTime"  ; st.BS_VertexTime= new TH1F(Name.c_str(), Name.c_str(),                   300, -150, 150); st.BS_VertexTime->Sumw2();
   Name = "BS_SegSep"  ; st.BS_SegSep= new TH1F(Name.c_str(), Name.c_str(),                   1000, 0, 10.0); st.BS_SegSep->Sumw2();
   Name = "BS_SegMinEtaSep_DT"  ; st.BS_SegMinEtaSep_DT= new TH1F(Name.c_str(), Name.c_str(),                   1000, -5, 5); st.BS_SegMinEtaSep_DT->Sumw2();
   Name = "BS_SegMinEtaSep_CSC"  ; st.BS_SegMinEtaSep_CSC= new TH1F(Name.c_str(), Name.c_str(),                   1000, -5, 5.0); st.BS_SegMinEtaSep_CSC->Sumw2();
   Name = "BS_SegMinEtaSep"  ; st.BS_SegMinEtaSep= new TH1F(Name.c_str(), Name.c_str(),                   1000, -5., 5.); st.BS_SegMinEtaSep->Sumw2();
   Name = "BS_SegMinPhiSep"  ; st.BS_SegMinPhiSep= new TH1F(Name.c_str(), Name.c_str(),                   1000, -3.3, 3.3); st.BS_SegMinPhiSep->Sumw2();
   Name = "BS_SegMinEtaSep_FailDz"  ; st.BS_SegMinEtaSep_FailDz= new TH1F(Name.c_str(), Name.c_str(),                   1000, -5., 5.); st.BS_SegMinEtaSep_FailDz->Sumw2();
   Name = "BS_SegMinEtaSep_PassDz"  ; st.BS_SegMinEtaSep_PassDz= new TH1F(Name.c_str(), Name.c_str(),                   1000, -5., 5.); st.BS_SegMinEtaSep_PassDz->Sumw2();
   Name = "BS_SegEta_MinEtaSep"  ; st.BS_SegEta_MinEtaSep= new TH2F(Name.c_str(), Name.c_str(), 1000, -2.1, 2.1, 1000, -5., 5.); st.BS_SegEta_MinEtaSep->Sumw2();
   Name = "BS_IsTracker"  ; st.BS_IsTracker = new TH1F(Name.c_str(), Name.c_str(),                   2, -0.5, 1.5); st.BS_IsTracker->Sumw2();
   Name = "BS_IsUpdated"  ; st.BS_IsUpdated = new TH1F(Name.c_str(), Name.c_str(),                   2, -0.5, 1.5); st.BS_IsUpdated->Sumw2();
   Name = "BS_MatchedStations"  ; st.BS_MatchedStations= new TH1F(Name.c_str(), Name.c_str(),                   8, -0.5, 7.5); st.BS_MatchedStations->Sumw2();
   Name = "BS_InvBetaErr"  ; st.BS_InvBetaErr= new TH1F(Name.c_str(), Name.c_str(),                   100, 0, 1); st.BS_InvBetaErr->Sumw2();
   Name = "BS_PV"  ; st.BS_PV = new TH1F(Name.c_str(), Name.c_str(),                   101, -0.5, 100.5); st.BS_PV->Sumw2();
   Name = "BS_dR_NVTrack"  ; st.BS_dR_NVTrack = new TH1F(Name.c_str(), Name.c_str(),                   100, 0, 1); st.BS_dR_NVTrack->Sumw2();

   Name = "BS_EtaP" ; st.BS_EtaP  = new TH2F(Name.c_str(), Name.c_str(),                   100,-3, 3, 150, 0, PtHistoUpperBound);
   Name = "BS_EtaPt"; st.BS_EtaPt = new TH2F(Name.c_str(), Name.c_str(),                   100,-3, 3, 150, 0, PtHistoUpperBound);
   Name = "BS_EtaTOF" ; st.BS_EtaTOF  = new TH2F(Name.c_str(), Name.c_str(),                   100,-3, 3, 100, -2, 5);
   Name = "BS_EtaTime" ; st.BS_EtaTime  = new TH2F(Name.c_str(), Name.c_str(),                   100,-3, 3, 100, -200, 200);
   Name = "BS_EtaDz" ; st.BS_EtaDz  = new TH2F(Name.c_str(), Name.c_str(),  100,-3, 3, 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound);
   Name = "BS_PhiTime" ; st.BS_PhiTime  = new TH2F(Name.c_str(), Name.c_str(),                   50,-3.14, 3.14, 100, -200, 200);
   Name = "BS_DzTime" ; st.BS_DzTime  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzTime_DT" ; st.BS_DzTime_DT  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzTime_CSC" ; st.BS_DzTime_CSC  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzPt" ; st.BS_DzPt  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 150, 0, PtHistoUpperBound);
   Name = "BS_PtTOF" ; st.BS_PtTOF= new TH2F(Name.c_str(), Name.c_str(),                   150, 0, PtHistoUpperBound, 100, 0.5, 1.5);

   Name = "H_A"; st.H_A = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A->Sumw2();
   Name = "H_B"; st.H_B = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B->Sumw2();
   Name = "H_C"; st.H_C = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C->Sumw2();
   Name = "H_D"; st.H_D = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D->Sumw2();

   Name = "H_A_Cen"; st.H_A_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen->Sumw2();
   Name = "H_B_Cen"; st.H_B_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen->Sumw2();
   Name = "H_C_Cen"; st.H_C_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen->Sumw2();
   Name = "H_D_Cen"; st.H_D_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen->Sumw2();

   Name = "H_A_For"; st.H_A_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For->Sumw2();
   Name = "H_B_For"; st.H_B_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For->Sumw2();
   Name = "H_C_For"; st.H_C_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For->Sumw2();
   Name = "H_D_For"; st.H_D_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For->Sumw2();

   Name = "H_A_Low"; st.H_A_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Low->Sumw2();
   Name = "H_B_Low"; st.H_B_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Low->Sumw2();
   Name = "H_C_Low"; st.H_C_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Low->Sumw2();
   Name = "H_D_Low"; st.H_D_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Low->Sumw2();

   Name = "H_A_Cen_Low"; st.H_A_Cen_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen_Low->Sumw2();
   Name = "H_B_Cen_Low"; st.H_B_Cen_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen_Low->Sumw2();
   Name = "H_C_Cen_Low"; st.H_C_Cen_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen_Low->Sumw2();
   Name = "H_D_Cen_Low"; st.H_D_Cen_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen_Low->Sumw2();

   Name = "H_A_For_Low"; st.H_A_For_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For_Low->Sumw2();
   Name = "H_B_For_Low"; st.H_B_For_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For_Low->Sumw2();
   Name = "H_C_For_Low"; st.H_C_For_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For_Low->Sumw2();
   Name = "H_D_For_Low"; st.H_D_For_Low = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For_Low->Sumw2();

   for(int i=0; i<DzRegions; i++) {
     Name = "H_A_Syst_"+RegionNames[i]; st.H_A_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Syst[i]->Sumw2();
     Name = "H_B_Syst_"+RegionNames[i]; st.H_B_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Syst[i]->Sumw2();
     Name = "H_C_Syst_"+RegionNames[i]; st.H_C_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Syst[i]->Sumw2();
     Name = "H_D_Syst_"+RegionNames[i]; st.H_D_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Syst[i]->Sumw2();

     Name = "H_A_Cen_Syst_"+RegionNames[i]; st.H_A_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen_Syst[i]->Sumw2();
     Name = "H_B_Cen_Syst_"+RegionNames[i]; st.H_B_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen_Syst[i]->Sumw2();
     Name = "H_C_Cen_Syst_"+RegionNames[i]; st.H_C_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen_Syst[i]->Sumw2();
     Name = "H_D_Cen_Syst_"+RegionNames[i]; st.H_D_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen_Syst[i]->Sumw2();

     Name = "H_A_For_Syst_"+RegionNames[i]; st.H_A_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For_Syst[i]->Sumw2();
     Name = "H_B_For_Syst_"+RegionNames[i]; st.H_B_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For_Syst[i]->Sumw2();
     Name = "H_C_For_Syst_"+RegionNames[i]; st.H_C_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For_Syst[i]->Sumw2();
     Name = "H_D_For_Syst_"+RegionNames[i]; st.H_D_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For_Syst[i]->Sumw2();
   }
   Name = "H_DzCounts"; st.H_DzCounts = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts->Sumw2();
   Name = "H_DzCounts_DT"; st.H_DzCounts_DT = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts_DT->Sumw2();
   Name = "H_DzCounts_CSC"; st.H_DzCounts_CSC = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts_CSC->Sumw2();

   Name = "CtrlPt_S1_TOF"; st.CtrlPt_S1_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlPt_S1_TOF->Sumw2();
   Name = "CtrlPt_S2_TOF"; st.CtrlPt_S2_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlPt_S2_TOF->Sumw2();
   Name = "CtrlPt_S3_TOF"; st.CtrlPt_S3_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlPt_S3_TOF->Sumw2();
   Name = "CtrlPt_S4_TOF"; st.CtrlPt_S4_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlPt_S4_TOF->Sumw2();

   Name = "CtrlCen_Pt_S1_TOF"; st.CtrlCen_Pt_S1_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlCen_Pt_S1_TOF->Sumw2();
   Name = "CtrlCen_Pt_S2_TOF"; st.CtrlCen_Pt_S2_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlCen_Pt_S2_TOF->Sumw2();
   Name = "CtrlCen_Pt_S3_TOF"; st.CtrlCen_Pt_S3_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlCen_Pt_S3_TOF->Sumw2();
   Name = "CtrlCen_Pt_S4_TOF"; st.CtrlCen_Pt_S4_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlCen_Pt_S4_TOF->Sumw2();

   Name = "CtrlFor_Pt_S1_TOF"; st.CtrlFor_Pt_S1_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlFor_Pt_S1_TOF->Sumw2();
   Name = "CtrlFor_Pt_S2_TOF"; st.CtrlFor_Pt_S2_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlFor_Pt_S2_TOF->Sumw2();
   Name = "CtrlFor_Pt_S3_TOF"; st.CtrlFor_Pt_S3_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlFor_Pt_S3_TOF->Sumw2();
   Name = "CtrlFor_Pt_S4_TOF"; st.CtrlFor_Pt_S4_TOF = new TH1D(Name.c_str(), Name.c_str(),200,-2,7); st.CtrlFor_Pt_S4_TOF->Sumw2();

   Name = "CtrlTOF_S1_Pt"; st.CtrlTOF_S1_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlTOF_S1_Pt->Sumw2();
   Name = "CtrlTOF_S2_Pt"; st.CtrlTOF_S2_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlTOF_S2_Pt->Sumw2();
   Name = "CtrlTOF_S3_Pt"; st.CtrlTOF_S3_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlTOF_S3_Pt->Sumw2();
   Name = "CtrlTOF_S4_Pt"; st.CtrlTOF_S4_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlTOF_S4_Pt->Sumw2();

   Name = "CtrlCen_TOF_S1_Pt"; st.CtrlCen_TOF_S1_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlCen_TOF_S1_Pt->Sumw2();
   Name = "CtrlCen_TOF_S2_Pt"; st.CtrlCen_TOF_S2_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlCen_TOF_S2_Pt->Sumw2();
   Name = "CtrlCen_TOF_S3_Pt"; st.CtrlCen_TOF_S3_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlCen_TOF_S3_Pt->Sumw2();
   Name = "CtrlCen_TOF_S4_Pt"; st.CtrlCen_TOF_S4_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlCen_TOF_S4_Pt->Sumw2();

   Name = "CtrlFor_TOF_S1_Pt"; st.CtrlFor_TOF_S1_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlFor_TOF_S1_Pt->Sumw2();
   Name = "CtrlFor_TOF_S2_Pt"; st.CtrlFor_TOF_S2_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlFor_TOF_S2_Pt->Sumw2();
   Name = "CtrlFor_TOF_S3_Pt"; st.CtrlFor_TOF_S3_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlFor_TOF_S3_Pt->Sumw2();
   Name = "CtrlFor_TOF_S4_Pt"; st.CtrlFor_TOF_S4_Pt = new TH1D(Name.c_str(), Name.c_str(),200,80,1580); st.CtrlFor_TOF_S4_Pt->Sumw2();

   if(!SkipTree) {
   st.Tree = new TTree("HscpCandidates", "HscpCandidates");
   st.Tree->SetDirectory(0);
   st.Tree->Branch("Run"     ,&st.Tree_Run       ,"Run/i");
   st.Tree->Branch("Event"   ,&st.Tree_Event     ,"Event/i");
   st.Tree->Branch("Hscp"    ,&st.Tree_Hscp      ,"Hscp/i");
   st.Tree->Branch("Pt"      ,&st.Tree_Pt        ,"Pt/F");
   st.Tree->Branch("TOF"     ,&st.Tree_TOF       ,"TOF/F");
   }
   HistoFile->cd();
}


void stPlots_InitFromFile(TFile* HistoFile, stPlots& st, std::string BaseName, TFile* InputFile)
{
   st.Name = BaseName;
   std::string Name;
   Name = BaseName;

   st.Directory = new TDirectory((Name+"ReadFromFile").c_str(), (Name+"ReadFromFile").c_str());
   st.Directory->cd();
   TDirectory::AddDirectory(kTRUE);
   TH1::AddDirectory(kTRUE);

   st.TotalE            = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TotalE");
   st.TotalEPU          = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TotalEPU");
   st.TotalTE           = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TotalTE");
   st.Total             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Total");
   st.Reconstructed               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Reconstructed");
   st.TriggerMatch               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TriggerMatch");
   st.tofFound               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/tofFound");
   st.Preselected               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Preselected");
   st.Preselected_DT               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Preselected_DT");
   st.Preselected_CSC              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Preselected_CSC");
   st.Eta               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Eta");
   st.NVTrack               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/NVTrack");
   st.MinPt               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MinPt");
   st.Stations               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Stations");
   st.Dxy               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Dxy");
   st.Dz               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Dz");
   st.nDof              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/nDof");
   st.tofError               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/tofError");
   st.Pt                = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Pt");
   st.TOF               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TOF");
   st.TOFExclusive               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TOFExclusive");
   st.HSCPE             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE");
   st.SegSep               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/SegSep");
   st.FailDz               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/FailDz");
   st.FailDz_DT               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/FailDz_DT");
   st.FailDz_CSC               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/FailDz_CSC");

   st.HSCPE             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE");
   st.HSCPE_SystT             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystT");
   st.HSCPE_SystP             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystP");

   st.Eta_RegionA           = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Eta_RegionA");
   st.Eta_RegionB           = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Eta_RegionB");
   st.Eta_RegionC           = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Eta_RegionC");
   st.Eta_RegionD           = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Eta_RegionD");

   st.Beta_Gen          = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Gen");
   st.Beta_Triggered    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Triggered");
   st.Beta_Matched      = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Matched");
   st.Beta_PreselectedC = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_PreselectedC");
   st.Beta_SelectedP    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_SelectedP");
   st.Beta_SelectedT    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_SelectedT");

   st.Pt_Gen          = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Pt_Gen");

   st.DistToGen    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/DistToGen");
   st.DistTrigger    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/DistTrigger");

   st.BS_Dz_FailDxy    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_FailDxy");
   st.BS_Dz_PassDxy    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_PassDxy");
   st.BS_Eta_FailDxy    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta_FailDxy");
   st.BS_Eta_PassDxy    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta_PassDxy");
   st.BS_Eta_FailDz    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta_FailDz");
   st.BS_Eta_PassDz    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta_PassDz");
   st.BS_V3D    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_V3D");
   st.BS_Dxy    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dxy");
   st.BS_Dz    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz");
   st.BS_Dz_CSC    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_CSC");
   st.BS_Dz_DT    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_DT");
   st.BS_Dz_FailPhi    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_FailPhi");
   st.BS_Dz_PassPhi    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_PassPhi");
   st.BS_Dz_NoZed    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dz_NoZed");
   st.BS_Dxy_Dz    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Dxy_Dz");
   st.BS_Chi2   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Chi2");
   st.BS_TNOH   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOH");
   st.BS_TNOH_Barrel   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOH_Barrel");
   st.BS_TNOH_Endcap   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOH_Endcap");
   st.BS_Eta    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta");
   st.BS_Phi    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Phi");
   st.BS_nDof   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_nDof");
   st.BS_Pterr  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtErr");
   st.BS_PterrSq= (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtErrSq");
   st.BS_Pt_All  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_All");
   st.BS_Pt_Global  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_Global");
   st.BS_TOF_All  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_All");

   st.BS_Pt_FailDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_FailDz");
   st.BS_Pt_FailDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_FailDz_DT");
   st.BS_Pt_FailDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_FailDz_CSC");
   st.BS_Pt_PassDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_PassDz");
   st.BS_Pt_PassDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_PassDz_DT");
   st.BS_Pt_PassDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_PassDz_CSC");
   st.BS_Pt_FailDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_FailDxy");
   st.BS_Pt_PassDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_PassDxy");
   st.BS_GlobalPt_FailDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_GlobalPt_FailDxy");
   st.BS_GlobalPt_PassDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_GlobalPt_PassDxy");
   st.BS_Pt_FailPhi  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_FailPhi");
   st.BS_Pt_PassPhi  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_PassPhi");
   st.BS_TOF_FailDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_FailDz");
   st.BS_TOF_PassDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_PassDz");
   st.BS_TOF_FailDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_FailDz_CSC");
   st.BS_TOF_PassDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_PassDz_CSC");
   st.BS_TOF_FailDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_FailDz_DT");
   st.BS_TOF_PassDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_PassDz_DT");
   st.BS_Time_FailDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_FailDz");
   st.BS_Time_PassDz  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_PassDz");
   st.BS_Time_FailDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_FailDz_CSC");
   st.BS_Time_PassDz_CSC  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_PassDz_CSC");
   st.BS_Time_FailDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_FailDz_DT");
   st.BS_Time_PassDz_DT  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Time_PassDz_DT");
   st.BS_TOF_FailDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_FailDxy");
   st.BS_TOF_PassDxy  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_PassDxy");
   st.BS_TOF_FailPhi  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_FailPhi");
   st.BS_TOF_PassPhi  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_PassPhi");

   st.BS_BetaDiff  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_BetaDiff");
   st.BS_VertexTime = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_VertexTime");
   st.BS_SegSep = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegSep");
   st.BS_SegMinEtaSep_CSC = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinEtaSep_CSC");
   st.BS_SegMinEtaSep_DT = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinEtaSep_DT");
   st.BS_SegMinPhiSep = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinPhiSep");
   st.BS_SegMinEtaSep = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinEtaSep");
   st.BS_SegMinEtaSep_FailDz = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinEtaSep_FailDz");
   st.BS_SegMinEtaSep_PassDz = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegMinEtaSep_PassDz");
   st.BS_SegEta_MinEtaSep = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_SegEta_MinEtaSep");
   st.BS_IsTracker = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_IsTracker");
   st.BS_IsUpdated = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_IsUpdated");
   st.BS_MatchedStations = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_MatchedStations");
   st.BS_InvBetaErr = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_InvBetaErr");
   st.BS_PV = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PV");
   st.BS_dR_NVTrack = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_dR_NVTrack");

   st.BS_P      = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_P");
   st.BS_Pt     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt");
   st.BS_Pt_Bar     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_Bar");
   st.BS_Pt_For     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt_For");
   st.BS_QoverPt     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_QoverPt");
   st.BS_TOF    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF");
   st.BS_TOF_Bar    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_Bar");
   st.BS_TOF_For    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_For");
   st.BS_TOF_DT    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_DT");
   st.BS_TOF_CSC    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_CSC");

   st.BS_EtaP   = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaP");
   st.BS_EtaPt  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaPt");
   st.BS_EtaTOF  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaTOF");
   st.BS_EtaTime  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaTime");
   st.BS_EtaDz  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaDz");
   st.BS_PhiTime  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PhiTime");
   st.BS_DzTime  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_DzTime");
   st.BS_DzTime_DT  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_DzTime_DT");
   st.BS_DzTime_CSC  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_DzTime_CSC");
   st.BS_DzPt  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_DzPt");
   st.BS_PtTOF   = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtTOF");

   st.H_A  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A");
   st.H_B  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B");
   st.H_C  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C");
   st.H_D  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D");

   st.H_A_Cen  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_Cen");
   st.H_B_Cen  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_Cen");
   st.H_C_Cen  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_Cen");
   st.H_D_Cen  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_Cen");

   st.H_A_For  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_For");
   st.H_B_For  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_For");
   st.H_C_For  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_For");
   st.H_D_For  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_For");

   st.H_A_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_Low");
   st.H_B_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_Low");
   st.H_C_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_Low");
   st.H_D_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_Low");

   st.H_A_Cen_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_Cen_Low");
   st.H_B_Cen_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_Cen_Low");
   st.H_C_Cen_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_Cen_Low");
   st.H_D_Cen_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_Cen_Low");

   st.H_A_For_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_For_Low");
   st.H_B_For_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_For_Low");
   st.H_C_For_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_For_Low");
   st.H_D_For_Low  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_For_Low");

   for(int i=0; i<DzRegions; i++) {
     string Region = RegionNames[i];
     st.H_A_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_Syst_"+Region);
     st.H_B_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_Syst_"+Region);
     st.H_C_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_Syst_"+Region);
     st.H_D_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_Syst_"+Region);

     st.H_A_Cen_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_Cen_Syst_"+Region);
     st.H_B_Cen_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_Cen_Syst_"+Region);
     st.H_C_Cen_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_Cen_Syst_"+Region);
     st.H_D_Cen_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_Cen_Syst_"+Region);

     st.H_A_For_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_A_For_Syst_"+Region);
     st.H_B_For_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_B_For_Syst_"+Region);
     st.H_C_For_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_C_For_Syst_"+Region);
     st.H_D_For_Syst[i]  = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_D_For_Syst_"+Region);
   }
   st.H_DzCounts      = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_DzCounts");
   st.H_DzCounts_DT      = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_DzCounts_DT");
   st.H_DzCounts_CSC      = (TH1D*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/H_DzCounts_CSC");

   HistoFile->cd();
}

void stPlots_Clear(stPlots& st, bool WriteFirst=false, bool skipTree=true)
{
   if(WriteFirst){
     if(!skipTree) st.Tree->SetDirectory(st.Directory);
      st.Directory->Write();
   }
   delete st.Directory;
}


void stPlots_FillTree(stPlots& st, unsigned int Run, unsigned int Event, unsigned int Hscp, double Pt, double TOF, int MaxEntries=20000){
   if(MaxEntries>0 && st.Tree->GetEntries()>=MaxEntries)return;
   st.Tree_Run   = Run;
   st.Tree_Event = Event;
   st.Tree_Hscp  = Hscp;
   st.Tree_Pt    = Pt;
   st.Tree_TOF   = TOF;
   st.Tree->Fill();
}


void stPlots_Dump(stPlots& st, FILE* pFile, int CutIndex){
   fprintf(pFile,"---------- %10s ----------\n",st.Name.c_str());
   fprintf(pFile,"#Events                       = %4.2E\n",st.TotalE->GetBinContent(1       ));
   fprintf(pFile,"#Triggered Events             = %4.2E Eff=%4.3E\n",st.TotalTE->GetBinContent(1     ),st.TotalTE->GetBinContent(1      )/st.TotalE->GetBinContent(1       ));
   fprintf(pFile,"#Tracks                       = %4.2E Eff=%4.3E\n",st.Reconstructed->GetBinContent(1), st.Reconstructed->GetBinContent(1)/(2*st.TotalTE->GetBinContent(1)));
   fprintf(pFile,"#Tracks matched to trigger    = %4.2E Eff=%4.3E\n",st.TriggerMatch->GetBinContent(1), st.TriggerMatch->GetBinContent(1)/st.Reconstructed->GetBinContent(1));
   fprintf(pFile,"#Tracks with matched NV track = %4.2E Eff=%4.3E\n",st.NVTrack->GetBinContent(1), st.NVTrack->GetBinContent(1)/st.TriggerMatch->GetBinContent(1));
   fprintf(pFile,"#Tracks passing eta cut       = %4.2E Eff=%4.3E\n",st.Eta->GetBinContent(1), st.Eta->GetBinContent(1)/st.NVTrack->GetBinContent(1));
   fprintf(pFile,"#Tracks passing min pt  cut   = %4.2E Eff=%4.3E\n",st.MinPt->GetBinContent(1), st.MinPt->GetBinContent(1)/st.Eta->GetBinContent(1));
   fprintf(pFile,"#Tracks passing stations cut  = %4.2E Eff=%4.3E\n",st.Stations->GetBinContent(1), st.Stations->GetBinContent(1)/st.MinPt->GetBinContent(1));
   fprintf(pFile,"#Tracks passing nDOF cut      = %4.2E Eff=%4.3E\n",st.nDof->GetBinContent(1), st.nDof->GetBinContent(1)/st.Stations->GetBinContent(1));
   fprintf(pFile,"#Tracks passing tof error cut = %4.2E Eff=%4.3E\n",st.tofError->GetBinContent(1), st.tofError->GetBinContent(1)/st.nDof->GetBinContent(1));
   fprintf(pFile,"#Tracks passing Eta Sep cut   = %4.2E Eff=%4.3E\n",st.SegSep->GetBinContent(1), st.SegSep->GetBinContent(1)/st.tofError->GetBinContent(1));
   fprintf(pFile,"#Tracks passing Dxy cut        = %4.2E Eff=%4.3E\n",st.Dxy->GetBinContent(1), st.Dxy->GetBinContent(1)/st.SegSep->GetBinContent(1));
   fprintf(pFile,"#Tracks passing Dz cut        = %4.2E Eff=%4.3E\n",st.Dz->GetBinContent(1), st.Dz->GetBinContent(1)/st.Dxy->GetBinContent(1));
   fprintf(pFile,"#Tracks passing preselection  = %4.2E Eff=%4.3E\n",st.Preselected->GetBinContent(1), st.Preselected->GetBinContent(1)/st.Reconstructed->GetBinContent(1));

   fprintf(pFile,"#Tracks failing Dz     cuts           = %4.2E Eff=%4.3E\n",st.FailDz->GetBinContent(1), st.FailDz->GetBinContent(1)/st.Dxy->GetBinContent(1));

   fprintf(pFile,"#Tracks passing Pt      cuts   = %4.2E Eff=%4.3E\n",st.Pt   ->GetBinContent(CutIndex+1), st.Pt->GetBinContent(CutIndex+1) /st.Preselected->GetBinContent(1 ));
   fprintf(pFile,"#Tracks passing TOF Only cuts  = %4.2E Eff=%4.3E\n",st.TOFExclusive->GetBinContent(CutIndex+1), st.TOFExclusive->GetBinContent(CutIndex+1) /st.Preselected->GetBinContent(1));
   fprintf(pFile,"#Tracks passing TOF     cuts   = %4.2E Eff=%4.3E\n",st.TOF  ->GetBinContent(CutIndex+1), st.TOF->GetBinContent(CutIndex+1) /st.Pt->GetBinContent(CutIndex+1));
   fprintf(pFile,"#Tracks passing selection     = %4.2E Eff=%4.3E\n",st.TOF  ->GetBinContent(CutIndex+1), st.TOF->GetBinContent(CutIndex+1) /st.Reconstructed->GetBinContent(1)); 
   fprintf(pFile,"--------------------\n");
   fprintf(pFile,"HSCP Detection Efficiency Before Trigger                           Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(2*st.TotalE ->GetBinContent(1       )));
   fprintf(pFile,"HSCP Detection Efficiency After  Trigger                           Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(2*st.TotalTE->GetBinContent(1       )));
   fprintf(pFile,"HSCPTrack per HSCPEvent (with at least one HSCPTrack)              Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(  st.HSCPE  ->GetBinContent(CutIndex+1)));
   fprintf(pFile,"Events passing cuts = %4.2E                                     Eff=%4.3E\n",st.HSCPE->GetBinContent(CutIndex+1), st.HSCPE->GetBinContent(CutIndex+1) /(  st.TotalE->GetBinContent(1)));

   fprintf(pFile,"--------------------\n");
}


void stPlots_Draw(stPlots& st, std::string SavePath, std::string LegendTitle, unsigned int CutIndex)
{

   TH2** Histos = new TH2*[10];
   std::vector<std::string> legend;
   TCanvas* c1;

   char CutIndexStr[255];sprintf(CutIndexStr,"_%03i",CutIndex);

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_EtaP;                  legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p (GeV/c)", 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaP_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_EtaPt;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p_{T} (GeV/c)", -2.1,2.1, 0,1000, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaPt_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_EtaTOF;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "1/#beta", -2.2,2.2, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaTOF_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_EtaTime;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "Vertex Time", -2.2,2.2, -40,40, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaTime_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_EtaDz;                 legend.push_back("Before Cut");
   Histos[0]->Rebin2D(1,10);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "Dz [cm]", -2.2,2.2, -300,300, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaDz_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_PhiTime;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#phi", "Vertex Time", 0,0, -100,100, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PhiTime_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_DzTime;                 legend.push_back("Before Cut");
   Histos[0]->Rebin2D(10,1);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "dz [cm]", "Vertex Time", -300,300, -100,100, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DzTime_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_DzTime_DT;                 legend.push_back("Before Cut");
   Histos[0]->Rebin2D(10,1);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "dz [cm]", "Vertex Time", -150,150, -40,40, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DzTime_DT_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_DzTime_CSC;                 legend.push_back("Before Cut");
   Histos[0]->Rebin2D(10,1);
   //DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "dz [cm]", "Vertex Time", -300,300, -100,100, false);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "dz [cm]", "Vertex Time", -120,120, -40,40, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DzTime_CSC_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_DzPt;                 legend.push_back("Before Cut");
   Histos[0]->Rebin2D(10,1);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "dz [cm]", "p_{T} [GeV]", -300,300, 0,800, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DzPt_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_PtTOF;                  legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV/c)", "1/#beta", 0,1000, 0.5,1.5, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PtTOF_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_Dxy_Dz->Clone();                  legend.push_back("Before Cut");
   Histos[0]->Rebin2D(8,8);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "Dxy", "Dz", -100,100, -100,100, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DxyDz_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH2*)st.BS_SegEta_MinEtaSep->Clone();                  legend.push_back("Before Cut");
   Histos[0]->Rebin2D(20,20);
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "Eta", "dEta", -100,100, -100,100, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegEta_MinEtaSep_BS", true);
   delete c1;

   TH1** Histos1D = new TH1*[10];
   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.Eta_RegionA->GetXaxis()->SetRange(CutIndex+1,CutIndex+1); st.Eta_RegionA->RebinY(1);
   string nameA = "EtaA"+SavePath;
   Histos1D[0] = (TH1*)st.Eta_RegionA->ProjectionY(nameA.c_str(),CutIndex+1,CutIndex+1); legend.push_back("Region A");
   if(Histos1D[0]->Integral()>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral());
   st.Eta_RegionB->GetXaxis()->SetRange(CutIndex+1,CutIndex+1); st.Eta_RegionB->RebinY(1);
   string nameB = "EtaB"+SavePath;
   Histos1D[1] = (TH1*)st.Eta_RegionB->ProjectionY(nameB.c_str(),CutIndex+1,CutIndex+1); legend.push_back("Region B");
   if(Histos1D[1]->Integral()>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral());
   st.Eta_RegionC->GetXaxis()->SetRange(CutIndex+1,CutIndex+1); st.Eta_RegionC->RebinY(1);
   string nameC = "EtaC"+SavePath;
   Histos1D[2] = (TH1*)st.Eta_RegionC->ProjectionY(nameC.c_str(),CutIndex+1,CutIndex+1); legend.push_back("Region C");
   if(Histos1D[2]->Integral()>0) Histos1D[2]->Scale(1.0/Histos1D[2]->Integral());
   //st.Eta_RegionD->GetXaxis()->SetRange(CutIndex+1,CutIndex+1); st.Eta_RegionD->RebinY(1);
   //string nameD = "EtaD"+SavePath;
   //Histos1D[3] = (TH1*)st.Eta_RegionD->ProjectionY(nameD.c_str(),CutIndex+1,CutIndex+1); legend.push_back("Region D");
   //if(Histos1D[3]->Integral()>0) Histos1D[3]->Scale(1.0/Histos1D[3]->Integral());
   DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "#eta", "", 0,0, 0,0, false);
   DrawPreliminary(IntegratedLuminosity);
   DrawLegend((TObject**)Histos1D,legend,LegendTitle,"P", true);
   SaveCanvas(c1,SavePath,"Eta_Regions", false);
   delete Histos1D[0]; delete Histos1D[1]; delete Histos1D[2]; //delete Histos1D[3];
   delete c1;


      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
      Histos1D[0] = (TH1*)st.Beta_Gen;                                                  legend.push_back("Gen");
      Histos1D[1] = (TH1*)st.Beta_Triggered;                                            legend.push_back("Triggered");
      DrawSuperposedHistos((TH1**)Histos1D, legend,"HIST E1",  "#beta", "# HSCP", 0,0, 0,0);
      DrawLegend((TObject**)Histos1D,legend,"","P", 0.36, 0.92, 0.20, 0.04);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"Gen_Beta", true);
      delete c1;

      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
      Histos1D[0] = (TH1*)st.Beta_Triggered;                                            legend.push_back("Triggered");
      Histos1D[1] = (TH1*)st.Beta_Matched;                                              legend.push_back("Reconstructed");
      Histos1D[2] = (TH1*)st.Beta_PreselectedC;                                         legend.push_back("Preselected");
      Histos1D[3] = (TH1*)st.Beta_SelectedP->ProjectionY("A",CutIndex+1,CutIndex+1);    legend.push_back("p_{T}>Cut");
      Histos1D[4] = (TH1*)st.Beta_SelectedT->ProjectionY("B",CutIndex+1,CutIndex+1);    legend.push_back("ToF>Cut");
      DrawSuperposedHistos((TH1**)Histos1D, legend,"HIST E1",  "#beta", "# HSCP", 0,0, 0,0);
      DrawLegend((TObject**)Histos1D,legend,LegendTitle,"P", 0.36, 0.92, 0.20, 0.025);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"Beta", true);
      delete Histos1D[3]; delete Histos1D[4]; 
      delete c1;

      if(st.Name=="Cosmic") {
      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
      Histos1D[0] = (TH1*)st.BS_Pt_FailDxy->Clone();  Histos1D[0]->Rebin(4);              legend.push_back("abs(dxy)>10");
      if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
      Histos1D[1] = (TH1*)st.BS_Pt_PassDxy->Clone();  Histos1D[1]->Rebin(4);              legend.push_back("abs(dxy)<10");
      if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
      DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,400, 0.0005,1, false, true, false);
      DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"_Pt_Dxy_Comp", true);
      delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_GlobalPt_FailDxy->Clone();  Histos1D[0]->Rebin(4);              legend.push_back("abs(dxy)>10");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_GlobalPt_PassDxy->Clone();  Histos1D[1]->Rebin(4);              legend.push_back("abs(dxy)<10");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,600, 0,0, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_GlobalPt_Dxy_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Pt_FailDz->Clone(); Histos1D[0]->Rebin(4);            legend.push_back("abs(dz)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Pt_PassDz->Clone();  Histos1D[1]->Rebin(4);               legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,600, 0.0005,1, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Pt_Dz_Comp", false);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Pt_FailDz_CSC->Clone(); Histos1D[0]->Rebin(4);            legend.push_back("abs(dz)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Pt_PassDz_CSC->Clone();  Histos1D[1]->Rebin(4);               legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,600, 0.0005,1, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Pt_Dz_Comp_CSC", false);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Pt_FailDz_DT->Clone(); Histos1D[0]->Rebin(4);            legend.push_back("abs(dz)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Pt_PassDz_DT->Clone();  Histos1D[1]->Rebin(4);               legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,600, 0.0005,1, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Pt_Dz_Comp_DT", false);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Pt_FailPhi->Clone();  Histos1D[0]->Rebin(2);              legend.push_back("Sep<0.2");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Pt_PassPhi->Clone();  Histos1D[1]->Rebin(2);              legend.push_back("Sep>0.2");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Pt [GeV]", "arbitrary units", 0,600, 0.0005,1, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.69, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Pt_Sep_Comp", true);
     delete c1;

      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
      Histos1D[0] = (TH1*)st.BS_TOF_FailDxy->Clone();                                    legend.push_back("abs(dxy)>10");
      if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
      Histos1D[1] = (TH1*)st.BS_TOF_PassDxy->Clone();                                    legend.push_back("abs(dxy)<10");
      if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
      DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "1/#beta", "arbitrary units", -2, 5, 0.0005, 1, false, true, false);
      DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.3, 0.1);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"_TOF_Dxy_Comp", true);
      delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_TOF_FailDz->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_TOF_PassDz->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "1/#beta", "arbitrary units", -2, 4,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_TOF_Dz_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_TOF_FailDz_CSC->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_TOF_PassDz_CSC->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "1/#beta", "arbitrary units", -2, 4,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_TOF_Dz_CSC_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_TOF_FailDz_DT->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_TOF_PassDz_DT->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "1/#beta", "arbitrary units", -2, 4,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_TOF_Dz_DT_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Time_FailDz->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Time_PassDz->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Time", "arbitrary units", 0, 0,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Time_Dz_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Time_FailDz_CSC->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Time_PassDz_CSC->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Time", "arbitrary units", 0, 0,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Time_Dz_CSC_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Time_FailDz_DT->Clone();                                         legend.push_back("abs(z)>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Time_PassDz_DT->Clone();                                    legend.push_back("abs(dz)<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Time", "arbitrary units", 0, 0,0.0005 ,1 ,false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.89, 0.92, 0.2, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Time_Dz_DT_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_TOF_FailPhi->Clone();                                    legend.push_back("Sep<0.2");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_TOF_PassPhi->Clone();                                    legend.push_back("Sep>0.2");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "1/#beta", "arbitrary units", -2, 5, 0.0005, 1, false, true, false);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.3, 0.1);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_TOF_Sep_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Eta_FailDxy->Clone();  Histos1D[0]->Rebin(1);              legend.push_back("dxy>10");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Eta_PassDxy->Clone();  Histos1D[1]->Rebin(1);              legend.push_back("dxy<10");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, true);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.15, 0.08);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Eta_Dxy_Comp", true);
     delete c1;

     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos1D[0] = (TH1*)st.BS_Eta_FailDz->Clone();  Histos1D[0]->Rebin(1);              legend.push_back("dz>35");
     if(Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX()+1));
     Histos1D[1] = (TH1*)st.BS_Eta_PassDz->Clone();  Histos1D[1]->Rebin(1);              legend.push_back("dz<35");
     if(Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos1D, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, true);
     DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.15, 0.08);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"_Eta_Dz_Comp", true);
     delete c1;
      }
      char tmp[2048];
      c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
      Histos1D[0] = (TH1*)st.BS_Pt_Bar->Clone(); legend.push_back("|Eta|<0.9");
      Histos1D[0]->Rebin(1); if(Histos1D[0]->Integral(1, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(1, Histos1D[0]->GetNbinsX()+1));
      Histos1D[1] = (TH1*)st.BS_Pt_For->Clone(); legend.push_back("|Eta|>0.9");
      Histos1D[1]->Rebin(1); if(Histos1D[1]->Integral(1, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(1, Histos1D[1]->GetNbinsX()+1));
      sprintf(tmp,"Fraction of tracks/%2.0f GeV/c",Histos1D[0]->GetBinWidth(1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "",  "p_{T} (GeV/c)", tmp, 0, 1500, 0, 0, false, true, false);
      DrawLegend((TObject**)Histos1D,legend,LegendTitle,"P");
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"Pt_EtaComp_BS");
      delete c1;

      c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
      Histos1D[0] = (TH1*)st.BS_TOF_Bar->Clone(); legend.push_back("|Eta|<0.9");
      Histos1D[0]->Rebin(4); if(Histos1D[0]->Integral(1, Histos1D[0]->GetNbinsX()+1)>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral(1, Histos1D[0]->GetNbinsX()+1));
      Histos1D[1] = (TH1*)st.BS_TOF_For->Clone(); legend.push_back("|Eta|>0.9");
      Histos1D[1]->Rebin(4); if(Histos1D[1]->Integral(1, Histos1D[1]->GetNbinsX()+1)>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral(1, Histos1D[1]->GetNbinsX()+1));
      sprintf(tmp,"Fraction of tracks/%0.2f",Histos1D[0]->GetBinWidth(1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "",  "1/#beta", tmp, 0, 1500, 0, 0, false, true, false);
      DrawLegend((TObject**)Histos1D,legend,LegendTitle,"P");
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,"TOF_EtaComp_BS");
      delete c1;
}

void stPlots_DrawComparison(std::string SavePath, std::string LegendTitle, unsigned int CutIndex, stPlots* st1, stPlots* st2=NULL, stPlots* st3=NULL, stPlots* st4=NULL, stPlots* st5=NULL, stPlots* st6=NULL, stPlots* st7=NULL)
{ 
   char CutIndexStr[255];sprintf(CutIndexStr,"_%03i",CutIndex);

   char tmp[2048];
  std::vector<std::string> lg;
  std::vector<stPlots*> st;
  if(st1)st.push_back(st1); 
  if(st2)st.push_back(st2);   
  if(st3)st.push_back(st3);   
  if(st4)st.push_back(st4);
  if(st5)st.push_back(st5);
  if(st6)st.push_back(st6);
  if(st7)st.push_back(st7);

  std::vector<stSignal> signals;
  GetSignalDefinition(signals);
  for(unsigned int i=0;i<st.size();i++){
     int Index = -1;
     for(unsigned int s=0;s<signals.size();s++){
        if(signals[s].Name==st[i]->Name){Index=s;break;}
     }
     if(st[i]->Name=="MCTr"){lg.push_back("MC");}
     else if(st[i]->Name=="Data_Track"){lg.push_back("Global Muon");}
     else if(st[i]->Name=="Data_NoTrack"){lg.push_back("SA Only");}
     else if(st[i]->Name=="Data_Control"){lg.push_back("Cosmic Tagged");}
     else if(Index==-1){lg.push_back(st[i]->Name);}
     else{lg.push_back(signals[Index].Legend);}
  }
   TH1** Histos = new TH1*[10];
   std::vector<std::string> legend;
   TCanvas* c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->DistToGen->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dist To Gen HSCP", "arbitrary units", 0,0, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DistToGen", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->DistTrigger->Clone();        legend.push_back(lg[i]);  Histos[i]->Rebin(8);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));}
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dist To Trigger", "arbitrary units", 0, 5, 0,0, false, true);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"DistTrigger", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_FailDxy->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(1);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -70,70, 0.0001,2, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   //DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_FailDxy_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_PassDxy->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(1);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -70,70, 0.0001,2, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   //DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_PassDxy_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_V3D->Clone();      legend.push_back(lg[i]);  Histos[i]->Rebin(8); if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); 
   }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "V3D (cm)", "arbitrary units", 0, IPLimit, 0, 0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"V3D_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   int plots=0;
   for(unsigned int i=0;i<st.size();i++){
     if(lg[i]=="#tilde{g} 1000") {
       Histos[0] = (TH1*)st[i]->BS_Dz->Clone(); 
       if(Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)>0) Histos[0]->Scale(1.0/Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)); plots++;
  }
     if(lg[i]=="Cosmic") {
       Histos[1] = (TH1*)st[i]->BS_Dz->Clone(); if(Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)>0) Histos[1]->Scale(1.0/Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)); plots++;  }
   }
   if(plots==2) {
     double *glEff = new double[Histos[0]->GetNbinsX()/2+1];
     double *cosmicEff = new double[Histos[0]->GetNbinsX()/2 + 1];
     double *glEffErr = new double[Histos[0]->GetNbinsX()/2 + 1];
     double *cosmicEffErr = new double[Histos[0]->GetNbinsX()/2 + 1];

     glEff[0]=Histos[0]->GetBinContent(Histos[0]->GetNbinsX()/2)+Histos[0]->GetBinContent(Histos[0]->GetNbinsX()/2+1);
     cosmicEff[0]=Histos[1]->GetBinContent(Histos[1]->GetNbinsX()/2)+Histos[1]->GetBinContent(Histos[1]->GetNbinsX()/2+1);
     glEffErr[0]=Histos[0]->GetBinError(Histos[0]->GetNbinsX()/2)+Histos[0]->GetBinError(Histos[0]->GetNbinsX()/2+1);
     cosmicEffErr[0]=Histos[1]->GetBinError(Histos[1]->GetNbinsX()/2)+Histos[1]->GetBinError(Histos[1]->GetNbinsX()/2+1);

     unsigned int count=0;
     for(int i=Histos[0]->GetNbinsX()/2-1; i>-1; i--) {
       count++;
       int oppSide=Histos[0]->GetNbinsX()-i+1;
       glEff[count]=Histos[0]->GetBinContent(i) + Histos[0]->GetBinContent(oppSide) + glEff[count-1];
       cosmicEff[count]=Histos[1]->GetBinContent(i) + Histos[1]->GetBinContent(oppSide) + cosmicEff[count-1];
       glEffErr[count]=0;
       cosmicEffErr[count]=0;
     }

     TGraphErrors *gr1 = new TGraphErrors (count, glEff, cosmicEff, glEffErr, cosmicEffErr);
     gr1->GetXaxis()->SetTitle("Signal efficiency");
     gr1->GetYaxis()->SetTitle("Cosmic efficiency");
     gr1->Draw("APL");
     c1->SetGridx(true);
     c1->SetGridy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,"Dz_Eff_BS", true);
     for(unsigned int i=0;i<2;i++){delete Histos[i];}
     delete[] glEff; delete[] cosmicEff;
   }
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dxy->Clone();      legend.push_back(lg[i]); Histos[i]->Rebin(2); 
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dxy (cm)", tmp, -70, 70, 0.0001,0.6, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dxy_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(4);  
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -250, 250, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_BS", false);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_CSC->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(8);  
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -1*IPLimit, IPLimit, 0.000001,1, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_CSC_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_DT->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(8);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -1*IPLimit, IPLimit, 0.000001,1, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_DT_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_NoZed->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(2);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -1*IPLimit, IPLimit, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   //DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_NoZed_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Pt_FailPhi->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(8);  
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, 0, 400, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   SaveCanvas(c1,SavePath,"Pt_FailPhi_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_FailPhi->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(8);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -1*IPLimit, IPLimit, 0.0001,1, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.79, 0.92, 0.25, 0.08);
   c1->SetLogy(true);
   //DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_FailPhi_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   for(unsigned int i=0;i<st.size();i++){
     if(st[i]->Name!="Data_Track")continue;
     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos[0] = (TH1*)st[i]->BS_Dz_FailPhi->Clone();  Histos[0]->Rebin(2);              legend.push_back("Sep<0.2");
     if(Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)>0) Histos[0]->Scale(1.0/Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1));
     Histos[1] = (TH1*)st[i]->BS_Dz_PassPhi->Clone();  Histos[1]->Rebin(2);              legend.push_back("Sep>0.2");
     if(Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)>0) Histos[1]->Scale(1.0/Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, false);
     DrawLegend((TObject**)Histos,legend,"","P");
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,st[i]->Name + "_Dz_Sep_Comp", true);
     delete c1;
   }

   for(unsigned int i=0;i<st.size();i++){
     if(st[i]->Name!="Data_Track" && st[i]->Name!="Cosmic")continue;
     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos[0] = (TH1*)st[i]->BS_Dz_FailDxy->Clone();  Histos[0]->Rebin(1);              legend.push_back("dxy>10");
     if(Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)>0) Histos[0]->Scale(1.0/Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1));
     Histos[1] = (TH1*)st[i]->BS_Dz_PassDxy->Clone();  Histos[1]->Rebin(1);              legend.push_back("dxy<10");
     if(Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)>0) Histos[1]->Scale(1.0/Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, true);
     DrawLegend((TObject**)Histos,legend,"","P", 0.79, 0.92, 0.15, 0.08);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,st[i]->Name + "_Dz_Dxy_Comp", true);
     delete c1;
   }

   for(unsigned int i=0;i<st.size();i++){
     if(st[i]->Name!="Cosmic")continue;
     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos[0] = (TH1*)st[i]->BS_Eta_FailDxy->Clone();  Histos[0]->Rebin(1);              legend.push_back("dxy>10");
     if(Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)>0) Histos[0]->Scale(1.0/Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1));
     Histos[1] = (TH1*)st[i]->BS_Eta_PassDxy->Clone();  Histos[1]->Rebin(1);              legend.push_back("dxy<10");
     if(Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)>0) Histos[1]->Scale(1.0/Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, true);
     DrawLegend((TObject**)Histos,legend,"","P", 0.79, 0.92, 0.15, 0.08);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,st[i]->Name + "_Eta_Dxy_Comp", true);
     delete c1;
   }

   for(unsigned int i=0;i<st.size();i++){
     if(st[i]->Name!="Cosmic")continue;
     c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
     Histos[0] = (TH1*)st[i]->BS_Eta_FailDz->Clone();  Histos[0]->Rebin(1);              legend.push_back("dz>35");
     if(Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1)>0) Histos[0]->Scale(1.0/Histos[0]->Integral(0, Histos[0]->GetNbinsX()+1));
     Histos[1] = (TH1*)st[i]->BS_Eta_PassDz->Clone();  Histos[1]->Rebin(1);              legend.push_back("dz<35");
     if(Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1)>0) Histos[1]->Scale(1.0/Histos[1]->Integral(0, Histos[1]->GetNbinsX()+1));
     DrawSuperposedHistos((TH1**)Histos, legend,"E1",  "Dz [cm]", "arbitrary units", -100,100, 0,0, false, true, true);
     DrawLegend((TObject**)Histos,legend,"","P", 0.79, 0.92, 0.15, 0.08);
     c1->SetLogy(true);
     DrawPreliminary(IntegratedLuminosity);
     SaveCanvas(c1,SavePath,st[i]->Name + "_Eta_Dz_Comp", true);
     delete c1;
   }

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz_PassPhi->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(2);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, -1*IPLimit, IPLimit, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.82, 0.96, 0.16, 0.03);
   c1->SetLogy(true);
   //DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_PassPhi_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Dz->Clone();  legend.push_back(lg[i]); Histos[i]->Rebin(4);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
   }
   sprintf(tmp,"Fraction of tracks/%2.0f [cm]",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Dz (cm)", tmp, 0, 120, 0,0, false, false, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Dz_OneSide_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Chi2->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#chi^{2}/ndof", "arbitrary units", 0,15, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Chi2_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TNOH->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#NOH", "arbitrary units", 0,0, 0,0, false, true);
   //DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOH_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TNOH_Endcap->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#NOH", "arbitrary units", 0,0, 0,0, false, true);
   //DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOH_Endcap_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TNOH_Barrel->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#NOH", "arbitrary units", 0,0, 0,0, false, true);
   //DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOH_Barrel_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Eta->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#eta", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   //c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Eta_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Phi->Clone();        legend.push_back(lg[i]);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#eta", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   //c1->SetLogy(true);                                                                                                                                                            
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Phi_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_nDof->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "TOF_{nDof}", "arbitrary units", 0,0, 0.0005,1);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"nDof_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_VertexTime->Clone(); Histos[i]->Rebin(2);  legend.push_back(lg[i]); if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Time at vertex [ns]", "arbitrary units", -100,100, 0, 0);
   //DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"VertexTime_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_PV->Clone(); Histos[i]->Rebin(1);  legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Primary Vertices", "arbitrary units", 0,0, 0.0001,2);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PV_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_dR_NVTrack->Clone(); Histos[i]->Rebin(1);  legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dR", "arbitrary units", 0,0.4, 0.0001, 2, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"dR_NVTrack_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegSep->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(4);  
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dR to opp side segment", "arbitrary units", 0,2.5, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegSep_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_CSC->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(20);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dR to opp side segment", "arbitrary units", 0,0, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinEtaSep_CSC_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_DT->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(20);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dR to opp side segment", "arbitrary units", 0,0, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinEtaSep_DT_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinPhiSep->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(20);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dPhi to opp side segment", "arbitrary units", 0,0, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinPhiSep_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(10);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dEta to opp side segment", "arbitrary units", 0,0, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinEtaSep_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_FailDz->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(10);  
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dR to opp side segment", "arbitrary units", -0.5,0.5, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinEtaSep_FailDz_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_PassDz->Clone();  legend.push_back(lg[i]);  Histos[i]->Rebin(2);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "dEta to opp side segment", "arbitrary units", -0.5,0.5, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"SegMinEtaSep_PassDz_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_IsTracker->Clone(); legend.push_back(lg[i]);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Has Track (0 is no, 1 is yes)", "arbitrary units", 0,0, 0,1.2);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"HasTrack_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_IsUpdated->Clone(); legend.push_back(lg[i]);
     if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Updated (0 is no, 1 is yes)", "arbitrary units", 0,0, 0,1.2);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"IsUpdated_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_MatchedStations->Clone(); legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Muon stations", "arbitrary units", 0,5, 0.01,2);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"MatchedStations_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_InvBetaErr->Clone(); Histos[i]->Rebin(1);  legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Inverse Beta Err", "arbitrary units", 0,0.3, 0.0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"InvBetaErr_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Pterr->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} Err / p_{T}", "arbitrary units", 0,2, 0.001,0.5);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pterr_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_PterrSq->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} Err / (p_{T}*p_{T})", "arbitrary units", 0,0.01, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PterrSq_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Pt_All->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T}", "arbitrary units", 0,2000, 0,0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pt_All_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Pt_Global->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T}", "arbitrary units", 0,1500, 0.000005,1, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pt_Global_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_TOF_All->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} Err / p_{T}", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.89, 0.92, 0.16, 0.05);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_All_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Pt_FailDz->Clone();  Histos[i]->Rebin(2);  legend.push_back(lg[i]);  if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T}", "arbitrary units", 0,1000, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pt_FailDz_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Pt->Clone(); legend.push_back(lg[i]);  
     Histos[i]->Rebin(1); if(Histos[i]->Integral(1, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(1, Histos[i]->GetNbinsX()+1));}
   sprintf(tmp,"Fraction of tracks/%2.0f GeV/c",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "",  "p_{T} (GeV/c)", tmp, 0, 1500, 0, 0, false, true, false);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pt_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_QoverPt->Clone(); legend.push_back(lg[i]);
     Histos[i]->Rebin(8); if(Histos[i]->Integral(1, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(1, Histos[i]->GetNbinsX()+1));}
   sprintf(tmp,"Fraction of tracks/%2.0f GeV/c",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "",  "((q/p_{T})_{Gen} - (q/p_{T})_{Reco}) / (q/p_{T})_{Gen})", "Arbitrary units", -3, 5, 0.0005, 0.3, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.4);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"QoverPt_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_TOF; legend.push_back(lg[i]);
   if(Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1)>0) Histos[i]->Scale(1.0/Histos[i]->Integral(0, Histos[i]->GetNbinsX()+1));
}
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0,4, 0,0, false, true, true);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_Bar; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   //char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_Bar_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_For; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   //char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_For_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_DT; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   //char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_DT_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_CSC; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   //char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_CSC_BS");
   delete c1;
}
