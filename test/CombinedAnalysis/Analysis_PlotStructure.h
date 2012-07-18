#include "Analysis_Samples.h"


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
   float        Tree_I;
   float        Tree_TOF;
   float        Tree_Mass;

   TH2F*  Mass;
   TH2F*  MassTOF;
   TH2F*  MassComb;
   TH2F*  MaxEventMass;

   TH2F*  Mass_SystP;
   TH2F*  MassTOF_SystP;
   TH2F*  MassComb_SystP;
   TH2F*  MaxEventMass_SystP;

   TH2F*  Mass_SystI;
   TH2F*  MassTOF_SystI;
   TH2F*  MassComb_SystI;
   TH2F*  MaxEventMass_SystI;

   TH2F*  Mass_SystM;
   TH2F*  MassTOF_SystM;
   TH2F*  MassComb_SystM;
   TH2F*  MaxEventMass_SystM;

   TH2F*  Mass_SystT;
   TH2F*  MassTOF_SystT;
   TH2F*  MassComb_SystT;
   TH2F*  MaxEventMass_SystT;

  TH2F*  Mass_SystPU;
  TH2F*  MassTOF_SystPU;
  TH2F*  MassComb_SystPU;
  TH2F*  MaxEventMass_SystPU;

   TH1F* TotalE;
   TH1F* TotalEPU; 
   TH1F* TotalTE;
   TH1F* Total;
   TH1F* V3D; 
   TH1F* Chi2;  
   TH1F* Qual; 
   TH1F* TNOH;
   TH1F* TNOM;
   TH1F* nDof;
   TH1F* Pterr;
   TH1F* MPt; 
   TH1F* MI; 
   TH1F* MTOF; 
   TH1F* TIsol;
   TH1F* EIsol;
   TH1F* Pt;	
   TH1F* I;	
   TH1F* TOF;
   TH1F* HSCPE;

   TH1F* HSCPE_SystP;
   TH1F* HSCPE_SystI;
   TH1F* HSCPE_SystM;
   TH1F* HSCPE_SystT;
   TH1F* HSCPE_SystPU;

   TH1F* Beta_Gen;
   TH1F* Beta_GenCharged;
   TH1F* Beta_Triggered;
   TH1F* Beta_Matched;
   TH1F* Beta_PreselectedA;
   TH1F* Beta_PreselectedB;
   TH1F* Beta_PreselectedC;
   TH2F* Beta_SelectedP;
   TH2F* Beta_SelectedI;
   TH2F* Beta_SelectedT;

   TH1F*  BS_V3D;
   TH1F*  BS_Chi2;
   TH1F*  BS_Qual;
   TH1F*  BS_TNOH;
   TH1F*  BS_TNOHFraction;
   TH1F*  BS_Eta;
   TH1F*  BS_TNOM;
   TH1F*  BS_nDof;
   TH1F*  BS_Pterr;
   TH1F*  BS_MPt; 
   TH1F*  BS_MIs; 
   TH1F*  BS_MIm; 
   TH1F*  BS_MTOF;
   TH1F*  BS_TIsol;
   TH1F*  BS_EIsol;

  TH2F* AS_Eta_RegionA;
  TH2F* AS_Eta_RegionB;
  TH2F* AS_Eta_RegionC;
  TH2F* AS_Eta_RegionD;
  TH2F* AS_Eta_RegionE;
  TH2F* AS_Eta_RegionF;
  TH2F* AS_Eta_RegionG;
  TH2F* AS_Eta_RegionH;

   TH1F*  BS_P; 	   TH2F*  AS_P;
   TH1F*  BS_Pt;	   TH2F*  AS_Pt;
   TH1F*  BS_Is;	   TH2F*  AS_Is;
   TH1F*  BS_Im;           TH2F*  AS_Im;
   TH1F*  BS_TOF;          TH2F*  AS_TOF;
   TH1F*  BS_TOF_DT;
   TH1F*  BS_TOF_CSC;

   TH2F*  BS_EtaIs;        //TH3F*  AS_EtaIs;
   TH2F*  BS_EtaIm;        //TH3F*  AS_EtaIm;
   TH2F*  BS_EtaP;	   //TH3F*  AS_EtaP;
   TH2F*  BS_EtaPt;	   //TH3F*  AS_EtaPt;
   TH2F*  BS_EtaTOF;       //TH3F*  AS_EtaTOF;

   TH2F*  BS_PIs;	   TH3F*  AS_PIs;
   TH2F*  BS_PIm;          TH3F*  AS_PIm;
   TH2F*  BS_PtIs;         TH3F*  AS_PtIs;
   TH2F*  BS_PtIm;         TH3F*  AS_PtIm;
   TH2F*  BS_TOFIs;        TH3F*  AS_TOFIs;  
   TH2F*  BS_TOFIm;        TH3F*  AS_TOFIm;   

  //Muon only plots
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
  TH1F* tofError;
  TH1F* Dxy;
  TH1F* Dz;
  TH1F* TOFExclusive;
  TH1F* SegSep;
  TH1F* FailDz;
  TH1F* FailDz_DT;
  TH1F* FailDz_CSC;

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
   TH1F*  BS_TNOH_Barrel;
   TH1F*  BS_TNOH_Endcap;
  TH1F*  BS_Phi;
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

  TH1F*  BS_Pt_Bar;
  TH1F*  BS_Pt_For;
  TH1F*  BS_QoverPt;
  TH1F*  BS_TOF_Bar;
  TH1F*  BS_TOF_For;

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
  TH1D* H_E;
  TH1D* H_F;
  TH1D* H_G;
  TH1D* H_H;

  TH1D* H_A_Cen;
  TH1D* H_B_Cen;
  TH1D* H_C_Cen;
  TH1D* H_D_Cen;
  TH1D* H_E_Cen;
  TH1D* H_F_Cen;
  TH1D* H_G_Cen;
  TH1D* H_H_Cen;

  TH1D* H_A_For;
  TH1D* H_B_For;
  TH1D* H_C_For;
  TH1D* H_D_For;
  TH1D* H_E_For;
  TH1D* H_F_For;
  TH1D* H_G_For;
  TH1D* H_H_For;

  TH1D*  HCuts_Pt;
  TH1D*  HCuts_Is;
  TH1D*  HCuts_TOF;

  TH1D*  Hist_Pt ;
  TH1D*  Hist_Is  ;
  TH1D*  Hist_TOF;

  TH3D*  Pred_EtaP ;
  TH2D*  Pred_I    ;
  TH2D*  Pred_TOF  ;
  TH2D*  Pred_EtaB;
  TH2D*  Pred_EtaS;
  TH2D*  Pred_EtaS2;

  TH2D*  RegionD_P;
  TH2D*  RegionD_I;
  TH2D*  RegionD_TOF;

  TH1D* H_A_Flip;
  TH1D* H_B_Flip;
  TH1D* H_C_Flip;
  TH1D* H_D_Flip;
  TH1D* H_E_Flip;
  TH1D* H_F_Flip;
  TH1D* H_G_Flip;
  TH1D* H_H_Flip;

  TH1D* H_A_Cen_Flip;
  TH1D* H_B_Cen_Flip;
  TH1D* H_C_Cen_Flip;
  TH1D* H_D_Cen_Flip;
  TH1D* H_E_Cen_Flip;
  TH1D* H_F_Cen_Flip;
  TH1D* H_G_Cen_Flip;
  TH1D* H_H_Cen_Flip;

  TH1D* H_A_For_Flip;
  TH1D* H_B_For_Flip;
  TH1D* H_C_For_Flip;
  TH1D* H_D_For_Flip;
  TH1D* H_E_For_Flip;
  TH1D* H_F_For_Flip;
  TH1D* H_G_For_Flip;
  TH1D* H_H_For_Flip;

  TH3D*  Pred_EtaP_Flip ;
  TH2D*  Pred_I_Flip    ;
  TH2D*  Pred_TOF_Flip  ;
  TH2D*  Pred_EtaB_Flip;
  TH2D*  Pred_EtaS_Flip;
  TH2D*  Pred_EtaS2_Flip;

  TH2D* H_D_DzSidebands;
  TH2D* H_D_DzSidebands_DT;
  TH2D* H_D_DzSidebands_CSC;
  /*
  TH1D* H_A_Syst[DzRegions];
  TH1D* H_B_Syst[DzRegions];
  TH1D* H_C_Syst[DzRegions];
  TH1D* H_D_Syst[DzRegions];
  TH1D* H_E_Syst[DzRegions];
  TH1D* H_F_Syst[DzRegions];
  TH1D* H_G_Syst[DzRegions];
  TH1D* H_H_Syst[DzRegions];

  TH1D* H_A_Cen_Syst[DzRegions];
  TH1D* H_B_Cen_Syst[DzRegions];
  TH1D* H_C_Cen_Syst[DzRegions];
  TH1D* H_D_Cen_Syst[DzRegions];
  TH1D* H_E_Cen_Syst[DzRegions];
  TH1D* H_F_Cen_Syst[DzRegions];
  TH1D* H_G_Cen_Syst[DzRegions];
  TH1D* H_H_Cen_Syst[DzRegions];

  TH1D* H_A_For_Syst[DzRegions];
  TH1D* H_B_For_Syst[DzRegions];
  TH1D* H_C_For_Syst[DzRegions];
  TH1D* H_D_For_Syst[DzRegions];
  TH1D* H_E_For_Syst[DzRegions];
  TH1D* H_F_For_Syst[DzRegions];
  TH1D* H_G_For_Syst[DzRegions];
  TH1D* H_H_For_Syst[DzRegions];

  TH1D* H_DzCounts;
  TH1D* H_DzCounts_DT;
  TH1D* H_DzCounts_CSC;
  */
  TH1D*  CtrlPt_S1_Is;
  TH1D*  CtrlPt_S2_Is;
  TH1D*  CtrlPt_S3_Is;
  TH1D*  CtrlPt_S4_Is;

  TH1D*  CtrlIs_S1_TOF;
  TH1D*  CtrlIs_S2_TOF;
  TH1D*  CtrlIs_S3_TOF;
  TH1D*  CtrlIs_S4_TOF;

  TH1D*  CtrlIm_S1_TOF;
  TH1D*  CtrlIm_S2_TOF;
  TH1D*  CtrlIm_S3_TOF;
  TH1D*  CtrlIm_S4_TOF;

  TH1D*  CtrlPt_S1_Im;
  TH1D*  CtrlPt_S2_Im;
  TH1D*  CtrlPt_S3_Im;
  TH1D*  CtrlPt_S4_Im;

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

};

void stPlots_Init(TFile* HistoFile, stPlots& st, std::string BaseName, unsigned int NCuts, bool SkipSelectionPlot=false, bool isSignal=false)
{
   st.SelPlot = !SkipSelectionPlot;
   st.Name = BaseName;
   st.NCuts = NCuts;

   std::string Name;
   Name = BaseName;               st.Directory = HistoFile->mkdir(Name.c_str(), Name.c_str()); 
   st.Directory->cd();

   Name = "TotalE";   st.TotalE  = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "TotalEPU"; st.TotalEPU= new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "TotalTE";  st.TotalTE = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Total";    st.Total   = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "V3D";      st.V3D     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Chi2";     st.Chi2    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Qual";     st.Qual    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "TNOH";     st.TNOH    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "TNOM";     st.TNOM    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "nDof";     st.nDof    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Pterr";    st.Pterr   = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "TIsol";    st.TIsol   = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "EIsol";    st.EIsol   = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "MPt";      st.MPt     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "MI";       st.MI      = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "MTOF";     st.MTOF    = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);     
   Name = "Pt";       st.Pt      = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     
   Name = "I";        st.I       = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     
   Name = "TOF";      st.TOF     = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     
   Name = "HSCPE";    st.HSCPE   = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);     

   Name = "HSCPE_SystP";    st.HSCPE_SystP  = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystI";    st.HSCPE_SystI  = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystM";    st.HSCPE_SystM  = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystT";    st.HSCPE_SystT  = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "HSCPE_SystPU";   st.HSCPE_SystPU = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);


   Name = "Mass";     st.Mass     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass    ->Sumw2();
   Name = "MassTOF";  st.MassTOF  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF ->Sumw2();
   Name = "MassComb"; st.MassComb = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb->Sumw2();
   Name = "MaxEventMass";     st.MaxEventMass     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MaxEventMass    ->Sumw2();

   Name = "Mass_SystP";     st.Mass_SystP     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass_SystP    ->Sumw2();
   Name = "MassTOF_SystP";  st.MassTOF_SystP  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF_SystP ->Sumw2();
   Name = "MassComb_SystP"; st.MassComb_SystP = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb_SystP->Sumw2();
   Name = "MaxEventMass_SystP";     st.MaxEventMass_SystP = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);st.MaxEventMass_SystP->Sumw2();

   Name = "Mass_SystI";     st.Mass_SystI     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass_SystI    ->Sumw2();
   Name = "MassTOF_SystI";  st.MassTOF_SystI  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF_SystI ->Sumw2();
   Name = "MassComb_SystI"; st.MassComb_SystI = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb_SystI->Sumw2();
   Name = "MaxEventMass_SystI";     st.MaxEventMass_SystI = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);st.MaxEventMass_SystI->Sumw2();

   Name = "Mass_SystM";     st.Mass_SystM     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass_SystM    ->Sumw2();
   Name = "MassTOF_SystM";  st.MassTOF_SystM  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF_SystM ->Sumw2();
   Name = "MassComb_SystM"; st.MassComb_SystM = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb_SystM->Sumw2();
   Name = "MaxEventMass_SystM";     st.MaxEventMass_SystM = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);st.MaxEventMass_SystM->Sumw2();

   Name = "Mass_SystT";     st.Mass_SystT     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass_SystT    ->Sumw2();
   Name = "MassTOF_SystT";  st.MassTOF_SystT  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF_SystT ->Sumw2();
   Name = "MassComb_SystT"; st.MassComb_SystT = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb_SystT->Sumw2();
   Name = "MaxEventMass_SystT";     st.MaxEventMass_SystT = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);st.MaxEventMass_SystT->Sumw2();

   Name = "Mass_SystPU";    st.Mass_SystPU     = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.Mass_SystPU    ->Sumw2();
   Name = "MassTOF_SystPU"; st.MassTOF_SystPU  = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassTOF_SystPU ->Sumw2();
   Name = "MassComb_SystPU";st.MassComb_SystPU = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);   st.MassComb_SystPU->Sumw2();
   Name = "MaxEventMass_SystPU";  st.MaxEventMass_SystPU = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts, MassNBins, 0, MassHistoUpperBound);st.MaxEventMass_SystPU->Sumw2();

   if(SkipSelectionPlot)return;
   Name = "Beta_Gen"         ; st.Beta_Gen         = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Gen         ->Sumw2();
   Name = "Beta_GenChaged"   ; st.Beta_GenCharged  = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_GenCharged  ->Sumw2();
   Name = "Beta_Triggered"   ; st.Beta_Triggered   = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Triggered   ->Sumw2();
   Name = "Beta_Matched"     ; st.Beta_Matched     = new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_Matched     ->Sumw2();
   Name = "Beta_PreselectedA"; st.Beta_PreselectedA= new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_PreselectedA->Sumw2();
   Name = "Beta_PreselectedB"; st.Beta_PreselectedB= new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_PreselectedB->Sumw2();
   Name = "Beta_PreselectedC"; st.Beta_PreselectedC= new TH1F(Name.c_str(), Name.c_str(),                 20, 0,  1);  st.Beta_PreselectedC->Sumw2();
   Name = "Beta_SelectedP"   ; st.Beta_SelectedP   = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts,   20, 0,  1);  st.Beta_SelectedP   ->Sumw2();
   Name = "Beta_SelectedI"   ; st.Beta_SelectedI   = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts,   20, 0,  1);  st.Beta_SelectedI   ->Sumw2();
   Name = "Beta_SelectedT"   ; st.Beta_SelectedT   = new TH2F(Name.c_str(), Name.c_str(),NCuts,0,NCuts,   20, 0,  1);  st.Beta_SelectedT   ->Sumw2();

   Name = "BS_V3D"  ; st.BS_V3D   = new TH1F(Name.c_str(), Name.c_str(),  20,  0,  5);                st.BS_V3D->Sumw2();
   Name = "BS_Chi2" ; st.BS_Chi2  = new TH1F(Name.c_str(), Name.c_str(),  20,  0,  5);                st.BS_Chi2->Sumw2();
   Name = "BS_Qual" ; st.BS_Qual  = new TH1F(Name.c_str(), Name.c_str(),  20,  0, 20);                st.BS_Qual->Sumw2();
   Name = "BS_TNOH" ; st.BS_TNOH  = new TH1F(Name.c_str(), Name.c_str(),  50,  0,  40);                st.BS_TNOH->Sumw2();
   Name = "BS_TNOHFraction" ; st.BS_TNOHFraction  = new TH1F(Name.c_str(), Name.c_str(),  50,  0,  1);                st.BS_TNOHFraction->Sumw2();
   Name = "BS_Eta" ; st.BS_Eta  = new TH1F(Name.c_str(), Name.c_str(),  50,  -2.6,  2.6);                st.BS_Eta->Sumw2();
   Name = "BS_TNOM" ; st.BS_TNOM  = new TH1F(Name.c_str(), Name.c_str(),  40,  0, 40);                st.BS_TNOM->Sumw2();
   Name = "BS_nDof" ; st.BS_nDof  = new TH1F(Name.c_str(), Name.c_str(),  20,  0, 40);                st.BS_nDof->Sumw2();
   Name = "BS_PtErr"; st.BS_Pterr = new TH1F(Name.c_str(), Name.c_str(),  40,  0,  1);                st.BS_Pterr->Sumw2();
   Name = "BS_MPt"  ; st.BS_MPt   = new TH1F(Name.c_str(), Name.c_str(),  50,  0, PtHistoUpperBound); st.BS_MPt->Sumw2();
   Name = "BS_MIs"  ; st.BS_MIs   = new TH1F(Name.c_str(), Name.c_str(),  50,  0, dEdxS_UpLim);       st.BS_MIs->Sumw2();
   Name = "BS_MIm"  ; st.BS_MIm   = new TH1F(Name.c_str(), Name.c_str(),  100,  0, dEdxM_UpLim);       st.BS_MIm->Sumw2();
   Name = "BS_MTOF" ; st.BS_MTOF  = new TH1F(Name.c_str(), Name.c_str(),  50, -2, 5);                 st.BS_MTOF->Sumw2();
   Name = "BS_TIsol"; st.BS_TIsol = new TH1F(Name.c_str(), Name.c_str(),  25,  0, 100);               st.BS_TIsol->Sumw2();
   Name = "BS_EIsol"; st.BS_EIsol = new TH1F(Name.c_str(), Name.c_str(),  25,  0, 1.5);               st.BS_EIsol->Sumw2();
   Name = "BS_P"    ; st.BS_P     = new TH1F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound); st.BS_P->Sumw2();
   Name = "BS_Pt"   ; st.BS_Pt    = new TH1F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound); st.BS_Pt->Sumw2();
   Name = "BS_Is"   ; st.BS_Is    = new TH1F(Name.c_str(), Name.c_str(),                   100, 0, dEdxS_UpLim);       st.BS_Is->Sumw2();
   Name = "BS_Im"   ; st.BS_Im    = new TH1F(Name.c_str(), Name.c_str(),                   100, 3, dEdxM_UpLim);       st.BS_Im->Sumw2();
   Name = "BS_TOF"  ; st.BS_TOF   = new TH1F(Name.c_str(), Name.c_str(),                   150, 1, 5);                 st.BS_TOF->Sumw2();
   Name = "BS_TOF_DT"  ; st.BS_TOF_DT   = new TH1F(Name.c_str(), Name.c_str(),                   150, 1, 5);                 st.BS_TOF_DT->Sumw2();
   Name = "BS_TOF_CSC"  ; st.BS_TOF_CSC   = new TH1F(Name.c_str(), Name.c_str(),                   150, 1, 5);                 st.BS_TOF_CSC->Sumw2();

   Name = "AS_Eta_RegionA" ; st.AS_Eta_RegionA  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionA->Sumw2();
   Name = "AS_Eta_RegionB" ; st.AS_Eta_RegionB  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionB->Sumw2();
   Name = "AS_Eta_RegionC" ; st.AS_Eta_RegionC  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionC->Sumw2();
   Name = "AS_Eta_RegionD" ; st.AS_Eta_RegionD  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionD->Sumw2();
   Name = "AS_Eta_RegionE" ; st.AS_Eta_RegionE  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionE->Sumw2();
   Name = "AS_Eta_RegionF" ; st.AS_Eta_RegionF  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionF->Sumw2();
   Name = "AS_Eta_RegionG" ; st.AS_Eta_RegionG  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionG->Sumw2();
   Name = "AS_Eta_RegionH" ; st.AS_Eta_RegionH  = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts,  50,  -2.6,  2.6);           st.AS_Eta_RegionH->Sumw2();

   Name = "AS_P"    ; st.AS_P     = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound); st.AS_P->Sumw2();
   Name = "AS_Pt"   ; st.AS_Pt    = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound); st.AS_Pt->Sumw2();
   Name = "AS_Is"   ; st.AS_Is    = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, dEdxS_UpLim);       st.AS_Is->Sumw2();
   Name = "AS_Im"   ; st.AS_Im    = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, dEdxM_UpLim);       st.AS_Im->Sumw2();
   Name = "AS_TOF"  ; st.AS_TOF   = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 1, 5);                 st.AS_TOF->Sumw2();


   Name = "BS_EtaIs"; st.BS_EtaIs = new TH2F(Name.c_str(), Name.c_str(),                   50,-3, 3, 50, 0, dEdxS_UpLim);
   Name = "BS_EtaIm"; st.BS_EtaIm = new TH2F(Name.c_str(), Name.c_str(),                   50,-3, 3, 50, 2.8, dEdxM_UpLim);
   Name = "BS_EtaP" ; st.BS_EtaP  = new TH2F(Name.c_str(), Name.c_str(),                   50,-3, 3, 50, 0, PtHistoUpperBound);
   Name = "BS_EtaPt"; st.BS_EtaPt = new TH2F(Name.c_str(), Name.c_str(),                   50,-3, 3, 50, 0, PtHistoUpperBound);
   Name = "BS_EtaTOF" ; st.BS_EtaTOF  = new TH2F(Name.c_str(), Name.c_str(),                   50,-3, 3, 50, 0, 3);
   Name = "BS_PIs"  ; st.BS_PIs   = new TH2F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound, 50, 0, dEdxS_UpLim);
   Name = "BS_PIm"  ; st.BS_PIm   = new TH2F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound, 50, 0, dEdxM_UpLim);
   Name = "BS_PtIs" ; st.BS_PtIs  = new TH2F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound, 50, 0, dEdxS_UpLim);
   Name = "BS_PtIm" ; st.BS_PtIm  = new TH2F(Name.c_str(), Name.c_str(),                   50, 0, PtHistoUpperBound, 50, 0, dEdxM_UpLim);
   Name = "BS_TOFIs"; st.BS_TOFIs = new TH2F(Name.c_str(), Name.c_str(),                   100, 1, 5, 100, 0, dEdxS_UpLim);
   Name = "BS_TOFIm"; st.BS_TOFIm = new TH2F(Name.c_str(), Name.c_str(),                   100, 1, 5, 100, 0, dEdxM_UpLim);

//   Name = "AS_EtaIs"; st.AS_EtaIs = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, dEdxS_UpLim);
//   Name = "AS_EtaIm"; st.AS_EtaIm = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, dEdxM_UpLim);
//   Name = "AS_EtaP" ; st.AS_EtaP  = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
//   Name = "AS_EtaPt"; st.AS_EtaPt = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
//   Name = "AS_EtaTOF"; st.AS_EtaTOF = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, 3);
   Name = "AS_PIs"  ; st.AS_PIs   = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound, 50, 0, dEdxS_UpLim);
   Name = "AS_PIm"  ; st.AS_PIm   = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound, 50, 0, dEdxM_UpLim);
   Name = "AS_PtIs" ; st.AS_PtIs  = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound, 50, 0, dEdxS_UpLim);
   Name = "AS_PtIm" ; st.AS_PtIm  = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 0, PtHistoUpperBound, 50, 0, dEdxM_UpLim);
   Name = "AS_TOFIs"; st.AS_TOFIs = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 1, 5, 50, 0, dEdxS_UpLim);
   Name = "AS_TOFIm"; st.AS_TOFIm = new TH3F(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50, 1, 5, 50, 0, dEdxM_UpLim);

   //Muon only plots
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
   Name = "TOFExclusive"; st.TOFExclusive     = new TH1F(Name.c_str(), Name.c_str(),  NCuts, 0,  NCuts);
   Name = "SegSep";      st.SegSep     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz";      st.FailDz     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz_DT";      st.FailDz_DT     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);
   Name = "FailDz_CSC";      st.FailDz_CSC     = new TH1F(Name.c_str(), Name.c_str(),  1    , 0,  1);

   Name = "Pt_Gen"         ; st.Pt_Gen         = new TH1F(Name.c_str(), Name.c_str(),                 100, 0,  1000);  st.Pt_Gen         ->Sumw2();
   Name = "DistToGen"  ; st.DistToGen  = new TH1F(Name.c_str(), Name.c_str(),  200,  0,  2);                st.DistToGen->Sumw2();
   Name = "DistTrigger"; st.DistTrigger  = new TH1F(Name.c_str(), Name.c_str(),  1000,  0,  10);                st.DistTrigger->Sumw2();
   Name = "BS_Dz_FailDxy"  ; st.BS_Dz_FailDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -100,  100);                st.BS_Dz_FailDxy->Sumw2();
   Name = "BS_Dz_PassDxy"  ; st.BS_Dz_PassDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -100,  100);                st.BS_Dz_PassDxy->Sumw2();
   Name = "BS_Eta_FailDxy"  ; st.BS_Eta_FailDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_FailDxy->Sumw2();
   Name = "BS_Eta_PassDxy"  ; st.BS_Eta_PassDxy   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_PassDxy->Sumw2();
   Name = "BS_Eta_FailDz"  ; st.BS_Eta_FailDz   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_FailDz->Sumw2();
   Name = "BS_Eta_PassDz"  ; st.BS_Eta_PassDz   = new TH1F(Name.c_str(), Name.c_str(),  100,  -2.1,  2.1);                st.BS_Eta_PassDz->Sumw2();
   Name = "BS_Dxy";st.BS_Dxy   = new TH1F(Name.c_str(), Name.c_str(),  2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dxy->Sumw2();
   Name = "BS_Dz"; st.BS_Dz   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz->Sumw2();
   Name = "BS_Dz_CSC"; st.BS_Dz_CSC = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_CSC->Sumw2();
   Name = "BS_Dz_DT"; st.BS_Dz_DT=new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_DT->Sumw2();
   Name = "BS_Dz_FailPhi"; st.BS_Dz_FailPhi   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_FailPhi->Sumw2();
   Name = "BS_Dz_PassPhi"; st.BS_Dz_PassPhi   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_PassPhi->Sumw2();
   Name = "BS_Dz_NoZed"; st.BS_Dz_NoZed   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_NoZed->Sumw2();
   Name = "BS_Dz_Def"; st.BS_Dz_Def   = new TH1F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,  -IPHistoUpperBound,  IPHistoUpperBound); st.BS_Dz_Def->Sumw2();
   Name = "BS_Dxy_Dz";st.BS_Dxy_Dz = new TH2F(Name.c_str(), Name.c_str(),2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound,2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound); st.BS_Dxy_Dz->Sumw2();
   Name = "BS_TNOH_Barrel" ; st.BS_TNOH_Barrel  = new TH1F(Name.c_str(), Name.c_str(),  61,  -0.5,  60.5);                st.BS_TNOH_Barrel->Sumw2();
   Name = "BS_TNOH_Endcap" ; st.BS_TNOH_Endcap  = new TH1F(Name.c_str(), Name.c_str(),  61,  -0.5,  60.5);                st.BS_TNOH_Endcap->Sumw2();
   Name = "BS_Phi" ; st.BS_Phi  = new TH1F(Name.c_str(), Name.c_str(),  50,  -3.14,  3.14);              st.BS_Phi->Sumw2();
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

   Name = "BS_Pt_Bar"   ; st.BS_Pt_Bar    = new TH1F(Name.c_str(), Name.c_str(),                   160, 70, PtHistoUpperBound); st.BS_Pt_Bar->Sumw2();
   Name = "BS_Pt_For"   ; st.BS_Pt_For    = new TH1F(Name.c_str(), Name.c_str(),                   160, 70, PtHistoUpperBound); st.BS_Pt_For->Sumw2();
   Name = "BS_QoverPt"   ; st.BS_QoverPt    = new TH1F(Name.c_str(), Name.c_str(),                   300, -5, 3); st.BS_QoverPt->Sumw2();
   Name = "BS_TOF_Bar"  ; st.BS_TOF_Bar   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_Bar->Sumw2();
   Name = "BS_TOF_For"  ; st.BS_TOF_For   = new TH1F(Name.c_str(), Name.c_str(),                   150, 0, 2);                 st.BS_TOF_For->Sumw2();

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

   Name = "BS_EtaTime" ; st.BS_EtaTime  = new TH2F(Name.c_str(), Name.c_str(),                   100,-3, 3, 100, -200, 200);
   Name = "BS_EtaDz" ; st.BS_EtaDz  = new TH2F(Name.c_str(), Name.c_str(),  100,-3, 3, 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound);
   Name = "BS_PhiTime" ; st.BS_PhiTime  = new TH2F(Name.c_str(), Name.c_str(),                   50,-3.14, 3.14, 100, -200, 200);
   Name = "BS_DzTime" ; st.BS_DzTime  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzTime_DT" ; st.BS_DzTime_DT  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzTime_CSC" ; st.BS_DzTime_CSC  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 100, -200, 200);
   Name = "BS_DzPt" ; st.BS_DzPt  = new TH2F(Name.c_str(), Name.c_str(), 2*IPHistoUpperBound,-IPHistoUpperBound, IPHistoUpperBound, 150, 0, PtHistoUpperBound);
   Name = "BS_PtTOF" ; st.BS_PtTOF= new TH2F(Name.c_str(), Name.c_str(),                   150, 0, PtHistoUpperBound, 100, 0.5, 1.5);

   //Background prediction histograms don't need to be made for signal
   if(!isSignal) {
   Name = "H_A"; st.H_A = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A->Sumw2();
   Name = "H_B"; st.H_B = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B->Sumw2();
   Name = "H_C"; st.H_C = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C->Sumw2();
   Name = "H_D"; st.H_D = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D->Sumw2();
   Name = "H_E"; st.H_E = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E->Sumw2();
   Name = "H_F"; st.H_F = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F->Sumw2();
   Name = "H_G"; st.H_G = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G->Sumw2();
   Name = "H_H"; st.H_H = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H->Sumw2();

   Name = "H_A_Cen"; st.H_A_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen->Sumw2();
   Name = "H_B_Cen"; st.H_B_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen->Sumw2();
   Name = "H_C_Cen"; st.H_C_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen->Sumw2();
   Name = "H_D_Cen"; st.H_D_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen->Sumw2();
   Name = "H_E_Cen"; st.H_E_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_Cen->Sumw2();
   Name = "H_F_Cen"; st.H_F_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_Cen->Sumw2();
   Name = "H_G_Cen"; st.H_G_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_Cen->Sumw2();
   Name = "H_H_Cen"; st.H_H_Cen = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_Cen->Sumw2();

   Name = "H_A_For"; st.H_A_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For->Sumw2();
   Name = "H_B_For"; st.H_B_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For->Sumw2();
   Name = "H_C_For"; st.H_C_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For->Sumw2();
   Name = "H_D_For"; st.H_D_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For->Sumw2();
   Name = "H_E_For"; st.H_E_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_For->Sumw2();
   Name = "H_F_For"; st.H_F_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_For->Sumw2();
   Name = "H_G_For"; st.H_G_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_For->Sumw2();
   Name = "H_H_For"; st.H_H_For = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_For->Sumw2();

   Name = "Hist_Is"; st.Hist_Is = new TH1D(Name.c_str(), Name.c_str() ,200,0,dEdxS_UpLim); st.Hist_Is->Sumw2();
   Name = "Hist_Pt"; st.Hist_Pt = new TH1D(Name.c_str(), Name.c_str() ,200,0,PtHistoUpperBound); st.Hist_Pt->Sumw2();
   Name = "Hist_TOF"; st.Hist_TOF = new TH1D(Name.c_str(), Name.c_str() ,200,-10,20); st.Hist_TOF->Sumw2();
   Name = "Pred_I"; st.Pred_I = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 200,GlobalMinIm,dEdxM_UpLim); st.Pred_I->Sumw2();
   Name = "Pred_EtaB"; st.Pred_EtaB = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaB->Sumw2();
   Name = "Pred_EtaS"; st.Pred_EtaS = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaS->Sumw2();
   Name = "Pred_EtaS2"; st.Pred_EtaS2 = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaS2->Sumw2();
   Name = "Pred_EtaP"; st.Pred_EtaP = new TH3D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50, -3, 3, 200,GlobalMinPt,PtHistoUpperBound); st.Pred_EtaP->Sumw2();
   Name = "Pred_TOF"; st.Pred_TOF = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts,   200,GlobalMinTOF,5); st.Pred_TOF->Sumw2();
   Name = "RegionD_I"; st.RegionD_I = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 200,GlobalMinIm,dEdxM_UpLim); st.RegionD_I->Sumw2();
   Name = "RegionD_P"; st.RegionD_P = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 200,GlobalMinPt,PtHistoUpperBound); st.RegionD_P->Sumw2();
   Name = "RegionD_TOF"; st.RegionD_TOF = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 200,GlobalMinTOF,5); st.RegionD_TOF->Sumw2();
   
   Name = "H_A_Flip"; st.H_A_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Flip->Sumw2();
   Name = "H_B_Flip"; st.H_B_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Flip->Sumw2();
   Name = "H_C_Flip"; st.H_C_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Flip->Sumw2();
   Name = "H_D_Flip"; st.H_D_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Flip->Sumw2();
   Name = "H_E_Flip"; st.H_E_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_Flip->Sumw2();
   Name = "H_F_Flip"; st.H_F_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_Flip->Sumw2();
   Name = "H_G_Flip"; st.H_G_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_Flip->Sumw2();
   Name = "H_H_Flip"; st.H_H_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_Flip->Sumw2();

   Name = "H_A_Cen_Flip"; st.H_A_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen_Flip->Sumw2();
   Name = "H_B_Cen_Flip"; st.H_B_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen_Flip->Sumw2();
   Name = "H_C_Cen_Flip"; st.H_C_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen_Flip->Sumw2();
   Name = "H_D_Cen_Flip"; st.H_D_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen_Flip->Sumw2();
   Name = "H_E_Cen_Flip"; st.H_E_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_Cen_Flip->Sumw2();
   Name = "H_F_Cen_Flip"; st.H_F_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_Cen_Flip->Sumw2();
   Name = "H_G_Cen_Flip"; st.H_G_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_Cen_Flip->Sumw2();
   Name = "H_H_Cen_Flip"; st.H_H_Cen_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_Cen_Flip->Sumw2();

   Name = "H_A_For_Flip"; st.H_A_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For_Flip->Sumw2();
   Name = "H_B_For_Flip"; st.H_B_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For_Flip->Sumw2();
   Name = "H_C_For_Flip"; st.H_C_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For_Flip->Sumw2();
   Name = "H_D_For_Flip"; st.H_D_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For_Flip->Sumw2();
   Name = "H_E_For_Flip"; st.H_E_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_For_Flip->Sumw2();
   Name = "H_F_For_Flip"; st.H_F_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_For_Flip->Sumw2();
   Name = "H_G_For_Flip"; st.H_G_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_For_Flip->Sumw2();
   Name = "H_H_For_Flip"; st.H_H_For_Flip = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_For_Flip->Sumw2();

   Name = "Pred_I_Flip"; st.Pred_I_Flip = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 200,GlobalMinIm,dEdxM_UpLim); st.Pred_I_Flip->Sumw2();
   Name = "Pred_EtaB_Flip"; st.Pred_EtaB_Flip = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaB_Flip->Sumw2();
   Name = "Pred_EtaS_Flip"; st.Pred_EtaS_Flip = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaS_Flip->Sumw2();
   Name = "Pred_EtaS2_Flip"; st.Pred_EtaS2_Flip = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50,-3,3); st.Pred_EtaS2_Flip->Sumw2();
   Name = "Pred_EtaP_Flip"; st.Pred_EtaP_Flip = new TH3D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, 50, -3, 3, 200,GlobalMinPt,PtHistoUpperBound); st.Pred_EtaP_Flip->Sumw2();
   Name = "Pred_TOF_Flip"; st.Pred_TOF_Flip = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts,   200,GlobalMinTOF,5); st.Pred_TOF_Flip->Sumw2();

   Name = "H_D_DzSidebands"; st.H_D_DzSidebands = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, DzRegions, 0, DzRegions); st.H_D_DzSidebands->Sumw2();
   Name = "H_D_DzSidebands_DT"; st.H_D_DzSidebands_DT = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, DzRegions, 0, DzRegions); st.H_D_DzSidebands_DT->Sumw2();
   Name = "H_D_DzSidebands_CSC"; st.H_D_DzSidebands_CSC = new TH2D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts, DzRegions, 0, DzRegions); st.H_D_DzSidebands_CSC->Sumw2();
   /*
   for(int i=0; i<DzRegions; i++) {
     Name = "H_A_Syst_"+RegionNames[i]; st.H_A_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Syst[i]->Sumw2();
     Name = "H_B_Syst_"+RegionNames[i]; st.H_B_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Syst[i]->Sumw2();
     Name = "H_C_Syst_"+RegionNames[i]; st.H_C_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Syst[i]->Sumw2();
     Name = "H_D_Syst_"+RegionNames[i]; st.H_D_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Syst[i]->Sumw2();
     Name = "H_E_Syst_"+RegionNames[i]; st.H_E_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_Syst[i]->Sumw2();
     Name = "H_F_Syst_"+RegionNames[i]; st.H_F_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_Syst[i]->Sumw2();
     Name = "H_G_Syst_"+RegionNames[i]; st.H_G_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_Syst[i]->Sumw2();
     Name = "H_H_Syst_"+RegionNames[i]; st.H_H_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_Syst[i]->Sumw2();

     Name = "H_A_Cen_Syst_"+RegionNames[i]; st.H_A_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_Cen_Syst[i]->Sumw2();
     Name = "H_B_Cen_Syst_"+RegionNames[i]; st.H_B_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_Cen_Syst[i]->Sumw2();
     Name = "H_C_Cen_Syst_"+RegionNames[i]; st.H_C_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_Cen_Syst[i]->Sumw2();
     Name = "H_D_Cen_Syst_"+RegionNames[i]; st.H_D_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_Cen_Syst[i]->Sumw2();
     Name = "H_E_Cen_Syst_"+RegionNames[i]; st.H_E_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_Cen_Syst[i]->Sumw2();
     Name = "H_F_Cen_Syst_"+RegionNames[i]; st.H_F_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_Cen_Syst[i]->Sumw2();
     Name = "H_G_Cen_Syst_"+RegionNames[i]; st.H_G_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_Cen_Syst[i]->Sumw2();
     Name = "H_H_Cen_Syst_"+RegionNames[i]; st.H_H_Cen_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_Cen_Syst[i]->Sumw2();

     Name = "H_A_For_Syst_"+RegionNames[i]; st.H_A_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_A_For_Syst[i]->Sumw2();
     Name = "H_B_For_Syst_"+RegionNames[i]; st.H_B_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_B_For_Syst[i]->Sumw2();
     Name = "H_C_For_Syst_"+RegionNames[i]; st.H_C_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_C_For_Syst[i]->Sumw2();
     Name = "H_D_For_Syst_"+RegionNames[i]; st.H_D_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_D_For_Syst[i]->Sumw2();
     Name = "H_E_For_Syst_"+RegionNames[i]; st.H_E_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_E_For_Syst[i]->Sumw2();
     Name = "H_F_For_Syst_"+RegionNames[i]; st.H_F_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_F_For_Syst[i]->Sumw2();
     Name = "H_G_For_Syst_"+RegionNames[i]; st.H_G_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_G_For_Syst[i]->Sumw2();
     Name = "H_H_For_Syst_"+RegionNames[i]; st.H_H_For_Syst[i] = new TH1D(Name.c_str(), Name.c_str() ,NCuts,0,NCuts); st.H_H_For_Syst[i]->Sumw2();
   }

   Name = "H_DzCounts"; st.H_DzCounts = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts->Sumw2();
   Name = "H_DzCounts_DT"; st.H_DzCounts_DT = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts_DT->Sumw2();
   Name = "H_DzCounts_CSC"; st.H_DzCounts_CSC = new TH1D(Name.c_str(), Name.c_str() ,DzRegions,0,DzRegions); st.H_DzCounts_CSC->Sumw2();
   */

   Name = "CtrlPt_S1_Is"; st.CtrlPt_S1_Is = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxS_UpLim); st.CtrlPt_S1_Is->Sumw2();
   Name = "CtrlPt_S2_Is"; st.CtrlPt_S2_Is = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxS_UpLim); st.CtrlPt_S2_Is->Sumw2();
   Name = "CtrlPt_S3_Is"; st.CtrlPt_S3_Is = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxS_UpLim); st.CtrlPt_S3_Is->Sumw2();
   Name = "CtrlPt_S4_Is"; st.CtrlPt_S4_Is = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxS_UpLim); st.CtrlPt_S4_Is->Sumw2();

   Name = "CtrlPt_S1_Im"; st.CtrlPt_S1_Im = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxM_UpLim); st.CtrlPt_S1_Im->Sumw2();
   Name = "CtrlPt_S2_Im"; st.CtrlPt_S2_Im = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxM_UpLim); st.CtrlPt_S2_Im->Sumw2();
   Name = "CtrlPt_S3_Im"; st.CtrlPt_S3_Im = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxM_UpLim); st.CtrlPt_S3_Im->Sumw2();
   Name = "CtrlPt_S4_Im"; st.CtrlPt_S4_Im = new TH1D(Name.c_str(), Name.c_str(),200,0,dEdxM_UpLim); st.CtrlPt_S4_Im->Sumw2();

   Name = "CtrlIs_S1_TOF"; st.CtrlIs_S1_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIs_S1_TOF->Sumw2();
   Name = "CtrlIs_S2_TOF"; st.CtrlIs_S2_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIs_S2_TOF->Sumw2();
   Name = "CtrlIs_S3_TOF"; st.CtrlIs_S3_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIs_S3_TOF->Sumw2();
   Name = "CtrlIs_S4_TOF"; st.CtrlIs_S4_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIs_S4_TOF->Sumw2();

   Name = "CtrlIm_S1_TOF"; st.CtrlIm_S1_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIm_S1_TOF->Sumw2();
   Name = "CtrlIm_S2_TOF"; st.CtrlIm_S2_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIm_S2_TOF->Sumw2();
   Name = "CtrlIm_S3_TOF"; st.CtrlIm_S3_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIm_S3_TOF->Sumw2();
   Name = "CtrlIm_S4_TOF"; st.CtrlIm_S4_TOF = new TH1D(Name.c_str(), Name.c_str(),200,0,5); st.CtrlIm_S4_TOF->Sumw2();

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
   }

   st.Tree = new TTree("HscpCandidates", "HscpCandidates");
   st.Tree->SetDirectory(0);
   st.Tree->Branch("Run"     ,&st.Tree_Run       ,"Run/i");
   st.Tree->Branch("Event"   ,&st.Tree_Event     ,"Event/i");
   st.Tree->Branch("Hscp"    ,&st.Tree_Hscp      ,"Hscp/i");
   st.Tree->Branch("Pt"      ,&st.Tree_Pt        ,"Pt/F");
   st.Tree->Branch("I"       ,&st.Tree_I         ,"I/F");
   st.Tree->Branch("TOF"     ,&st.Tree_TOF       ,"TOF/F");
   st.Tree->Branch("Mass"    ,&st.Tree_Mass      ,"Mass/F");

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
   st.V3D               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/V3D");
   st.Chi2              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Chi2");
   st.Qual              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Qual");
   st.TNOH              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TNOH");
   st.TNOM              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TNOM");
   st.nDof              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/nDof");
   st.Pterr             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Pterr");
   st.TIsol             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TIsol");
   st.EIsol             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/EIsol");
   st.MPt               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MPt");
   st.MI                = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MI");
   st.MTOF              = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MTOF");
   st.Pt                = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Pt");
   st.I                 = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/I");
   st.TOF               = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/TOF");
   st.HSCPE             = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE");

   st.HSCPE_SystP       = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystP");
   st.HSCPE_SystI       = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystI");
   st.HSCPE_SystM       = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystM");
   st.HSCPE_SystT       = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystT");
   st.HSCPE_SystPU      = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/HSCPE_SystPU");

   st.Mass              = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Mass");
   st.MassTOF           = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassTOF");
   st.MassComb          = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassComb");
   st.MaxEventMass      = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MaxEventMass");

   st.Mass_SystP        = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Mass_SystP");
   st.MassTOF_SystP     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassTOF_SystP");
   st.MassComb_SystP    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassComb_SystP");
   st.MaxEventMass_SystP= (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MaxEventMass_SystP");

   st.Mass_SystI        = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Mass_SystI");
   st.MassTOF_SystI     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassTOF_SystI");
   st.MassComb_SystI    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassComb_SystI");
   st.MaxEventMass_SystI    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MaxEventMass_SystI");

   st.Mass_SystT        = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Mass_SystT");
   st.MassTOF_SystT     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassTOF_SystT");
   st.MassComb_SystT    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassComb_SystT");
   st.MaxEventMass_SystT    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MaxEventMass_SystT");

   st.Mass_SystPU        = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Mass_SystPU");
   st.MassTOF_SystPU     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassTOF_SystPU");
   st.MassComb_SystPU    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MassComb_SystPU");
   st.MaxEventMass_SystPU= (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/MaxEventMass_SystPU");

   st.Beta_Gen          = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Gen");
   st.Beta_GenCharged   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_GenCharged");
   st.Beta_Triggered    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Triggered");
   st.Beta_Matched      = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_Matched");
   st.Beta_PreselectedA = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_PreselectedA");
   st.Beta_PreselectedB = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_PreselectedB");
   st.Beta_PreselectedC = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_PreselectedC");
   st.Beta_SelectedP    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_SelectedP");
   st.Beta_SelectedI    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_SelectedI");
   st.Beta_SelectedT    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/Beta_SelectedT");

   st.BS_V3D    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_V3D");
   st.BS_Chi2   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Chi2");
   st.BS_Qual   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Qual");
   st.BS_TNOH   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOH");
   st.BS_TNOHFraction   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOHFraction");
   st.BS_Eta    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Eta");
   st.BS_TNOM   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TNOM");
   st.BS_nDof   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_nDof");
   st.BS_Pterr  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtErr");
   st.BS_MPt    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_MPt");
   st.BS_MIm    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_MIm");
   st.BS_MIs    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_MIs");
   st.BS_MTOF   = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_MTOF");
   st.BS_TIsol  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TIsol");
   st.BS_EIsol  = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EIsol");
   st.BS_P      = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_P");
   st.AS_P      = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_P");
   st.BS_Pt     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Pt");
   st.AS_Pt     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Pt");
   st.BS_Im     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Im");
   st.AS_Im     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Im");
   st.BS_Is     = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_Is");
   st.AS_Is     = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Is");
   st.BS_TOF    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF");
   st.BS_TOF_DT    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_DT");
   st.BS_TOF_CSC    = (TH1F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOF_CSC");
   st.AS_TOF    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_TOF");
   st.BS_EtaIs  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaIs");
   //st.AS_EtaIs  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_EtaIs");
   st.BS_EtaIm  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaIm");
   //st.AS_EtaIm  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_EtaIm");
   st.BS_EtaP   = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaP");
   //st.AS_EtaP   = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_EtaP");
   st.BS_EtaPt  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaPt");
   //st.AS_EtaPt  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_EtaPt");
   st.BS_EtaTOF  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_EtaTOF");
   //st.AS_EtaTOF  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_EtaTOF");
   st.BS_PIs    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PIs");
   st.AS_PIs    = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_PIs");
   st.BS_PIm    = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PIm");
   st.AS_PIm    = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_PIm");
   st.BS_PtIs   = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtIs");
   st.AS_PtIs   = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_PtIs");
   st.BS_PtIm   = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_PtIm");
   st.AS_PtIm   = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_PtIm");
   st.BS_TOFIs  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOFIs");
   st.AS_TOFIs  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_TOFIs");
   st.BS_TOFIm  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/BS_TOFIm");
   st.AS_TOFIm  = (TH3F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_TOFIm");

   st.AS_Eta_RegionA  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionA");
   st.AS_Eta_RegionB  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionB");
   st.AS_Eta_RegionC  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionC");
   st.AS_Eta_RegionD  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionD");
   st.AS_Eta_RegionE  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionE");
   st.AS_Eta_RegionF  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionF");
   st.AS_Eta_RegionG  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionG");
   st.AS_Eta_RegionH  = (TH2F*)GetObjectFromPath(st.Directory, InputFile,  BaseName + "/AS_Eta_RegionH");
   HistoFile->cd();
}

void stPlots_Clear(stPlots& st, bool WriteFirst=false)
{
   if(WriteFirst){
      st.Tree->SetDirectory(st.Directory);
      st.Directory->Write();
   }
   delete st.Directory;
}


void stPlots_FillTree(stPlots& st, unsigned int Run, unsigned int Event, unsigned int Hscp, double Pt, double I, double TOF, double Mass, int MaxEntry=20000){
   if(MaxEntry>0 && st.Tree->GetEntries()>=MaxEntry)return;
   st.Tree_Run   = Run;
   st.Tree_Event = Event;
   st.Tree_Hscp  = Hscp;
   st.Tree_Pt    = Pt;
   st.Tree_I     = I;
   st.Tree_TOF   = TOF;
   st.Tree_Mass  = Mass;
   st.Tree->Fill();
}


void stPlots_Dump(stPlots& st, FILE* pFile, int CutIndex){

   fprintf(pFile,"---------- %10s ----------\n",st.Name.c_str());
   fprintf(pFile,"#Events                       = %4.2E\n",st.TotalE->GetBinContent(1       ));
   fprintf(pFile,"#Triggered Events             = %4.2E Eff=%4.3E\n",st.TotalTE->GetBinContent(1     ),st.TotalTE->GetBinContent(1      )/st.TotalE->GetBinContent(1       ));
   fprintf(pFile,"#Tracks                       = %4.2E\n",st.Total->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing TNOH   cuts   = %4.2E Eff=%4.3E\n",st.TNOH ->GetBinContent(1       ), st.TNOH ->GetBinContent(1       ) /st.Total->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing TNOM   cuts   = %4.2E Eff=%4.3E\n",st.TNOM ->GetBinContent(1       ), st.TNOM ->GetBinContent(1       ) /st.TNOH ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing nDof   cuts   = %4.2E Eff=%4.3E\n",st.nDof ->GetBinContent(1       ), st.nDof ->GetBinContent(1       ) /st.TNOM ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Qual   cuts   = %4.2E Eff=%4.3E\n",st.Qual ->GetBinContent(1       ), st.Qual ->GetBinContent(1       ) /st.nDof ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Chi2   cuts   = %4.2E Eff=%4.3E\n",st.Chi2 ->GetBinContent(1       ), st.Chi2 ->GetBinContent(1       ) /st.Qual ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Min Pt cuts   = %4.2E Eff=%4.3E\n",st.MPt  ->GetBinContent(1       ), st.MPt  ->GetBinContent(1       ) /st.Chi2 ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Min I  cuts   = %4.2E Eff=%4.3E\n",st.MI   ->GetBinContent(1       ), st.MI   ->GetBinContent(1       ) /st.MPt  ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Min TOFcuts   = %4.2E Eff=%4.3E\n",st.MTOF ->GetBinContent(1       ), st.MTOF ->GetBinContent(1       ) /st.MI   ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing V3D    cuts   = %4.2E Eff=%4.3E\n",st.V3D  ->GetBinContent(1       ), st.V3D  ->GetBinContent(1       ) /st.MI   ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing TIsol  cuts   = %4.2E Eff=%4.3E\n",st.TIsol->GetBinContent(1       ), st.TIsol->GetBinContent(1       ) /st.V3D  ->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing EIsol  cuts   = %4.2E Eff=%4.3E\n",st.EIsol->GetBinContent(1       ), st.EIsol->GetBinContent(1       ) /st.TIsol->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing PtErr  cuts   = %4.2E Eff=%4.3E\n",st.Pterr->GetBinContent(1       ), st.Pterr->GetBinContent(1       ) /st.EIsol->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Basic  cuts   = %4.2E Eff=%4.3E\n",st.Pterr->GetBinContent(1       ), st.Pterr->GetBinContent(1       ) /st.Total->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing Pt     cuts   = %4.2E Eff=%4.3E\n",st.Pt   ->GetBinContent(CutIndex+1), st.Pt   ->GetBinContent(CutIndex+1) /st.Pterr->GetBinContent(1       ));
   fprintf(pFile,"#Tracks passing I      cuts   = %4.2E Eff=%4.3E\n",st.I    ->GetBinContent(CutIndex+1), st.I    ->GetBinContent(CutIndex+1) /st.Pt   ->GetBinContent(CutIndex+1));
   fprintf(pFile,"#Tracks passing TOF    cuts   = %4.2E Eff=%4.3E\n",st.TOF  ->GetBinContent(CutIndex+1), st.TOF  ->GetBinContent(CutIndex+1) /st.I    ->GetBinContent(CutIndex+1));
   fprintf(pFile,"#Tracks passing selection     = %4.2E Eff=%4.3E\n",st.TOF  ->GetBinContent(CutIndex+1), st.TOF  ->GetBinContent(CutIndex+1) /st.Total->GetBinContent(1       ));   
   fprintf(pFile,"--------------------\n");
   fprintf(pFile,"HSCP Detection Efficiency Before Trigger                           Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(2*st.TotalE ->GetBinContent(1       )));
   fprintf(pFile,"HSCP Detection Efficiency After  Trigger                           Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(2*st.TotalTE->GetBinContent(1       )));
   fprintf(pFile,"#HSCPTrack per HSCPEvent (with at least one HSCPTrack)             Eff=%4.3E\n",st.TOF->GetBinContent(CutIndex+1) /(  st.HSCPE  ->GetBinContent(CutIndex+1)));
   fprintf(pFile,"--------------------\n");
}


void stPlots_Draw(stPlots& st, std::string SavePath, std::string LegendTitle, unsigned int CutIndex)
{
   TObject** Histos = new TObject*[10];
   std::vector<std::string> legend;
   TCanvas* c1;

   char CutIndexStr[255];sprintf(CutIndexStr,"_%03i",CutIndex);


   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_EtaIs;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaIs_BS", true);
   delete c1;

//   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
//   st.AS_EtaIs->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);   
//   Histos[0] = (TH1*)st.AS_EtaIs->Project3D("zy"); legend.push_back("After Cut");
//   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", dEdxS_Legend.c_str(), 0,0, 0,0, false);
//   c1->SetLogz(true);
//   DrawPreliminary(IntegratedLuminosity);
//   SaveCanvas(c1,SavePath,std::string("EtaIs_AS")+CutIndexStr, true);
//   delete Histos[0];
//   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_EtaIm;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", dEdxM_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaIm_BS", true);
   delete c1;

//   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
//   st.AS_EtaIm->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
//   Histos[0] = (TH1*)st.AS_EtaIm->Project3D("zy");legend.push_back("After Cut");
//   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", dEdxM_Legend.c_str(), 0,0, 0,0, false);
//   c1->SetLogz(true);
//   DrawPreliminary(IntegratedLuminosity);
//   SaveCanvas(c1,SavePath,std::string("EtaIm_AS")+CutIndexStr, true);
//   delete Histos[0];
//   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_EtaP;                  legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p (GeV/c)", 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaP_BS", true);
   delete c1;

//   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
//   st.AS_EtaP->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
//   Histos[0] = (TH1*)st.AS_EtaP->Project3D("zy"); legend.push_back("After Cut");
//   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p (GeV/c)", 0,0, 0,0, false);
//   c1->SetLogz(true);
//   DrawPreliminary(IntegratedLuminosity);
//   SaveCanvas(c1,SavePath,std::string("EtaP_AS")+CutIndexStr, true);
//   delete Histos[0];
//   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_EtaPt;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p_{T} (GeV/c)", 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaPt_BS", true);
   delete c1;

//   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
//   st.AS_EtaPt->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
//   Histos[0] = (TH1*)st.AS_EtaPt->Project3D("zy");legend.push_back("After Cut");
//   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "p_{T} (GeV/c)", 0,0, 0,0, false);
//   c1->SetLogz(true);
//   DrawPreliminary(IntegratedLuminosity);
//   SaveCanvas(c1,SavePath,std::string("EtaPt_AS")+CutIndexStr, true);
//   delete Histos[0];
//   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_EtaTOF;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "1/#beta", 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"EtaTOF_BS", true);
   delete c1;

//   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
//   st.AS_EtaTOF->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
//   Histos[0] = (TH1*)st.AS_EtaTOF->Project3D("zy");legend.push_back("After Cut");
//   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", "1/#beta", 0,0, 0,0, false);
//   c1->SetLogz(true);
//   DrawPreliminary(IntegratedLuminosity);
//   SaveCanvas(c1,SavePath,std::string("EtaTOF_AS")+CutIndexStr, true);
//   delete Histos[0];
//   delete c1;


   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_PIs;                   legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV/c)", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PIs_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_PIm;                   legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV/c)", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PIm_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_PtIs;                  legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV/c)", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PtIs_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_PtIm;                  legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV/c)", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"PtIm_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_PIs->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_PIs->Project3D("zy");  legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV/c)", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("PIs_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_PIm->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_PIm->Project3D("zy");  legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV/c)", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("PIm_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_PtIs->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_PtIs->Project3D("zy"); legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV/c)", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("PtIs_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_PtIm->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_PtIm->Project3D("zy"); legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV/c)", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("PtIm_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;


   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_TOFIs;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOFIs_BS", true);
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)st.BS_TOFIm;                 legend.push_back("Before Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOFIm_BS", true);
   delete c1;


   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_TOFIs->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_TOFIs->Project3D("zy");legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("TOFIs_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   st.AS_TOFIm->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   Histos[0] = (TH1*)st.AS_TOFIm->Project3D("zy");legend.push_back("After Cut");
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxM_Legend.c_str(), 0,0, 0,15, false);
   c1->SetLogz(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("TOFIm_AS")+CutIndexStr, true);
   delete Histos[0];
   delete c1;


   TH1** Histos1D = new TH1*[10];
   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos1D[0] = (TH1*)st.AS_Eta_RegionA->ProjectionY((st.Name+"A").c_str(),CutIndex+1,CutIndex+1); legend.push_back("A");  
   if(Histos1D[0]->Integral()>0) Histos1D[0]->Scale(1.0/Histos1D[0]->Integral());
   Histos1D[1] = (TH1*)st.AS_Eta_RegionB->ProjectionY((st.Name+"B").c_str(),CutIndex+1,CutIndex+1); legend.push_back("B");
   if(Histos1D[1]->Integral()>0) Histos1D[1]->Scale(1.0/Histos1D[1]->Integral());
   //Histos1D[2] = (TH1*)st.AS_Eta_RegionC->ProjectionY((st.Name+"C").c_str(),CutIndex+1,CutIndex+1); legend.push_back("C");
   //if(Histos1D[2]->Integral()>0) Histos1D[2]->Scale(1.0/Histos1D[2]->Integral());
   //Histos1D[3] = (TH1*)st.AS_Eta_RegionD->ProjectionY((st.Name+"D").c_str(),CutIndex+1,CutIndex+1); legend.push_back("D");
   //if(Histos1D[3]->Integral()>0) Histos1D[3]->Scale(1.0/Histos1D[3]->Integral());
   Histos1D[2] = (TH1*)st.AS_Eta_RegionE->ProjectionY((st.Name+"E").c_str(),CutIndex+1,CutIndex+1); legend.push_back("E");
   if(Histos1D[2]->Integral()>0) Histos1D[2]->Scale(1.0/Histos1D[2]->Integral());
   Histos1D[3] = (TH1*)st.AS_Eta_RegionF->ProjectionY((st.Name+"F").c_str(),CutIndex+1,CutIndex+1); legend.push_back("F");
   if(Histos1D[3]->Integral()>0) Histos1D[3]->Scale(1.0/Histos1D[3]->Integral());
   //Histos1D[6] = (TH1*)st.AS_Eta_RegionG->ProjectionY((st.Name+"G").c_str(),CutIndex+1,CutIndex+1); legend.push_back("G");
   //if(Histos1D[6]->Integral()>0) Histos1D[6]->Scale(1.0/Histos1D[6]->Integral());
   //Histos1D[7] = (TH1*)st.AS_Eta_RegionH->ProjectionY((st.Name+"H").c_str(),CutIndex+1,CutIndex+1); legend.push_back("H");
   //if(Histos1D[7]->Integral()>0) Histos1D[7]->Scale(1.0/Histos1D[7]->Integral());
   DrawSuperposedHistos((TH1**)Histos1D, legend, "E1",  "p_{T} (GeV/c)", "arbitrary units", 0, 0, 0, 0);
   DrawLegend((TObject**)Histos1D,legend,LegendTitle,"P");
   c1->SetLogy(false);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("EtaRegions_AS")+CutIndexStr);
   //for(unsigned int i=0;i<8;i++){delete Histos1D[i];}
   delete c1;
}

void stPlots_DrawComparison(std::string SavePath, std::string LegendTitle, unsigned int CutIndex, stPlots* st1, stPlots* st2=NULL, stPlots* st3=NULL, stPlots* st4=NULL, stPlots* st5=NULL, stPlots* st6=NULL, stPlots* st7=NULL)
{ 
   char CutIndexStr[255];sprintf(CutIndexStr,"_%03i",CutIndex);

   bool IsTkOnly = (SavePath.find("Type0",0)<std::string::npos);

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
     if(st[i]->Name=="MCTr"){lg.push_back("MC - SM");}
     else if(st[i]->Name=="Data"){lg.push_back("Observed");}
     else if(Index==-1){lg.push_back(st[i]->Name);}else{lg.push_back(signals[Index].Legend);}
  }
   
   TH1** Histos = new TH1*[10];
   std::vector<std::string> legend;
   TCanvas* c1;

   for(unsigned int i=0;i<st.size();i++){
//      if(st[i]->Name=="Data")continue;
      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
      Histos[0] = (TH1*)st[i]->Beta_Gen;                                                  legend.push_back("Gen");
//      Histos[1] = (TH1*)st[i]->Beta_GenCharged;                                           legend.push_back("Charged Gen");
      Histos[1] = (TH1*)st[i]->Beta_Triggered;                                            legend.push_back("Triggered");
      DrawSuperposedHistos((TH1**)Histos, legend,"HIST E1",  "#beta", "# HSCP", 0,0, 0,0);
      DrawLegend((TObject**)Histos,legend,"","P", 0.36, 0.92, 0.20, 0.04);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,st[i]->Name + "_GenBeta", true);
      delete c1;
   }

   for(unsigned int i=0;i<st.size();i++){
//      if(st[i]->Name=="Data")continue;
      c1 = new TCanvas("c1","c1,",600,600);                                               legend.clear();
//    Histos[0] = (TH1*)st[i]->Beta_Gen;                                                  legend.push_back("Gen");
//    Histos[1] = (TH1*)st[i]->Beta_GenCharged;                                           legend.push_back("Charged Gen");
      Histos[0] = (TH1*)st[i]->Beta_Triggered;                                            legend.push_back("Triggered");
      Histos[1] = (TH1*)st[i]->Beta_Matched;                                              legend.push_back("Reconstructed");
//    Histos[0] = (TH1*)st[i]->Beta_PreselectedA;                                         legend.push_back("PreselectedA");
//    Histos[0] = (TH1*)st[i]->Beta_PreselectedB;                                         legend.push_back("PreselectedB");
      Histos[2] = (TH1*)st[i]->Beta_PreselectedC;                                         legend.push_back("Preselected");
      Histos[3] = (TH1*)st[i]->Beta_SelectedP->ProjectionY("A",CutIndex+1,CutIndex+1);    legend.push_back("p_{T}>Cut");
      Histos[4] = (TH1*)st[i]->Beta_SelectedI->ProjectionY("B",CutIndex+1,CutIndex+1);    legend.push_back("I  >Cut");
      if(!IsTkOnly){Histos[5] = (TH1*)st[i]->Beta_SelectedT->ProjectionY("C",CutIndex+1,CutIndex+1);    legend.push_back("ToF>Cut");}
      DrawSuperposedHistos((TH1**)Histos, legend,"HIST E1",  "#beta", "# HSCP", 0,0, 0,0);
      DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.36, 0.92, 0.20, 0.025);
      c1->SetLogy(true);
      DrawPreliminary(IntegratedLuminosity);
      SaveCanvas(c1,SavePath,st[i]->Name + "_Beta");
      delete Histos[3]; delete Histos[4]; 
      if(!IsTkOnly)delete Histos[5];
      delete c1;
   }



   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_V3D->Clone();      legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral());   }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "V3D (cm)", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"V3D_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Chi2->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#chi^{2}/ndof", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Chi2_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Qual->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "quality", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Quality_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TNOH->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#NOH", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOH_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_TNOHFraction->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Fraction of hits", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P",0.49);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOHFraction_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
     Histos[i] = (TH1*)st[i]->BS_Eta->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#eta", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Eta_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TNOM->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "#NOM", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"NOM_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_nDof->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "TOF_{nDof}", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"nDof_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Pterr->Clone();       legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} Err / p_{T}", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pterr_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_MPt->Clone();         legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV/c)", "arbitrary units", 0,1250, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"MPt_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_MIs->Clone();         legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxS_Legend.c_str(), "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.79, 0.19);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"MIs_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_MIm->Clone();         legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxM_Legend.c_str(), "arbitrary units", 0,20, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"MIm_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_MTOF->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", "arbitrary units", -2,5, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"MTOF_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TIsol->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); } 
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Isolation: Track SumPt (GeV/c)", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"IsolT_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_EIsol->Clone();        legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Isolation: (Ecal + Hcal) Energy / p", "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"IsolE_BS", true);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Is; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxS_Legend.c_str(), tmp, 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxS_Legend.c_str(), tmp, 0,1, 0,0, false, true);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Is_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Im; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxM_Legend.c_str(), "arbitrary units", 0,20, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Im_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)(st[i]->AS_Is->ProjectionY((st[i]->Name+"AA").c_str(),CutIndex+1,CutIndex+1)); legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxS_Legend.c_str(), "arbitrary units", 0,0, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.79, 0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("Is_AS")+CutIndexStr);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->AS_Im->ProjectionY((st[i]->Name+"BB").c_str(),CutIndex+1,CutIndex+1); legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  dEdxM_Legend.c_str(), "arbitrary units", 0,20, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("Im_AS")+CutIndexStr);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_Pt; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   sprintf(tmp,"Fraction of tracks/%2.0f GeV/#font[12]{c}",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV/#font[12]{c})", tmp, 0,1250, 0.000000001, 1.2);
   if(IsTkOnly) DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.45, 0.42, 0.26, 0.05);
   else DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.51, 0.39, 0.33, 0.05);
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV/#font[12]{c})", tmp, 0,1250, 0.000000001, 1.2, false, true);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"Pt_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->AS_Pt->ProjectionY((st[i]->Name+"CC").c_str(),CutIndex+1,CutIndex+1); legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV/c)", "arbitrary units", 0,1250, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("Pt_AS")+CutIndexStr);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   //char tmp[2048];
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_DT; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P", 0.85);//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_DT_BS");
   delete c1;

   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->BS_TOF_CSC; legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   sprintf(tmp,"Fraction of tracks/%0.2f",Histos[0]->GetBinWidth(1));
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", tmp, 0.5, 1.5, 0, 0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//,0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,"TOF_CSC_BS");
   delete c1;
   
   c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   for(unsigned int i=0;i<st.size();i++){
   Histos[i] = (TH1*)st[i]->AS_TOF->ProjectionY((st[i]->Name+"DD").c_str(),CutIndex+1,CutIndex+1); legend.push_back(lg[i]);  if(Histos[i]->Integral()>0) Histos[i]->Scale(1.0/Histos[i]->Integral()); }
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "1/#beta", "arbitrary units", 1,4, 0,0);
   DrawLegend((TObject**)Histos,legend,LegendTitle,"P");//, 0.35);
   c1->SetLogy(true);
   DrawPreliminary(IntegratedLuminosity);
   SaveCanvas(c1,SavePath,std::string("TOF_AS")+CutIndexStr);
   for(unsigned int i=0;i<st.size();i++){delete Histos[i];}
   delete c1;

}
