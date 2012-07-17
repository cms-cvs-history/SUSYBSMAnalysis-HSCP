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

using namespace std;

#include "Analysis_Global.h"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"
#include "Analysis_Samples.h"

/////////////////////////// CODE PARAMETERS /////////////////////////////

void Analysis_Step4(string MODE_="COMPILE", string InputPattern="")
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

   printf("Step4: Doing final computations\n");

   string Input;
   if(MODE=="ANALYSE_DATA") Input     = InputPattern + "Histos_Data.root";
   else Input     = InputPattern + "Histos_MC.root";
   string SavePath  = InputPattern;
   MakeDirectories(SavePath);

   TFile* InputFile = new TFile(Input.c_str(), "UPDATE");
   TH1D*  HCuts_Pt       = (TH1D*)GetObjectFromPath(InputFile, "HCuts_Pt");
   TH1D*  HCuts_I        = (TH1D*)GetObjectFromPath(InputFile, "HCuts_I");
   TH1D*  HCuts_TOF      = (TH1D*)GetObjectFromPath(InputFile, "HCuts_TOF");

   TH3D*  Pred_EtaP            = (TH3D*)GetObjectFromPath(InputFile, "Data/Pred_EtaP");
   TH2D*  Pred_I            = (TH2D*)GetObjectFromPath(InputFile, "Data/Pred_I");
   TH2D*  Pred_TOF            = (TH2D*)GetObjectFromPath(InputFile, "Data/Pred_TOF");
   TH2D*  Pred_EtaB            = (TH2D*)GetObjectFromPath(InputFile, "Data/Pred_EtaB");
   TH2D*  Pred_EtaS            = (TH2D*)GetObjectFromPath(InputFile, "Data/Pred_EtaS");
   TH2D*  Pred_EtaS2            = (TH2D*)GetObjectFromPath(InputFile, "Data/Pred_EtaS2");

   TH1D*  H_A            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A");
   TH1D*  H_B            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B");
   TH1D*  H_C            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C");
   TH1D*  H_D            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D");
   TH1D*  H_E            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_E");
   TH1D*  H_F            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_F");
   TH1D*  H_G            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_G");
   TH1D*  H_H            = (TH1D*)GetObjectFromPath(InputFile, "Data/H_H");
   TH1D*  H_P = new TH1D("H_P" ,"H_P" ,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX());

   TH1D*  H_A_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_Cen");
   TH1D*  H_B_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_Cen");
   TH1D*  H_C_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_Cen");
   TH1D*  H_D_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_Cen");
   TH1D*  H_E_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_E_Cen");
   TH1D*  H_F_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_F_Cen");
   TH1D*  H_G_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_G_Cen");
   TH1D*  H_H_Cen     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_H_Cen");
   TH1D*  H_P_Cen = new TH1D("H_P_Cen" ,"H_P_Cen" ,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX());

   TH1D*  H_A_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_A_For");
   TH1D*  H_B_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_B_For");
   TH1D*  H_C_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_C_For");
   TH1D*  H_D_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_D_For");
   TH1D*  H_E_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_E_For");
   TH1D*  H_F_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_F_For");
   TH1D*  H_G_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_G_For");
   TH1D*  H_H_For     = (TH1D*)GetObjectFromPath(InputFile, "Data/H_H_For");
   TH1D*  H_P_For = new TH1D("H_P_For" ,"H_P_For" ,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX());

   char Name   [1024];
   sprintf(Name,"Pred_Mass");
   TH2D* Pred_Mass = new TH2D(Name,Name,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX(),MassNBins,0,MassHistoUpperBound);
   Pred_Mass->Sumw2();

   sprintf(Name,"Pred_MassTOF");
   TH2D* Pred_MassTOF = new TH2D(Name,Name,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX(), MassNBins,0,MassHistoUpperBound);
   Pred_MassTOF->Sumw2();

   sprintf(Name,"Pred_MassComb");
   TH2D* Pred_MassComb = new TH2D(Name,Name,HCuts_Pt->GetNbinsX(),0,HCuts_Pt->GetNbinsX(),MassNBins,0,MassHistoUpperBound);
   Pred_MassComb->Sumw2();

   //////////////////////////////////////////////////      MAKING THE PREDICTION
   for(int CutIndex=0;CutIndex<100;CutIndex++){
     //for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){

      const double& A=H_A->GetBinContent(CutIndex+1);
      const double& B=H_B->GetBinContent(CutIndex+1);
      const double& C=H_C->GetBinContent(CutIndex+1);
      const double& D=H_D->GetBinContent(CutIndex+1);
      const double& E=H_E->GetBinContent(CutIndex+1);
      const double& F=H_F->GetBinContent(CutIndex+1);
      const double& G=H_G->GetBinContent(CutIndex+1);
      const double& H=H_H->GetBinContent(CutIndex+1);

      const double& A_Cen=H_B_Cen->GetBinContent(CutIndex+1);
      const double& B_Cen=H_D_Cen->GetBinContent(CutIndex+1);
      const double& C_Cen=H_F_Cen->GetBinContent(CutIndex+1);
      const double& D_Cen=H_H_Cen->GetBinContent(CutIndex+1);
      const double& E_Cen=H_B_Cen->GetBinContent(CutIndex+1);
      const double& F_Cen=H_D_Cen->GetBinContent(CutIndex+1);
      const double& G_Cen=H_F_Cen->GetBinContent(CutIndex+1);
      const double& H_Cen=H_H_Cen->GetBinContent(CutIndex+1);

      const double& A_For=H_A_For->GetBinContent(CutIndex+1);
      const double& B_For=H_B_For->GetBinContent(CutIndex+1);
      const double& C_For=H_C_For->GetBinContent(CutIndex+1);
      const double& D_For=H_D_For->GetBinContent(CutIndex+1);
      const double& E_For=H_E_For->GetBinContent(CutIndex+1);
      const double& F_For=H_F_For->GetBinContent(CutIndex+1);
      const double& G_For=H_G_For->GetBinContent(CutIndex+1);
      const double& H_For=H_H_For->GetBinContent(CutIndex+1);

      double P=0;
      double Perr=0;
      double P_Cen=0;
      double Perr_Cen=0;
      double P_For=0;
      double Perr_For=0;

      printf("%4i --> Pt>%7.2f  I>%6.2f  TOF>%+5.2f --> A=%6.2E B=%6.E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E\n",CutIndex,HCuts_Pt->GetBinContent(CutIndex+1), HCuts_I->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1), A, B, C, D, E, F, G, H );

      if(E>0){
         P    = (A*F*G)/(E*E);
         Perr = sqrt( ((pow(F*G,2)* A + pow(A*G,2)*F + pow(A*F,2)*G)/pow(E,4)) + (pow((2*A*F*G)/pow(E,3),2)*E));
      }else if(F>0){
	P_Cen=((H_Cen*B_Cen)/F_Cen);
        Perr_Cen = sqrt( (pow(B_Cen/F_Cen,2)*H_Cen) + (pow(H_Cen/F_Cen,2)*B_Cen) + (pow((B_Cen*(H_Cen)/(F_Cen*F_Cen)),2)*F_Cen) );
	P_For=((H_For*B_For)/F_For);
        Perr_For = sqrt( (pow(B_For/F_For,2)*H_For) + (pow(H_For/F_For,2)*B_For) + (pow((B_For*(H_For)/(F_For*F_For)),2)*F_For) );
	P    = P_Cen + P_For;
	Perr = sqrt(Perr_Cen*Perr_Cen + Perr_For*Perr_For);
      }else if(A>0){
         P    = ((C*B)/A);
         Perr = sqrt( (pow(B/A,2)*C) + (pow(C/A,2)*B) + (pow((B*(C)/(A*A)),2)*A) );
      }

      H_P->SetBinContent(CutIndex+1,P);
      H_P->SetBinError  (CutIndex+1,Perr);
      if(P==0 || isnan((float)P))continue; //Skip this CutIndex --> No Prediction possible

      printf("%4i --> Pt>%7.2f  I>%6.2f  TOF>%+5.2f --> D=%6.2E vs Pred = %6.2E +- %6.2E (%6.2E%%)\n", CutIndex,HCuts_Pt->GetBinContent(CutIndex+1), HCuts_I->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1),D, P,  Perr, 100.0*Perr/P );

      TH1D* Pred_EtaB_Proj = Pred_EtaB->ProjectionY("ProjEtaB",CutIndex+1,CutIndex+1);  // Pred_EtaB_Proj->Scale(1.0/Pred_EtaB_Proj->Integral());
      TH1D* Pred_EtaS_Proj = Pred_EtaS->ProjectionY("ProjEtaS",CutIndex+1,CutIndex+1); //  Pred_EtaS_Proj->Scale(1.0/Pred_EtaS_Proj->Integral());
      TH1D* Pred_EtaS2_Proj = Pred_EtaS2->ProjectionY("ProjEtaS2",CutIndex+1,CutIndex+1);//   Pred_EtaS2_Proj->Scale(1.0/Pred_EtaS2_Proj->Integral())
      TH1D* Pred_EtaB_Proj_PE  = (TH1D*)Pred_EtaB_Proj->Clone("Pred_EtaB_Proj_PE");   Pred_EtaB_Proj_PE->Reset();
      TH1D* Pred_EtaS_Proj_PE  = (TH1D*)Pred_EtaS_Proj->Clone("Pred_EtaS_Proj_PE");   Pred_EtaS_Proj_PE->Reset();
      TH1D* Pred_EtaS2_Proj_PE = (TH1D*)Pred_EtaS2_Proj->Clone("Pred_EtaS2_Proj_PE"); Pred_EtaS2_Proj_PE->Reset();
      Pred_EtaP->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
      TH2D* Pred_EtaPWeighted = (TH2D*)Pred_EtaP->Project3D("zy");
      TH2D* Pred_EtaPWeighted_PE = (TH2D*)Pred_EtaPWeighted->Clone("Pred_EtaPWeightedPE");   Pred_EtaPWeighted_PE->Reset();

      /*
      for(int x=0;x<=Pred_EtaPWeighted->GetXaxis()->GetNbins();x++){
         double WeightP = 0.0;
         if(Pred_EtaB_Proj->GetBinContent(x)>0){
            WeightP = Pred_EtaS_Proj->GetBinContent(x)/Pred_EtaB_Proj->GetBinContent(x);
            if(TypeMode==2)WeightP*= Pred_EtaS2_Proj->GetBinContent(x)/Pred_EtaB_Proj->GetBinContent(x);
         }

         for(int y=0;y<=Pred_EtaPWeighted->GetYaxis()->GetNbins();y++){
            Pred_EtaPWeighted->SetBinContent(x,y,Pred_EtaPWeighted->GetBinContent(x,y)*WeightP);
         }
      }
*/

//      TH1D* Pred_P_Proj = Pred_EtaPWeighted->ProjectionY("ProjP");
      TH1D* Pred_I_Proj = Pred_I->ProjectionY("ProjI",CutIndex+1,CutIndex+1);
      TH1D* Pred_T_Proj = Pred_TOF->ProjectionY("ProjT",CutIndex+1,CutIndex+1);
      TH1D* Pred_I_ProjPE = (TH1D*) Pred_I_Proj->Clone("Pred_I_ProjPE"); Pred_I_ProjPE->Reset();
      TH1D* Pred_T_ProjPE = (TH1D*) Pred_T_Proj->Clone("Pred_T_ProjPE"); Pred_T_ProjPE->Reset();


//      TH1D* Pred_P_PDF = GetPDF(Pred_P_Proj);
//      TH1D* Pred_I_PDF = GetPDF(Pred_I_Proj);
//      TH1D* Pred_T_PDF = GetPDF(Pred_T_Proj);

      TH2D* Pred_Prof_Mass     =  new TH2D("Pred_Prof_Mass"    ,"Pred_Prof_Mass"    ,MassNBins,0,MassHistoUpperBound, 100, 0, 100); 
      TH2D* Pred_Prof_MassTOF  =  new TH2D("Pred_Prof_MassTOF" ,"Pred_Prof_MassTOF" ,MassNBins,0,MassHistoUpperBound, 100, 0, 100);  
      TH2D* Pred_Prof_MassComb =  new TH2D("Pred_Prof_MassComb","Pred_Prof_MassComb",MassNBins,0,MassHistoUpperBound, 100, 0, 100);


    for(int x=0;x<Pred_Mass->GetNbinsY()+1;x++){
       for(unsigned int pe=0;pe<100;pe++){
          Pred_Prof_Mass    ->SetBinContent(x, pe, 0);
          Pred_Prof_MassTOF ->SetBinContent(x, pe, 0);
          Pred_Prof_MassComb->SetBinContent(x, pe, 0);
       }
    }



      TRandom3* RNG = new TRandom3();
      printf("Predicting (%4i / %4i)     :",CutIndex+1,(int)HCuts_Pt->GetNbinsX());
      int TreeStep = 100/50;if(TreeStep==0)TreeStep=1;
      for(unsigned int pe=0;pe<100;pe++){    
      if(pe%TreeStep==0){printf(".");fflush(stdout);}

      TH1D* tmpH_Mass     =  new TH1D("tmpH_Mass"    ,"tmpH_Mass"    ,MassNBins,0,MassHistoUpperBound);
      TH1D* tmpH_MassTOF  =  new TH1D("tmpH_MassTOF" ,"tmpH_MassTOF" ,MassNBins,0,MassHistoUpperBound);
      TH1D* tmpH_MassComb =  new TH1D("tmpH_MassComb","tmpH_MassComb",MassNBins,0,MassHistoUpperBound);


      double PE_A=RNG->Poisson(A);
      double PE_B=RNG->Poisson(B);
      double PE_C=RNG->Poisson(C);
      //double PE_D=RNG->Poisson(D);
      double PE_E=RNG->Poisson(E);
      double PE_F=RNG->Poisson(F);
      double PE_G=RNG->Poisson(G);
      //double PE_H=RNG->Poisson(H);
      double PE_P = 0;

      if(E>0){
         PE_P    = (PE_E>0 ? (PE_A*PE_F*PE_G)/(PE_E*PE_E) : 0);
      }else if(A>0){
         PE_P    = (PE_A>0 ? ((PE_C*PE_B)/PE_A) : 0);
      }

      for(int i=0;i<Pred_EtaB_Proj_PE->GetNbinsX()+1;i++){Pred_EtaB_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaB_Proj->GetBinContent(i)) );}    Pred_EtaB_Proj_PE->Scale(1.0/Pred_EtaB_Proj_PE->Integral());
      for(int i=0;i<Pred_EtaS_Proj_PE->GetNbinsX()+1;i++){Pred_EtaS_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaS_Proj->GetBinContent(i)) );}    Pred_EtaS_Proj_PE->Scale(1.0/Pred_EtaS_Proj_PE->Integral());
      for(int i=0;i<Pred_EtaS2_Proj_PE->GetNbinsX()+1;i++){Pred_EtaS2_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaS2_Proj->GetBinContent(i)) );} Pred_EtaS2_Proj_PE->Scale(1.0/Pred_EtaS2_Proj_PE->Integral());


      for(int i=0;i<Pred_EtaPWeighted_PE->GetNbinsX()+1;i++){
      for(int j=0;j<Pred_EtaPWeighted_PE->GetNbinsY()+1;j++){
         Pred_EtaPWeighted_PE->SetBinContent(i,j,RNG->Poisson(Pred_EtaPWeighted->GetBinContent(i,j)));
      }}

      double WeightP = 0.0;
      for(int x=0;x<=Pred_EtaPWeighted_PE->GetXaxis()->GetNbins();x++){
         WeightP = 0.0;
         if(Pred_EtaB_Proj_PE->GetBinContent(x)>0){
                           WeightP = Pred_EtaS_Proj_PE ->GetBinContent(x)/Pred_EtaB_Proj_PE->GetBinContent(x);
            if(TypeMode==2)WeightP*= Pred_EtaS2_Proj_PE->GetBinContent(x)/Pred_EtaB_Proj_PE->GetBinContent(x);
         }

         for(int y=0;y<=Pred_EtaPWeighted_PE->GetYaxis()->GetNbins();y++){
            Pred_EtaPWeighted_PE->SetBinContent(x,y,Pred_EtaPWeighted_PE->GetBinContent(x,y)*WeightP);
         }
      }
      TH1D* Pred_P_ProjPE = Pred_EtaPWeighted_PE->ProjectionY("Pred_P_ProjPE");                                                        Pred_P_ProjPE->Scale(1.0/Pred_P_ProjPE->Integral());
      for(int i=0;i<Pred_I_ProjPE->GetNbinsX()+1;i++){Pred_I_ProjPE->SetBinContent(i,RNG->Poisson(Pred_I_Proj->GetBinContent(i)) );}   Pred_I_ProjPE->Scale(1.0/Pred_I_ProjPE->Integral());
      for(int i=0;i<Pred_T_ProjPE->GetNbinsX()+1;i++){Pred_T_ProjPE->SetBinContent(i,RNG->Poisson(Pred_T_Proj->GetBinContent(i)) );}   Pred_T_ProjPE->Scale(1.0/Pred_T_ProjPE->Integral());

      double Proba, MI, MComb;//, MT=0, ProbaT=0;
      for(int x=0;x<Pred_P_ProjPE->GetNbinsX()+1;x++){    if(Pred_P_ProjPE->GetBinContent(x)<=0.0){continue;}  const double& p = Pred_P_ProjPE->GetBinCenter(x);
      for(int y=0;y<Pred_I_ProjPE->GetNbinsX()+1;y++){    if(Pred_I_ProjPE->GetBinContent(y)<=0.0){continue;}  const double& i = Pred_I_ProjPE->GetBinCenter(y);
         Proba = Pred_P_ProjPE->GetBinContent(x) * Pred_I_ProjPE->GetBinContent(y);  if(Proba<=0 || isnan((float)Proba))continue;
         MI = GetMass(p,i);
         MComb = MI;
         tmpH_Mass->Fill(MI,Proba);

//         if(TypeMode==2){
//         for(int z=0;z<Pred_T_ProjPE->GetNbinsX()+1;z++){   if(Pred_T_ProjPE->GetBinContent(z)<=0.0){continue;}   const double& t = Pred_T_ProjPE->GetBinCenter(z);
//            ProbaT = Proba * Pred_T_ProjPE->GetBinContent(z);  if(ProbaT<=0 || isnan(ProbaT))continue;
//            MT = GetTOFMass(p,t);
//            tmpH_MassTOF->Fill(MT,ProbaT);
//            MComb = GetMassFromBeta(p, (GetIBeta(i) + (1/t))*0.5 );        
//            tmpH_MassComb->Fill(MComb,ProbaT);
//         }}else{
            tmpH_MassComb->Fill(MComb,Proba);
//         }
      }}

//      printf("PE_P = %f\n",PE_P);

      for(int x=0;x<tmpH_Mass->GetNbinsX()+1;x++){
         //const double& M = tmpH_Mass->GetXaxis()->GetBinCenter(x);
         Pred_Prof_Mass    ->SetBinContent(x, pe, tmpH_Mass    ->GetBinContent(x) * PE_P);
         Pred_Prof_MassTOF ->SetBinContent(x, pe, tmpH_MassTOF ->GetBinContent(x) * PE_P);
         Pred_Prof_MassComb->SetBinContent(x, pe, tmpH_MassComb->GetBinContent(x) * PE_P);
         if(isnan((float)(tmpH_Mass    ->GetBinContent(x) * PE_P))){printf("%f x %f\n",tmpH_Mass    ->GetBinContent(x),PE_P); fflush(stdout);exit(0);}
      }
     
      delete Pred_P_ProjPE;
      delete tmpH_Mass;
      delete tmpH_MassTOF;
      delete tmpH_MassComb;
     }printf("\n");

    for(int x=0;x<Pred_Mass->GetNbinsY()+1;x++){
//       Pred_Mass    ->SetBinContent(CutIndex+1,x,Pred_Prof_Mass    ->GetBinContent(x)); Pred_Mass      ->SetBinError(CutIndex+1,x,sqrt(pow(Pred_Prof_Mass    ->GetBinError(x),2) + Pred_Prof_Mass    ->GetBinContent(x)*(Perr/P)));
//       Pred_MassTOF ->SetBinContent(CutIndex+1,x,Pred_Prof_MassTOF ->GetBinContent(x)); Pred_MassTOF   ->SetBinError(CutIndex+1,x,sqrt(pow(Pred_Prof_MassTOF ->GetBinError(x),2) + Pred_Prof_MassTOF ->GetBinContent(x)*(Perr/P)));
//       Pred_MassComb->SetBinContent(CutIndex+1,x,Pred_Prof_MassComb->GetBinContent(x)); Pred_MassComb  ->SetBinError(CutIndex+1,x,sqrt(pow(Pred_Prof_MassComb->GetBinError(x),2) + Pred_Prof_MassComb->GetBinContent(x)*(Perr/P)));

       double Mean=0, MeanTOF=0, MeanComb=0;
       for(unsigned int pe=0;pe<100;pe++){
	 //if(CutIndex==4){printf("Bin=%4i pe=%3i --> BinCOntent=%f\n",x,pe,Pred_Prof_Mass    ->GetBinContent(x, pe));}
          Mean     += Pred_Prof_Mass    ->GetBinContent(x, pe);
          MeanTOF  += Pred_Prof_MassTOF ->GetBinContent(x, pe);
          MeanComb += Pred_Prof_MassComb->GetBinContent(x, pe);
       }Mean/=100.0; MeanTOF/=100.0;  MeanComb/=100.0;

       //if(CutIndex==4){printf("MEAN = %f\n",Mean);}


       double Err=0, ErrTOF=0, ErrComb=0;
       for(unsigned int pe=0;pe<100;pe++){
	  //if(CutIndex==4){printf("Bin=%4i pe=%3i --> DeltaM=%f\n",x,pe,sqrt(pow(Mean     - Pred_Prof_Mass    ->GetBinContent(x, pe),2)));}
          Err     += pow(Mean     - Pred_Prof_Mass    ->GetBinContent(x, pe),2);
          ErrTOF  += pow(MeanTOF  - Pred_Prof_MassTOF ->GetBinContent(x, pe),2);
          ErrComb += pow(MeanComb - Pred_Prof_MassComb->GetBinContent(x, pe),2);
       }Err=sqrt(Err/99.0); ErrTOF=sqrt(ErrTOF/99.0);  ErrComb=sqrt(ErrComb/99.0);
       //if(CutIndex==4){printf("ERROR = %f\n",Err);}


       Pred_Mass    ->SetBinContent(CutIndex+1,x,Mean    ); Pred_Mass      ->SetBinError(CutIndex+1,x,Err    );
       Pred_MassTOF ->SetBinContent(CutIndex+1,x,MeanTOF ); Pred_MassTOF   ->SetBinError(CutIndex+1,x,ErrTOF );
       Pred_MassComb->SetBinContent(CutIndex+1,x,MeanComb); Pred_MassComb  ->SetBinError(CutIndex+1,x,ErrComb);
    }
//    printf("MassInt %f\n",Pred_Prof_Mass->Integral());


    delete Pred_EtaB_Proj_PE;
    delete Pred_EtaS_Proj_PE;
    delete Pred_EtaS2_Proj_PE;

    delete Pred_Prof_Mass;
    delete Pred_Prof_MassTOF;
    delete Pred_Prof_MassComb;
    delete Pred_EtaPWeighted_PE;
    delete Pred_I_ProjPE;
    delete Pred_T_ProjPE;

//    delete Pred_P_PDF;
//    delete Pred_I_PDF;
//    delete Pred_T_PDF;
//    delete Pred_P_Proj;
    delete Pred_I_Proj;
    delete Pred_T_Proj;
    delete Pred_EtaB_Proj;
    delete Pred_EtaS_Proj;
    delete Pred_EtaS2_Proj;
    delete Pred_EtaPWeighted;

   }


   //////////////////////////////////////////////////     DUMP USEFUL INFORMATION
   char Buffer[2048];
   if(MODE=="ANALYSE_DATA"){sprintf(Buffer,"%s/Info.txt",SavePath.c_str());
   }else{                    sprintf(Buffer,"%s/Info_MC.txt",SavePath.c_str());}
   FILE* pFile = fopen(Buffer,"w");
   fprintf(pFile,"Selection      = %s\n",dEdxS_Label.c_str());
   fprintf(pFile,"Mass           = %s\n",dEdxM_Label.c_str());
   fprintf(pFile,"TOF            = %s\n",TOF_Label.c_str());
   fprintf(pFile,"|eta|          < %f\n",GlobalMaxEta);
   fprintf(pFile,"#Hit           > %02i\n",GlobalMinNOH);
   fprintf(pFile,"#dEdx Hit      > %02i\n",GlobalMinNOM);
   fprintf(pFile,"nDoF           > %02i\n",GlobalMinNOH);
   fprintf(pFile,"Chi2/ndf       < %6.2f\n",GlobalMaxChi2);
   fprintf(pFile,"SumPt          < %6.2f\n",GlobalMaxTIsol);
   fprintf(pFile,"E/p            < %6.2f\n",GlobalMaxEIsol);

   for(unsigned int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){
      const double& A=H_A->GetBinContent(CutIndex+1);
      const double& B=H_B->GetBinContent(CutIndex+1);
      const double& C=H_C->GetBinContent(CutIndex+1);
      const double& D=H_D->GetBinContent(CutIndex+1);
      const double& E=H_E->GetBinContent(CutIndex+1);
      const double& F=H_F->GetBinContent(CutIndex+1);
      const double& G=H_G->GetBinContent(CutIndex+1);
      const double& H=H_H->GetBinContent(CutIndex+1);

      //fprintf(pFile  ,"CutIndex=%4i --> (Pt>%6.2f I>%6.3f TOF>%6.3f) Ndata=%+6.2E  NPred=%6.3E+-%6.3E <--> A=%6.2E B=%6.E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E\n",CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1), HCuts_I  ->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1), D,H_P->GetBinContent(CutIndex+1),H_P->GetBinError(CutIndex+1) ,A, B, C, D, E, F, G, H);
   }
   fprintf(pFile,"--------------------\n");
   fclose(pFile);
   //////////////////////////////////////////////////     CREATE EFFICIENCY FILE

   fflush(stdout);
   InputFile->cd("Data");
   H_P->Write();
   H_P_Cen->Write();
   H_P_For->Write();
   Pred_Mass->Write();
   Pred_MassTOF->Write();
   Pred_MassComb->Write();
   //InputFile->Write();
   InputFile->Close();
}
