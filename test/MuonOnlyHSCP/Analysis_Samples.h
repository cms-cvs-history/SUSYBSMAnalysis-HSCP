#ifndef HSCP_ANALYSIS_SAMPLE
#define HSCP_ANALYSIS_SAMPLE

#define SID_GL300     0
#define SID_GL400     1
#define SID_GL500     2
#define SID_GL600     3
#define SID_GL700     4
#define SID_GL800     5
#define SID_GL900     6
#define SID_GL1000    7
#define SID_GL1100    8
#define SID_GL300N    9
#define SID_GL400N    10
#define SID_GL500N    11
#define SID_GL600N    12
#define SID_GL700N    13
#define SID_GL800N    14
#define SID_GL900N    15
#define SID_GL1000N   16
#define SID_GL1100N   17
#define SID_ST130     18
#define SID_ST200     19
#define SID_ST300     20
#define SID_ST400     21
#define SID_ST500     22
#define SID_ST600     23
#define SID_ST700     24
#define SID_ST800     25
#define SID_ST130N    26
#define SID_ST200N    27
#define SID_ST300N    28
#define SID_ST400N    29
#define SID_ST500N    30
#define SID_ST600N    31
#define SID_ST700N    32
#define SID_ST800N    33
#define SID_GS100     34
#define SID_GS126     35
#define SID_GS156     36
#define SID_GS200     37
#define SID_GS247     38
#define SID_GS308     39
#define SID_GS370     40
#define SID_GS432     41
#define SID_GS494     42
#define SID_PS100     43
#define SID_PS126     44
#define SID_PS156     45
#define SID_PS200     46
#define SID_PS247     47
#define SID_PS308     48
#define SID_D08K100   49
#define SID_D08K121   50
#define SID_D08K182   51
#define SID_D08K242   52
#define SID_D08K302   53
#define SID_D08K350   54
#define SID_D08K370   55
#define SID_D08K390   56
#define SID_D08K395   57
#define SID_D08K400   58
#define SID_D08K410   59
#define SID_D08K420   60
#define SID_D08K500   61
#define SID_D12K100   62
#define SID_D12K182   63
#define SID_D12K302   64
#define SID_D12K500   65
#define SID_D12K530   66
#define SID_D12K570   67
#define SID_D12K590   68
#define SID_D12K595   69
#define SID_D12K600   70
#define SID_D12K610   71
#define SID_D12K620   72
#define SID_D12K700   73
#define SID_D16K100   74
#define SID_D16K182   75
#define SID_D16K302   76
#define SID_D16K500   77
#define SID_D16K700   78
#define SID_D16K730   79
#define SID_D16K770   80
#define SID_D16K790   81
#define SID_D16K795   82
#define SID_D16K800   83
#define SID_D16K820   84
#define SID_D16K900   85


int                  RunningPeriods = 2;
double               IntegratedLuminosity = 3608; //2410;//2125; //2080; //1912; //1947; //1631; //976.204518023; //705.273820; //342.603275; //204.160928; //191.04;
float                Event_Weight = 1;
int                  MaxEntry = -1;


class stSignal{
   public:
   std::string Type;
   std::string Name;
   std::string FileName;
   std::string Legend;
   double Mass;
   double XSec;
   bool   MakePlot;
   bool   IsS4PileUp;

   stSignal(); 
      stSignal(std::string Type_, std::string Name_, std::string FileName_, std::string Legend_, double Mass_, bool MakePlot_, bool IsS4PileUp_, double XSec_){Type=Type_; Name=Name_; FileName=FileName_; Legend=Legend_; Mass=Mass_; MakePlot=MakePlot_; IsS4PileUp=IsS4PileUp_;XSec=XSec_;}
};


void GetSignalDefinition(std::vector<stSignal>& signals, std::string File=""){
  if(File=="" || File=="Gluino300") signals.push_back(stSignal("Gluino", "Gluino300", "Gluino300"    , "#tilde{g} 300"                 , 300,  1, 1,  65.800000) ); //NLO
  if(File=="" || File=="Gluino400") signals.push_back(stSignal("Gluino", "Gluino400", "Gluino400"    , "#tilde{g} 400"                 , 400,  0, 1,   11.20000) ); //NLO
  if(File=="" || File=="Gluino500") signals.push_back(stSignal("Gluino", "Gluino500", "Gluino500"    , "#tilde{g} 500"                 , 500,  1, 1,   2.540000) ); //NLO
  if(File=="" || File=="Gluino600") signals.push_back(stSignal("Gluino", "Gluino600", "Gluino600"    , "#tilde{g} 600"                 , 600,  0, 1,   0.693000) ); //NLO
  if(File=="" || File=="Gluino700") signals.push_back(stSignal("Gluino", "Gluino700", "Gluino700"    , "#tilde{g} 700"                 , 700,  0, 1,   0.214000) ); //NLO
  if(File=="" || File=="Gluino800") signals.push_back(stSignal("Gluino", "Gluino800", "Gluino800"    , "#tilde{g} 800"                 , 800,  1, 1,   0.072500) ); //NLO
  if(File=="" || File=="Gluino900") signals.push_back(stSignal("Gluino", "Gluino900", "Gluino900"    , "#tilde{g} 900"                 , 900,  0, 1,   0.026200) ); //NLO
  if(File=="" || File=="Gluino1000") signals.push_back(stSignal("Gluino", "Gluino1000", "Gluino1000"  , "#tilde{g} 1000"                ,1000,  1, 1,   0.0098700) ); //NLO
  if(File=="" || File=="Gluino1100") signals.push_back(stSignal("Gluino", "Gluino1100", "Gluino1100"  , "#tilde{g} 1100"                ,1100,  1, 1,   0.0038600) ); //NLO 
  if(File=="" || File=="Gluino1200") signals.push_back(stSignal("Gluino", "Gluino1200", "Gluino1200"   , "#tilde{g} 1200"                ,1200,  1, 1,   0.0015400) ); //NLO

  if(File=="" || File=="GMStau100") signals.push_back(stSignal("Stau"  , "GMStau100", "stau_M-100"    , "GMSB #tilde{#tau}_{1} 100"     , 100,  0, 1,   1.3398) );
  if(File=="" || File=="GMStau126") signals.push_back(stSignal("Stau"  , "GMStau126", "stau_M-126"    , "GMSB #tilde{#tau}_{1} 126"     , 126,  0, 1,   0.274591) );
  if(File=="" || File=="GMStau156") signals.push_back(stSignal("Stau"  , "GMStau156", "stau_M-156"    , "GMSB #tilde{#tau}_{1} 156"     , 156,  0, 1,  0.0645953) );
  if(File=="" || File=="GMStau200") signals.push_back(stSignal("Stau"  , "GMStau200", "stau_M-200"    , "GMSB #tilde{#tau}_{1} 200"     , 200,  1, 1,   0.0118093) );
  if(File=="" || File=="GMStau247") signals.push_back(stSignal("Stau"  , "GMStau247", "stau_M-247"    , "GMSB #tilde{#tau}_{1} 247"     , 247,  0, 1,  0.00342512) );
  if(File=="" || File=="GMStau308") signals.push_back(stSignal("Stau"  , "GMStau308", "stau_M-308"    , "GMSB #tilde{#tau}_{1} 308"     , 308,  0, 1,  0.00098447 ) );
  if(File=="" || File=="GMStau370") signals.push_back(stSignal("Stau"  , "GMStau370", "stau_M-370"    , "GMSB #tilde{#tau}_{1} 370"     , 370,  0, 1,   0.000353388) );
  if(File=="" || File=="GMStau432") signals.push_back(stSignal("Stau"  , "GMStau432", "stau_M-432"    , "GMSB #tilde{#tau}_{1} 432"     , 432,  0, 1,   0.000141817) );
  if(File=="" || File=="GMStau494") signals.push_back(stSignal("Stau"  , "GMStau494", "stau_M-494"    , "GMSB #tilde{#tau}_{1} 494"     , 494,  1, 1,   0.00006177) );

  if(File=="" || File=="Stop130") signals.push_back(stSignal("Stop"  , "Stop130", "stop_M-130"      , "#tilde{t}_{1} 130"             , 130,  1, 1, 120.000000) ); //NLO
  if(File=="" || File=="Stop200") signals.push_back(stSignal("Stop"  , "Stop200", "stop_M-200"      , "#tilde{t}_{1} 200"             , 200,  0, 1,  13.000000) ); //NLO
  if(File=="" || File=="Stop300") signals.push_back(stSignal("Stop"  , "Stop300", "stop_M-300"      , "#tilde{t}_{1} 300"             , 300,  0, 1,   1.310000) ); //NLO
  if(File=="" || File=="Stop400") signals.push_back(stSignal("Stop"  , "Stop400", "stop_M-400"      , "#tilde{t}_{1} 400"             , 400,  0, 1,   0.218000) ); //NLO
  if(File=="" || File=="Stop500") signals.push_back(stSignal("Stop"  , "Stop500", "stop_M-500"      , "#tilde{t}_{1} 500"             , 500,  1, 1,  0.047800) ); //NLO
  if(File=="" || File=="Stop600") signals.push_back(stSignal("Stop"  , "Stop600", "stop_M-600"      , "#tilde{t}_{1} 600"             , 600,  0, 1,   0.012500) ); //NLO
  if(File=="" || File=="Stop700") signals.push_back(stSignal("Stop"  , "Stop700", "stop_M-700"      , "#tilde{t}_{1} 700"             , 700,  0, 1,   0.003560) ); //NLO
  if(File=="" || File=="Stop800") signals.push_back(stSignal("Stop"  , "Stop800", "stop_M-800"      , "#tilde{t}_{1} 800"             , 800,  1, 1,   0.001140) ); //NLO
}

struct stMC{
   std::string Name;
   double XSection;
   double MaxPtHat;
   double MaxEvent;
   bool   IsS4PileUp;

   stMC();
      stMC(std::string Name_, double XSection_, double MaxPtHat_, int MaxEvent_, bool IsS4PileUp_){Name = Name_; XSection = XSection_; MaxPtHat = MaxPtHat_; MaxEvent = MaxEvent_;IsS4PileUp = IsS4PileUp_;}
};

void GetMCDefinition(std::vector<stMC>& MC){

   MC.push_back(stMC("MC_DYToTauTau"            ,     1.300E3  , -1, -1, 0));
   MC.push_back(stMC("MC_DYToMuMu"              ,     1.300E3  , -1, -1, 0));
   MC.push_back(stMC("MC_WJetsToLNu"            ,     2.777E4  , -1, -1, 1));
   MC.push_back(stMC("MC_TTJets"                ,     9.400E1  , -1, -1, 1));
   MC.push_back(stMC("MC_QCD_Pt-15to30"         ,     8.16E8  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-30to50"         ,     5.310E7  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-50to80"         ,     6.360E6  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-80to120"        ,     7.840E5  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-120to170"       ,     1.150E5  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-170to300"       ,     2.430E4  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-300to470"       ,     1.170E3  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-470to600"       ,     7.020E1  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-600to800"       ,     1.560E1  , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-800to1000"      ,     1.84     , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-1000to1400"     ,     3.320E-1 , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-1400to1800"     ,     1.090E-2 , -1, -1, 0));
   MC.push_back(stMC("MC_QCD_Pt-1800"           ,     3.580E-4 , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-0to15"   ,     4.280E3  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-15to20"  ,     1.450E2  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-20to30"  ,     1.310E2  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-30to50"  ,     8.400E1  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-50to80"  ,     3.220E1  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-80to120" ,     9.98     , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-120to170",     2.73     , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-170to230",     7.21E-1  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-230to300",     1.94E-1  , -1, -1, 0));
   MC.push_back(stMC("MC_ZJetToMuMu_Pt-300"     ,     7.59E-2  , -1, -1, 0));
   MC.push_back(stMC("MC_ZZ"                    ,     4.287    , -1, -1, 1));
   MC.push_back(stMC("MC_WW"                    ,     2.783E1  , -1, -1, 1));
   MC.push_back(stMC("MC_WZ"                    ,     1.47E1   , -1, -1, 1));
}

void GetInputFiles(std::vector<std::string>& inputFiles, std::string SampleName, std::string Runs=""){
  std::string BaseDirectory = "dcache:/pnfs/cms/WAX/11/store/user/farrell3/30May2012HSCPEDMFiles/";
  //std::string BaseDirectory = "dcache:/pnfs/cms/WAX/11/store/user/venkat12/2012Data/";
  //std::string BaseDirectory = "/uscmst1b_scratch/lpc1/3DayLifetime/farrell/SAMuonOnlyHSCP/";
   if(SampleName=="Data"){
     inputFiles.push_back(BaseDirectory + "Data_" + Runs + ".root");
     //inputFiles.push_back("../BuildHSCParticles/Data/Merge.root");
   }else if(SampleName.find("Cosmic",0)<std::string::npos){
     //inputFiles.push_back("../BuildHSCParticles/Data/Merge.root");
     inputFiles.push_back(BaseDirectory + "Data_" + Runs + ".root");
   }else{
     inputFiles.push_back(BaseDirectory + SampleName + "BX1.root");
     //inputFiles.push_back(BaseDirectory + "Merge.root");
   }
}

#endif

