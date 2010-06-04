#void Analysis_Step2345(string MODE="COMPILE", double WP_Pt=-1.0, double WP_I=-1, int SplitMode_=2, int dEdxSel_=0, int dEdxMass_=0, int TypeMode_=0)
root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy = makeshared.ReplaceAll("-W ", "");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gSystem->Load("libDataFormatsVertexReco.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  .x Analysis_Step2345.C++("ANALYSE",-2.0,-2.5,2,5,0,0)
EOF

