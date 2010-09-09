#void Analysis_Step234(string MODE="COMPILE", double WP_Pt=-1.0, double WP_I=-1, int SplitMode_=2, int dEdxSel_=0, int dEdxMass_=0, int TypeMode_=0)
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
  .x Analysis_Step234.C++("COMPILE",-1.5,-1.5,2,11,3,0)
EOF

