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
  //.x Analysis_Step234.C++("CUTFINDER",-0.5,-0.5,-0.5 ,1,"dedxASmi","dedxCNPHarm2",1);
  .x Analysis_Step234.C++("PLOT"     ,-0.6,-0.6,-0.6 ,0,"dedxASmi","dedxCNPHarm2",2);
EOF

