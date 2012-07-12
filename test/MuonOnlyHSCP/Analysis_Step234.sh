root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy = makeshared.ReplaceAll("-W ", "");
  TString dummy = makeshared.ReplaceAll("-Wshadow ", " -std=c++0x ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gSystem->Load("libDataFormatsVertexReco.so");
  gSystem->Load("libDataFormatsTrackReco.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  gSystem->Load("libPhysicsToolsUtilities.so");
  gSystem->Load("libCondFormatsJetMETObjects.so");
  gSystem->Load("libJetMETCorrectionsObjects.so");
  gSystem->Load("libGeometryDTGeometry.so");
  .x Analysis_Step234.C++("ANALYSE_DATA", "195390_195749", 70.0, 2.4);
  //.x Analysis_Step234.C++("ANALYSE_SIGNAL", "Gluino1000", 70.0, 2.4);
  //.x Analysis_Step234.C++("ANALYSE_COSMIC", "192001_193751", 70.0, 2.4);
EOF

