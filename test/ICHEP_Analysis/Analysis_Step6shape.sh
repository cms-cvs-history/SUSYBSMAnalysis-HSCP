#! /bin/sh
####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

export SCRAM_ARCH=slc5_amd64_gcc462
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/nfs/soft/cms
cd /afs/cern.ch/work/q/querten/public/12_04_16_HSCP_2012/CMSSW_5_2_4/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis
eval `scramv1 runtime -sh`
root -l -b << EOF
   TString makeshared(gSystem->GetMakeSharedLib());
   TString dummy = makeshared.ReplaceAll("-W ", "");
   TString dummy = makeshared.ReplaceAll("-Wshadow ", " -std=c++0x ");
   gSystem->SetMakeSharedLib(makeshared);
   gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
//   .x /afs/cern.ch/work/q/querten/public/12_04_16_HSCP_2012/CMSSW_5_2_4/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/Analysis_Step6shape.C+("ANALYSE", "Results/dedxASmi/combined/Eta15/PtMin45/Type0/", "GMStau247", "GMStau247", 1.0, 1.0, 1.0, "")
   .x /afs/cern.ch/work/q/querten/public/12_04_16_HSCP_2012/CMSSW_5_2_4/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/Analysis_Step6shape.C+("Final", "Results/dedxASmi/combined/Eta15/PtMin45/Type0/", "GMStau247", "GMStau247", 1.0, 1.0, 1.0, "")
   .q
EOF

mv HscpLimits* /afs/cern.ch/work/q/querten/public/12_04_16_HSCP_2012/CMSSW_5_2_4/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/FARM/outputs/
