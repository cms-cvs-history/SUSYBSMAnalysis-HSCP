import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
  skipEvents = cms.untracked.uint32(0) ,
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
"file:/data1/arizzi/HSCP/CMSSW_2_2_3/src/cruzet-hscp.root"
    )
)

process.tofAnalysis  = cms.EDAnalyzer('CosmicTOFAnalyzer'
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)


process.p = cms.Path(process.tofAnalysis)
