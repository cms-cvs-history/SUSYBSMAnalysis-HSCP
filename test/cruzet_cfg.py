import FWCore.ParameterSet.Config as cms

process = cms.Process("HSCP")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:BE59E2FA-E072-DD11-8DC9-00304875A9E5.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('cruzet-hscp.root')
)

process.MessageLogger = cms.Service("MessageLogger")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)

process.prefer("GlobalTag")

process.load("SUSYBSMAnalysis.HSCP.betaFromTOF_cfi")

process.betaFromTOF.Muons = "STAMuonsBarrelOnly"

process.p = cms.Path(process.betaFromTOF)
process.e = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'CRUZET4_V2P::All'

