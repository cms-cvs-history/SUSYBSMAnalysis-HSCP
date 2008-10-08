import FWCore.ParameterSet.Config as cms

process = cms.Process("HSCP")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring (
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/DEB4851F-3E73-DD11-8375-003048769E6D.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/16272167-2473-DD11-A8B3-001A92810AD2.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/3E5B0667-3873-DD11-A4CE-001731AF66F1.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/9E723A39-D472-DD11-B2CA-001731AF685B.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/88F3F8BD-3973-DD11-9496-001731AF6943.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/C09A6175-1773-DD11-9554-0018F3D09600.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/9618753D-C572-DD11-A065-001731AF67EF.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/4C411E3D-C072-DD11-A3BE-001A92971B7E.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/9A9F8A46-F872-DD11-83F6-001731AF687F.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/328F340C-C072-DD11-B2FC-0018F3D09616.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/E6DA51C2-2373-DD11-B98F-001731AF66C2.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/2EC4B168-D372-DD11-8D19-001A92971B36.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/944F6914-E372-DD11-8704-0017313F01E4.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/F6FA3EB3-3D73-DD11-8F6E-00304876A0DB.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/90576D0C-C072-DD11-B6EF-0018F3D0962E.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0001/EAAC08DB-1E73-DD11-8C29-0018F3D096EE.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/C62D04E6-3D73-DD11-99BC-003048767ED1.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/20B5CB54-D672-DD11-868D-001A92811748.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/B090F928-2273-DD11-B5A5-001A92810AA4.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/88F55C62-2573-DD11-B955-00304875A9ED.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/326E9CD9-BF72-DD11-A073-001731AF66B7.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/EE2E3295-D372-DD11-BE02-001731AF677F.root',
'/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/3AF07B6C-2273-DD11-BB92-001A92971B9A.root'
)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('cruzet-hscp2.root'),
    outputCommands = cms.untracked.vstring( 'drop *', 'keep *_*_*_HSCP', 'keep recoTracks_*_*_*', 'keep recoMuons_*_*_*', 'keep recoTrackExtras_*_*_*')

)


#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('cruzet-hscp.root')
#)

process.MessageLogger = cms.Service("MessageLogger")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)

# Magnetic fiuld: force mag field to be 0.0 tesla
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")

# patch needed for CRUZET2 but not used in CRUZET3 cfg
process.SteppingHelixPropagatorAny.useInTeslaFromMagField = True
process.SteppingHelixPropagatorAlong.useInTeslaFromMagField = True
process.SteppingHelixPropagatorOpposite.useInTeslaFromMagField = True
process.SteppingHelixPropagatorAny.SetVBFPointer = True
process.SteppingHelixPropagatorAlong.SetVBFPointer = True
process.SteppingHelixPropagatorOpposite.SetVBFPointer = True
process.VolumeBasedMagneticFieldESProducer.label = 'VolumeBasedMagneticField'




process.prefer("GlobalTag")

from SUSYBSMAnalysis.HSCP.MuonSegmentMatcher_cff import *

process.betaFromTOF = cms.EDFilter("BetaFromTOF",
    MuonSegmentMatcher,
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny', 
            'PropagatorWithMaterial', 
            'PropagatorWithMaterialOpposite'),
        RPCLayers = cms.bool(True)
    ),
    DTsegments = cms.untracked.InputTag("dt4DSegments"),
    PruneCut = cms.double(0.1),
    HitsMin = cms.int32(3),
    debug = cms.bool(False),
    Muons = cms.untracked.InputTag("STAMuonsBarrelOnly")
)

process.load("SUSYBSMAnalysis.HSCP.ecalCosmicTrackTimingProducer_cfi")


process.p = cms.Path(process.betaFromTOF*process.ecalCosmicTrackTimingProducer)
process.e = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'CRUZET4_V2P::All'

