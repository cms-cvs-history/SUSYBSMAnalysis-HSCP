import FWCore.ParameterSet.Config as cms

####################################################################################
#   BEAMSPOT + TRAJECTORY BUILDERS
####################################################################################

from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
from RecoTracker.TrackProducer.TrackRefitters_cff import *


####################################################################################
#   DEDX ESTIMATORS 
####################################################################################

dedxHarm2 = cms.EDProducer("DeDxEstimatorProducer",
    tracks                     = cms.InputTag("TrackRefitter"),
    trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),

    estimator      = cms.string('generic'),
    exponent       = cms.double(-2.0),

    UseStrip       = cms.bool(True),
    UsePixel       = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*250),
    MeVperADCPixel = cms.double(3.61e-06),

    MisCalib_Mean      = cms.untracked.double(1.0),
    MisCalib_Sigma     = cms.untracked.double(0.00),

    UseCalibration  = cms.bool(False),
    calibrationPath = cms.string(""),
    ShapeTest       = cms.bool(False),
)

dedxTru40 = cms.EDProducer("DeDxEstimatorProducer",
    tracks                     = cms.InputTag("TrackRefitter"),
    trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),

    estimator      = cms.string('truncated'),
    fraction       = cms.double(0.4),

    UseStrip       = cms.bool(True),
    UsePixel       = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*250),
    MeVperADCPixel = cms.double(3.61e-06),

    MisCalib_Mean      = cms.untracked.double(1.0),
    MisCalib_Sigma     = cms.untracked.double(0.00),

    UseCalibration  = cms.bool(False),
    calibrationPath = cms.string(""),
    ShapeTest       = cms.bool(False),
)


dedxMed = cms.EDProducer("DeDxEstimatorProducer",
    tracks                     = cms.InputTag("TrackRefitter"),
    trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),

    estimator      = cms.string('median'),

    UseStrip       = cms.bool(True),
    UsePixel       = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*250),
    MeVperADCPixel = cms.double(3.61e-06),

    MisCalib_Mean      = cms.untracked.double(1.0),
    MisCalib_Sigma     = cms.untracked.double(0.00),

    UseCalibration  = cms.bool(False),
    calibrationPath = cms.string(""),
    ShapeTest       = cms.bool(False),
)

dedxNPHarm2                  = dedxHarm2.clone()
dedxNPHarm2.UsePixel         = cms.bool(False)

dedxNPTru40                  = dedxTru40.clone()
dedxNPTru40.UsePixel         = cms.bool(False)

dedxNPMed                    = dedxMed.clone()
dedxNPMed.UsePixel           = cms.bool(False)

dedxCHarm2                   = dedxHarm2.clone()
dedxCHarm2.UseCalibration    = cms.bool(True)
dedxCHarm2.calibrationPath   = cms.string("file:Gains.root")

dedxCTru40                   = dedxTru40.clone()
dedxCTru40.UseCalibration    = cms.bool(True)
dedxCTru40.calibrationPath   = cms.string("file:Gains.root")

dedxCMed                     = dedxMed.clone()
dedxCMed.UseCalibration      = cms.bool(True)
dedxCMed.calibrationPath     = cms.string("file:Gains.root")

dedxCNPHarm2                 = dedxNPHarm2.clone()
dedxCNPHarm2.UseCalibration  = cms.bool(True)
dedxCNPHarm2.calibrationPath = cms.string("file:Gains.root")

dedxCNPTru40                 = dedxNPTru40.clone()
dedxCNPTru40.UseCalibration  = cms.bool(True)
dedxCNPTru40.calibrationPath = cms.string("file:Gains.root")

dedxCNPMed                   = dedxNPMed.clone()
dedxCNPMed.UseCalibration    = cms.bool(True)
dedxCNPMed.calibrationPath   = cms.string("file:Gains.root")

dedxSTCNPHarm2                = dedxCNPHarm2.clone()
dedxSTCNPHarm2.ShapeTest      = cms.bool(True)

dedxSTCNPTru40                = dedxCNPTru40.clone()
dedxSTCNPTru40.ShapeTest      = cms.bool(True)

dedxSTCNPMed                  = dedxCNPMed.clone()
dedxSTCNPMed.ShapeTest        = cms.bool(True)

####################################################################################
#   DEDX DISCRIMINATORS 
####################################################################################

dedxProd               = cms.EDProducer("DeDxDiscriminatorProducer",
    tracks                     = cms.InputTag("TrackRefitter"),
    trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),

    Reccord            = cms.untracked.string("SiStripDeDxMip_3D_Rcd"),
    Formula            = cms.untracked.uint32(0),
#    ProbabilityMode    = cms.untracked.string("Integral"),
    ProbabilityMode    = cms.untracked.string("Accumulation"),


    UseStrip           = cms.bool(True),
    UsePixel           = cms.bool(True),
    MeVperADCStrip     = cms.double(3.61e-06*250),
    MeVperADCPixel     = cms.double(3.61e-06),

    MisCalib_Mean      = cms.untracked.double(1.0),
    MisCalib_Sigma     = cms.untracked.double(0.00),

    UseCalibration  = cms.bool(True),
    calibrationPath = cms.string("file:Gains.root"),
    ShapeTest          = cms.bool(False),

    MaxNrStrips        = cms.untracked.uint32(255)
)

dedxSmi = dedxProd.clone()
dedxSmi.Formula = cms.untracked.uint32(2)

dedxASmi = dedxProd.clone()
dedxASmi.Formula = cms.untracked.uint32(3)


dedxSTProd                  = dedxProd.clone()
dedxSTProd.ShapeTest        = cms.bool(True)

dedxSTSmi                   = dedxSmi.clone()
dedxSTSmi.ShapeTest         = cms.bool(True)

dedxSTASmi                  = dedxASmi.clone()
dedxSTASmi.ShapeTest        = cms.bool(True)



####################################################################################
#   MUON TIMING
####################################################################################

from RecoMuon.MuonIdentification.muonTiming_cfi import *
muontiming.MuonCollection = cms.InputTag("muons")

####################################################################################
#   HSCParticle Producer
####################################################################################

#ALL THIS IS NEEDED BY ECAL BETA CALCULATOR (TrackAssociator)
from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *
from TrackingTools.TrackAssociator.default_cfi import * 
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
from Geometry.CMSCommonData.cmsIdealGeometryXML_cfi import *
from Geometry.CaloEventSetup.CaloGeometry_cff import *
from Geometry.CaloEventSetup.CaloTopology_cfi import *
from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi import *
from Geometry.MuonNumbering.muonNumberingInitialization_cfi import *
from Geometry.DTGeometry.dtGeometry_cfi import *
from Geometry.RPCGeometry.rpcGeometry_cfi import *
from Geometry.CSCGeometry.cscGeometry_cfi import *
from Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi import *


from SUSYBSMAnalysis.HSCP.HSCPSelections_cff import *

HSCParticleProducer = cms.EDProducer("HSCParticleProducer",
   TrackAssociatorParameterBlock, #Needed for ECAL/Track Matching

   #DOES THE PRODUCER ACT AS AN EDFILTER?
   filter = cms.bool(True),

   #WHAT (BETA) INFORMATION TO COMPUTE
   useBetaFromTk      = cms.bool(True),
   useBetaFromMuon    = cms.bool(True),
   useBetaFromRpc     = cms.bool(True),
   useBetaFromEcal    = cms.bool(False),

   #TAG OF THE REQUIRED INPUT COLLECTION (ONLY ACTIVATED CALCULATOR)
   tracks             = cms.InputTag("TrackRefitter"),
   muons              = cms.InputTag("muons"),
   dedxEstimator1     = cms.InputTag("dedxCNPHarm2"),
   dedxEstimator2     = cms.InputTag("dedxCNPTru40"),
   dedxEstimator3     = cms.InputTag("dedxCNPMed"),
   dedxEstimator4     = cms.InputTag("dedxSTCNPHarm2"),
   dedxEstimator5     = cms.InputTag("dedxSTCNPTru40"),
   dedxEstimator6     = cms.InputTag("dedxSTCNPMed"),
   dedxDiscriminator1 = cms.InputTag("dedxProd"),
   dedxDiscriminator2 = cms.InputTag("dedxSmi"),
   dedxDiscriminator3 = cms.InputTag("dedxASmi"),
   dedxDiscriminator4 = cms.InputTag("dedxSTProd"),
   dedxDiscriminator5 = cms.InputTag("dedxSTSmi"),
   dedxDiscriminator6 = cms.InputTag("dedxSTASmi"),
   muontimingDt       = cms.InputTag("muontiming:dt"),
   muontimingCsc      = cms.InputTag("muontiming:csc"),
   muontimingCombined = cms.InputTag("muontiming:combined"),

   #TRACK SELECTION FOR THE HSCP SEED
   minMuP             = cms.double(5),
   minTkP             = cms.double(5),
   maxTkChi2          = cms.double(5),
   minTkHits          = cms.uint32(3),

   #MUON/TRACK MATCHING THRESHOLDS (ONLY IF NO MUON INNER TRACK)
   minDR              = cms.double(0.1),
   maxInvPtDiff       = cms.double(0.005),

   #SELECTION ON THE PRODUCED HSCP CANDIDATES (WILL STORE ONLY INTERESTING CANDIDATES)
   SelectionParameters = cms.VPSet(
      HSCPSelectionDefault,
   ),
)


####################################################################################
#   HSCParticle Selector  (Just an Example of what we can do)
####################################################################################

HSCParticleSelector = cms.EDFilter("HSCParticleSelector",
   source = cms.InputTag("HSCParticleProducer"),
   filter = cms.bool(True),

   SelectionParameters = cms.VPSet(
      HSCPSelectionHighdEdx, #THE OR OF THE TWO SELECTION WILL BE APPLIED
      HSCPSelectionHighTOF,
   ),
)

####################################################################################
#   HSCP Candidate Sequence
####################################################################################

HSCParticleProducerSeq = cms.Sequence(offlineBeamSpot + TrackRefitter + dedxCNPHarm2 + dedxCNPTru40 + dedxCNPMed + dedxSTCNPHarm2 + dedxSTCNPTru40 + dedxSTCNPMed + dedxProd + dedxSmi + dedxASmi + dedxSTProd + dedxSTSmi + dedxSTASmi + muontiming + HSCParticleProducer)


