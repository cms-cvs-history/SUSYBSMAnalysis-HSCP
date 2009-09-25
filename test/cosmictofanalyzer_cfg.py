import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/FED0633E-F199-DE11-A788-001731AF67E1.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/FEB11B6B-2299-DE11-B8CB-0018F3D09600.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/FABF433B-CB98-DE11-89A8-0030486790A6.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/F6E966F6-3999-DE11-89F6-001731AF6A8D.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/F6C88565-2299-DE11-AA43-0018F3D09658.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/F66D5D48-2299-DE11-8FFD-0018F3D0962E.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/F078C87D-2299-DE11-9379-001731AF66EF.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/F072EC4B-2299-DE11-AC15-001A92971B94.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/D876079A-2299-DE11-A692-001A92971B5E.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/D237AD57-2299-DE11-A727-0018F3D0962E.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/D20C9B55-4399-DE11-A44E-001731EF61B4.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/CA4F1A5D-2299-DE11-B187-0018F3D09616.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/C0562B68-2299-DE11-A2E3-0018F3D095EE.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/BE8C4F6F-2299-DE11-9EA8-001731AF687F.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/B67CC57F-2299-DE11-9BEF-001A92971BDA.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/B620BB3D-CB98-DE11-95E5-00304867BFAE.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/B21D0D69-2299-DE11-9A44-001A92971BA0.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/AA4CD86C-2299-DE11-B7E9-0018F3D0966C.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/A864536D-2299-DE11-9CB4-001A92810A94.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/A0A01298-2299-DE11-AA9E-001A92971B88.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/9C9B25F8-889B-DE11-B3A6-0018F3D095EC.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/9A178D7A-2299-DE11-B41B-0018F3D09628.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/9664565B-4399-DE11-8FDB-001731AF66C1.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/9629BCF8-3999-DE11-A445-001731AF67B9.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/8CC9338A-2299-DE11-953A-001A92971BDA.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/8A59EFDE-CA98-DE11-A7CE-001BFCDBD160.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/8822C26A-2299-DE11-96E5-0018F3D095EA.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/86BCCC54-2299-DE11-933F-0018F3D09612.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/7EDB747C-2299-DE11-B4BD-0018F3D09658.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/7E60B048-2299-DE11-A701-0018F3D096BA.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/74CC465A-2299-DE11-A06A-001A92971BB8.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/748B3F48-2299-DE11-90AB-001A92971B08.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/74647679-2299-DE11-92FA-001731AF66B9.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/72A65A91-2299-DE11-B66E-0018F3D09628.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/70DE6D27-CB98-DE11-85DB-003048679084.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/7042F764-2299-DE11-BAFC-0018F3D09636.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/701D927B-2299-DE11-A15A-001A92971AA8.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/6EEF9E85-2299-DE11-AE83-0018F3D0969C.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/6CBB7184-2299-DE11-979A-001A9281173E.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/689BD74E-4399-DE11-B1D3-003048679076.root',
        '/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0022/6822B871-2299-DE11-AD2B-001A92810A94.root'
    )
)

process.tofAnalysis  = cms.EDAnalyzer('CosmicTOFAnalyzer',
 ByRun = cms.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)


process.p = cms.Path(process.tofAnalysis)
