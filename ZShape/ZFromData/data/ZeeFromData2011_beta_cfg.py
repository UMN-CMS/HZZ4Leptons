import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("CondCore.DBCommon.CondDBSetup_cfi") # Alex removed
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi") # Alex removed
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi") # Alex removed
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff") # Alex added
#process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
#process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
#process.GlobalTag.globaltag = 'GR10_P_V10::All'
#process.GlobalTag.globaltag = 'GR_R_39X_V4::All'
process.GlobalTag.globaltag = 'GR_R_42_V12::All'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff") # Alex added
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi") # Alex added


process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring( "file:/local/cms/phedex/store/data/Run2011A/DoubleElectron/RECO/05Jul2011ReReco-HF/0000/0004FD54-68AA-E011-859D-003048678D52.root",)
)

#process.load("ZShape.ZFromData.GoodLumis_cfi")

#from ZShape.ZFromData.GoodLumis_cfi import *
#from ZShape.ZFromData.GoodLumis_July16ReReco_cfi import *
#from ZShape.ZFromData.GoodLumis_June14thReReco_cfi import *

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#    GoodLumis1+GoodLumis2+
#    GoodLumis3+
#    GoodLumis4+
#    GoodLumis5
#)

#print(process.source.lumisToProcess)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.TFileService = cms.Service("TFileService", 
        fileName = cms.string('ZFromData_output.root')
        )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        fileMode = cms.untracked.string('NOMERGE')
        )

process.TimerService = cms.Service("TimerService",
        useCPUtime = cms.untracked.bool(True)
        )

process.load("ZShape.ZFromData.ZFromDataElectrons_cfi") 

#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Photon15_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Photon15_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Photon15_Cleaned_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Photon15_Cleaned_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_L1R","","HLT")
#process.theHLT.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_CaloEleId_L1R","","HLT")
#process.theHLTGsf.hltTag = cms.untracked.InputTag("HLT_Ele15_SW_CaloEleId_L1R","","HLT")

process.vertFilter0 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(0), max_vert = cms.int32(0),)
process.vertFilter1 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(1), max_vert = cms.int32(1),)
process.vertFilter2 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(2), max_vert = cms.int32(2),)
process.vertFilter3 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(3), max_vert = cms.int32(3),)
process.vertFilter4 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(4), max_vert = cms.int32(4),)
process.vertFilter5 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(5), max_vert = cms.int32(5),)
process.vertFilter6 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(6), max_vert = cms.int32(6),)
process.vertFilter7 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(7), max_vert = cms.int32(7),)
process.vertFilter8 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(8), max_vert = cms.int32(8),)
process.vertFilter9 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(9), max_vert = cms.int32(9),)
process.vertFilter10 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(10), max_vert = cms.int32(10),)
process.vertFilter11 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(11), max_vert = cms.int32(11),)
process.vertFilter12 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(12), max_vert = cms.int32(12),)
process.vertFilter13 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(13), max_vert = cms.int32(13),)
process.vertFilter14 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(14), max_vert = cms.int32(14),)
process.vertFilter15 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(15), max_vert = cms.int32(15),)
process.vertFilter16 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(16), max_vert = cms.int32(16),)
process.vertFilter17 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(17), max_vert = cms.int32(17),)
process.vertFilter18 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(18), max_vert = cms.int32(18),)
process.vertFilter19 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(19), max_vert = cms.int32(19),)
process.vertFilter20 = cms.EDFilter("CutOnVerts", min_vert = cms.int32(20), max_vert = cms.int32(20),)

process.ZFromData = cms.EDAnalyzer("ZFromData",
        quiet = cms.untracked.bool(False),
        src = cms.untracked.InputTag('generator'),
        outHistogramsFile = cms.untracked.string('base_Jan23rd_1.root'),
        ExtraFromDataHistos = cms.untracked.bool(False),
        doMC = cms.untracked.bool(False),                               
        WeightsFile = cms.untracked.string('none'),                               

        dRMatchCut = cms.untracked.double(0.3),
        dPtMatchCut = cms.untracked.double(0.6),

        Effs = cms.untracked.VPSet(
            cms.PSet(
                name = cms.untracked.FileInPath('ZShape/EffAcc/data/Supercluster-Eta.txt'),
                effFile = cms.untracked.string('Supercluster-Eta.txt')
                #  ),
                #  cms.PSet(
                #         name = cms.untracked.string('GsfTrack-EtaDet'),
                #         effFile = cms.untracked.string('GsfTrack-EtaDet.txt')
                #  ),
                #  cms.PSet(
                #         name = cms.untracked.string('Iso-Pt'),
                #         effFile = cms.untracked.string('Iso-Pt.txt')
                #  ),
                #  cms.PSet(
                #         name = cms.untracked.string('ElectronId-EtaDet'),
                #         effFile = cms.untracked.string('ElectronId-EtaDet.txt')
                #  ),
                #  cms.PSet(
                #         name = cms.untracked.string('HLT-EtaDet'),
                #         effFile = cms.untracked.string('HLT-EtaDet.txt')
                )),

            #TagProbeProducer = cms.untracked.InputTag('tpMapGsfElectrons'),
            #TagProbeProducer = cms.untracked.InputTag('tpMapSuperClusters'),
            #TagProbeProducer = cms.untracked.InputTag('tpMapGsfAndHF'),
            #TagProbeProducer = cms.untracked.InputTag('tpMapWP95AndHF'),
            #TagProbeProducer = cms.untracked.InputTag('tpMapWP80AndSupers'), #Broken
            TagProbeProducer = cms.untracked.InputTag('tpMapWP90AndSupers'), #Broken
            #TagProbeProducer = cms.untracked.InputTag('tpMapWP80AndHF'),
            GsfProducer = cms.untracked.InputTag('theGsfElectrons'),

            CutNames          = cms.untracked.vstring("Supercluster-Eta", "GsfTrack-EtaDet",   "Iso-Pt",        "ElectronId-EtaDet", "HLT-EtaDet", "HFElectronId-EtaDet", "HFSuperCluster-Et","HFTightElectronId-EtaDet","EID95","ISO95","EID90","ISO90","EID85","ISO85","EID80","ISO80","EID70","ISO70","EID60","ISO60","HLT-GSF","ISO80Only","ISO80Conv","EID80Only","EID80Conv","WP95","WP90","WP85","WP80","NTLooseElectronId-EtaDet","NTTightElectronId-EtaDet" ),
            allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHFSuperClusters"), cms.InputTag("theHFSuperClusters"),cms.InputTag("theHFSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("theGsfElectrons"),cms.InputTag("NTElecLoose"),cms.InputTag("NTElecLoose")),
            passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHLT"), cms.InputTag("HFElectronID"), cms.InputTag("theHFSuperClusters"), cms.InputTag("HFElectronIDTight"), cms.InputTag("ElectronID95"), cms.InputTag("Iso95"), cms.InputTag("ElectronID90"), cms.InputTag("Iso90"), cms.InputTag("ElectronID85"), cms.InputTag("Iso85"), cms.InputTag("ElectronID80"), cms.InputTag("Iso80"), cms.InputTag("ElectronID70"), cms.InputTag("Iso70"), cms.InputTag("ElectronID60"), cms.InputTag("Iso60"), cms.InputTag("theHLTGsf"), cms.InputTag("Iso80Only"), cms.InputTag("Iso80WConv"), cms.InputTag("ElectronID80Only"), cms.InputTag("ElectronID80WConv"), cms.InputTag("ElectronID95"),cms.InputTag("ElectronID90"),cms.InputTag("ElectronID85"),cms.InputTag("ElectronID80"),cms.InputTag("NTElecLoose"),cms.InputTag("NTElecTight")),
            ExactMatch        = cms.untracked.vint32(1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),

            writeHistoBeforeEndJob = cms.untracked.bool(False),

            #Here is where the bulk of everything goes, into the Z definitions
            ZDefs = cms.untracked.VPSet( 
                #                cms.PSet(
                #                    name = cms.untracked.string('Maximal'),
                #                    e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                #                        "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
                #                    e2 = cms.untracked.vstring("ACC(ANY)","PT10"),						 
                #                    Z = cms.untracked.vstring(),
                #                    ptorder = cms.untracked.bool (False)
                #                    ),
                #                cms.PSet(
                #                    name = cms.untracked.string('AllInRange'),
                #                    Z = cms.untracked.vstring('m(60,120)'),
                #                    e1 = cms.untracked.vstring('PT0'),
                #                    e2 = cms.untracked.vstring('PT0')
                #                    ),
                #                cms.PSet(
                #                    name = cms.untracked.string('ECAL85-ECAL85'),
                #                    e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                #                        "GsfTrack-EtaDet","ISO85","EID85", "HLT-GSF"),
                #                    e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                #                        "GsfTrack-EtaDet","ISO85","EID85"),		
                #                    Z = cms.untracked.vstring('m(60,120)'),
                #                    ptorder = cms.untracked.bool (True)
                #                    ),
                #                cms.PSet(
                #                    name = cms.untracked.string('ECAL80-HF'),
                #                    e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                #                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
                #                    e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
                #                    Z = cms.untracked.vstring('m(60,120)'),
                #                    ptorder = cms.untracked.bool (False)
                #                    ),
                #                cms.PSet(
                #                    name = cms.untracked.string('ECAL80-HFTight'),
                #                    e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
                #                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
                #                    e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
                #                    Z = cms.untracked.vstring('m(60,120)'),
                #                    ptorder = cms.untracked.bool (False)
                #                    ),
                cms.PSet(
                    name = cms.untracked.string('ECAL_Efficiency_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","Supercluster-Eta", "GsfTrack-EtaDet","WP80"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('HF_Efficiency_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(HF)","PT20","HFSuperCluster-Et", "HFElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('Trackless_Efficiency_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","Supercluster-Eta","NTTightElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('ECAL_Background_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","!WP95"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('HF_Background_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(HF)","PT20","!HFElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('Trackless_Background_Chain'),
                    e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                    e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","!NTLooseElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('ECAL_Background_90_Chain'),
                    e1 = cms.untracked.vstring( "WP90","WP90","WP90","WP90"),
                    e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","!GsfTrack-EtaDet","!WP95"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('HF_Background_90_Chain'),
                    e1 = cms.untracked.vstring( "WP90","WP90","WP90"),
                    e2 = cms.untracked.vstring("ACC(HF)","PT20","!HFElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                cms.PSet(
                    name = cms.untracked.string('Trackless_Background_90_Chain'),
                    e1 = cms.untracked.vstring( "WP90","WP90","WP90"),
                    e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","!NTLooseElectronId-EtaDet"),		
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
                #                cms.PSet(
                #                    name = cms.untracked.string('ECAL70-HF'),
                #                    e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
            #                        "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
#                    e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
#                    Z = cms.untracked.vstring('m(60,120)'),
#                    ptorder = cms.untracked.bool (False)
#                    ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL70-HFTight'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),		
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL70-ECAL70'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","ISO70","EID70", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","ISO70","EID70"),		
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (True)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-HFTight-15'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT15",
#                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT15","HFTightElectronId-EtaDet"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-HFTight-17'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT17",
#                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-HFTight-20'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-HFTight-25'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
#                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-HFTight-30'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
#                            "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80EB-ECAL80EE'),
#                        e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80EE-ECAL80EE'),
#                        e1 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (True)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('MUON-MUON'),
#                        e1 = cms.untracked.vstring("ACC(MUON)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(MUON)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP95"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (True)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80EB-ECAL80EB'),
#                        e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80"),
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (True)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-NTLoose'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTLooseElectronId-EtaDet"),		
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        ),
#                cms.PSet(
#                        name = cms.untracked.string('ECAL80-NTTight'),
#                        e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
#                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
#                        e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTTightElectronId-EtaDet"),		
#                        Z = cms.untracked.vstring('m(60,120)'),
#                        ptorder = cms.untracked.bool (False)
#                        )
                )
            )

process.ZFromData0 = process.ZFromData.clone()
process.ZFromData1 = process.ZFromData.clone()
process.ZFromData2 = process.ZFromData.clone()
process.ZFromData3 = process.ZFromData.clone()
process.ZFromData4 = process.ZFromData.clone()
process.ZFromData5 = process.ZFromData.clone()
process.ZFromData6 = process.ZFromData.clone()
process.ZFromData7 = process.ZFromData.clone()
process.ZFromData8 = process.ZFromData.clone()
process.ZFromData9 = process.ZFromData.clone()
process.ZFromData10 = process.ZFromData.clone()
process.ZFromData11 = process.ZFromData.clone()
process.ZFromData12 = process.ZFromData.clone()
process.ZFromData13 = process.ZFromData.clone()
process.ZFromData14 = process.ZFromData.clone()
process.ZFromData15 = process.ZFromData.clone()
process.ZFromData16 = process.ZFromData.clone()
process.ZFromData17 = process.ZFromData.clone()
process.ZFromData18 = process.ZFromData.clone()
process.ZFromData19 = process.ZFromData.clone()
process.ZFromData20 = process.ZFromData.clone()

#process.load("ISpy/Analyzers/ISpy_Producer_cff")
#process.p2 = cms.Path(process.iSpy_sequence)

#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.vertFilter0 * (process.lepton_cands + process.ZFromData0))
process.p2 = cms.Path(process.vertFilter1 * (process.lepton_cands + process.ZFromData1))
process.p3 = cms.Path(process.vertFilter2 * (process.lepton_cands + process.ZFromData2))
process.p4 = cms.Path(process.vertFilter3 * (process.lepton_cands + process.ZFromData3))
process.p5 = cms.Path(process.vertFilter4 * (process.lepton_cands + process.ZFromData4))
process.p6 = cms.Path(process.vertFilter5 * (process.lepton_cands + process.ZFromData5))
process.p7 = cms.Path(process.vertFilter6 * (process.lepton_cands + process.ZFromData6))
process.p8 = cms.Path(process.vertFilter7 * (process.lepton_cands + process.ZFromData7))
process.p9 = cms.Path(process.vertFilter8 * (process.lepton_cands + process.ZFromData8))
process.p10 = cms.Path(process.vertFilter9 * (process.lepton_cands + process.ZFromData9))
process.p11 = cms.Path(process.vertFilter10 * (process.lepton_cands + process.ZFromData10))
process.p12 = cms.Path(process.vertFilter11 * (process.lepton_cands + process.ZFromData11))
process.p13 = cms.Path(process.vertFilter12 * (process.lepton_cands + process.ZFromData12))
process.p14 = cms.Path(process.vertFilter13 * (process.lepton_cands + process.ZFromData13))
process.p15 = cms.Path(process.vertFilter14 * (process.lepton_cands + process.ZFromData14))
process.p16 = cms.Path(process.vertFilter15 * (process.lepton_cands + process.ZFromData15))
process.p17 = cms.Path(process.vertFilter16 * (process.lepton_cands + process.ZFromData16))
process.p18 = cms.Path(process.vertFilter17 * (process.lepton_cands + process.ZFromData17))
process.p19 = cms.Path(process.vertFilter18 * (process.lepton_cands + process.ZFromData18))
process.p20 = cms.Path(process.vertFilter19 * (process.lepton_cands + process.ZFromData19))
process.p21 = cms.Path(process.vertFilter20 * (process.lepton_cands + process.ZFromData20))
