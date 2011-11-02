import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_1.root',
        #'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_2.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_3.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_4.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_5.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_6.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_7.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_8.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_9.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_10.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_11.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_12.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_13.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_14.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_15.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_16.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_17.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_18.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_19.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_20.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_21.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_22.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_23.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_24.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_25.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_26.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_27.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_28.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_29.root',
	#   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO/Zee_10TeV_HLT_RECO_30.root'
           'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_001.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_002.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_003.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_004.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_005.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_006.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_007.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_008.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_009.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_010.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_011.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_012.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_013.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_014.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_015.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_016.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_017.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_018.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_019.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_020.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_021.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_022.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_023.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_024.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_025.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_026.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_027.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_028.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_029.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_030.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_031.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_032.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_033.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_034.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_035.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_036.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_037.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_038.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_039.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_040.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_041.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_042.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_043.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_044.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_045.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_046.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_047.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_048.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_049.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_050.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_051.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_052.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_053.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_054.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_055.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_056.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_057.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_058.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_059.root',
	   #'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_060.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_061.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_062.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_063.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_064.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_065.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_066.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_067.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_068.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_069.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_070.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_071.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_072.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_073.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_074.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_075.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_076.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_077.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_078.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_079.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_080.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_080.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_081.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_082.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_083.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_084.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_085.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_086.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_087.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_088.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_089.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_090.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_091.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_092.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_093.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_094.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_095.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_096.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_097.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_098.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_099.root',
	   'file:/data/whybee0c/user/jmmans/Zee_10TeV_HLT_RECO_217/Zee_10TeV_HLT_RECO_100.root'
	)
)


process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(123200)
   #input = cms.untracked.int32(1000)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('histoJeremy09_only1hlt_TRIAL.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.TimerService = cms.Service("TimerService",
    useCPUtime = cms.untracked.bool(True)
)

process.load("ZShape.ZFromData.ZFromDataElectrons_cfi") 

process.ZFromData = cms.EDAnalyzer("ZFromData",
    quiet = cms.untracked.bool(False),
    outHistogramsFile = cms.untracked.string('base_Jan23rd_1.root'),
    ExtraFromDataHistos = cms.untracked.bool(True),

	dRMatchCut = cms.untracked.double(0.3),
	dPtMatchCut = cms.untracked.double(0.3),

	Effs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.string('Supercluster-Eta'),
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
	TagProbeProducer = cms.untracked.InputTag('tpMapGsfAndHF'),
	
	CutNames          = cms.untracked.vstring("Supercluster-Eta", "GsfTrack-EtaDet",   "Iso-Pt",        "ElectronId-EtaDet", "HLT-EtaDet", "HFElectronId-EtaDet", "EEEdge-Eta"),
	allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHFSuperClusters"),cms.InputTag("theSuperClusters")),
	passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHLT"), cms.InputTag("HFElectronID"), cms.InputTag("theEEHFGapSuperClusters")),
    ExactMatch        = cms.untracked.vint32(0,0,0,0,0),
	
	writeHistoBeforeEndJob = cms.untracked.bool(False),

	#Here is where the bulk of everything goes, into the Z definitions
	ZDefs = cms.untracked.VPSet(
	  cms.PSet(
	         name = cms.untracked.string('Golden-EB-EB'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         #"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),	
					 "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),						 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-EB-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet","HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),						 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Maximal'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ANY)","PT10"),						 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),		 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-VeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta"),		 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-AlmostVeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20"),			 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(HF)","PT10","PT20","HFElectronId-EtaDet"),		
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-EEHighEta'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20","EEEdge-Eta"),			 
			 z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (False)
	  )
	)
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(   #'drop *', 
        'keep *_*_*_Zefficiency'),
    fileName = cms.untracked.string('OutJeremy09_1.root')
)

#process.p1 = cms.Path(process.hfEMClusteringSequence * process.lepton_cands+process.TPEdm)
process.p1 = cms.Path(process.lepton_cands+process.ZFromData)
process.the_end = cms.EndPath(process.outpath)

