import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFROMDATA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( "")
)


process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(123200)
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
	TagProbeProducer = cms.untracked.InputTag('tpMapGsfAndHF'),
	GsfProducer = cms.untracked.InputTag('thePt20GsfElectrons'),
	
	CutNames          = cms.untracked.vstring("Supercluster-Eta", "GsfTrack-EtaDet",   "Iso-Pt",        "ElectronId-EtaDet", "HLT-EtaDet", "HFElectronId-EtaDet", "EEEdge-Eta", "HFSuperCluster-Et"),
	allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHFSuperClusters"),cms.InputTag("theSuperClusters"), cms.InputTag("theHFSuperClusters")),
	passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"),cms.InputTag("theGsfElectrons"),cms.InputTag("theIsolation"),cms.InputTag("theId"),cms.InputTag("theHLT"), cms.InputTag("HFElectronID"), cms.InputTag("theEEHFGapSuperClusters"), cms.InputTag("theHFSuperClusters")),
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
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-EB-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet","HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
										"GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),						 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Maximal'),
	         e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ANY)","PT10"),						 
			 Z = cms.untracked.vstring(),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),		 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-VeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta"),		 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-AlmostVeryLoose-ECAL'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20"),			 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-HF'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
			 Z = cms.untracked.vstring('m(70,110)','G(1,1)'),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-EEHighEta'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20","EEEdge-Eta"),			 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (False)
	  ),
	  cms.PSet(
	         name = cms.untracked.string('Tight-ECAL-Loose-EE'),
	         e1 = cms.untracked.vstring("ACC(ECAL)","PT10","Supercluster-Eta","PT20",
                                         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet", "HLT-EtaDet"),
			 e2 = cms.untracked.vstring("ACC(EE)","PT10","Supercluster-Eta","PT20",
				         "GsfTrack-EtaDet","Iso-Pt","ElectronId-EtaDet"),			 
			 Z = cms.untracked.vstring('m(70,110)','G(2,2)'),
		    ptorder = cms.untracked.bool (True)
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

