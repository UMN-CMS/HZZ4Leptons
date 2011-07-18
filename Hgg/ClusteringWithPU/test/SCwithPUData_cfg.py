#inspired to/by:
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatExamples/test/patTuple_42x_jec_cfg.py?revision=1.2&view=markup&sortby=date
# http://cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bicocca/PhysicsTools/MiBiCommonPAT/python/makeMiBiCommonNT_cff.py?revision=1.37&view=markup

import FWCore.ParameterSet.Config as cms

process = cms.Process("ZeeForPU")

from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *


# Setup the process
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V14::All'


# Source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/data/franzoni/data/Run2011A_DoubleElectron_AOD_PromptReco-v4_000_166_946_CE9FBCFF-4B98-E011-A6C3-003048F11C58.root')
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))


# Output
process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring(),
    fileName = cms.untracked.string('file:pippo.root'),
    )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("zeePU-histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


# cleaning: scraping filter
process.scrapingFilter = cms.EDFilter("FilterOutScraping",
                                      applyfilter = cms.untracked.bool(True),
                                      debugOn = cms.untracked.bool(False),
                                      numtrack = cms.untracked.uint32(10),
                                      thresh = cms.untracked.double(0.25)
                                      )



# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(3.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )


#------------------
#Load PAT sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.tools.pfTools")
postfix = "PFlow"

# THis is NOT MC => remove matching
removeMCMatching(process, ['All'])


# bugfix for DATA Run2011 (begin)
removeSpecificPATObjects( process, ['Taus'] )
process.patDefaultSequence.remove( process.patTaus )
# bugfix for DATA Run2011 (end)


# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. It is currently
# not possible to run PF2PAT+PAT and standart PAT at the same time
#from PhysicsTools.PatAlgos.tools.pfTools import *
# postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpPFlow.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False
 

#    
# the HCAL Noise Filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

#
# build met objects
#process.layer1METs = process.patMETs.clone(
#    metSource = cms.InputTag("tcMet","",""),
#    addTrigMatch = cms.bool(False),
#    addMuonCorrections = cms.bool(False),
#    addGenMET = cms.bool(False),
#    )
### specify here what you want to have on the plots! <===== MET THAT YOU WANT ON THE PLOTS  %%%%%%%
#myDesiredMetCollection = 'layer1METs'
## modify the sequence of the MET creation:
#process.makePatMETs = cms.Sequence(process.layer1METs)
##



## modify the final pat sequence: keep only electrons + METS (muons are needed for met corrections)
process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")

##
process.patElectrons.isoDeposits = cms.PSet()

process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),    
    )
##
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
##
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons)
# process.makePatMuons may be needed depending on how you calculate the MET
process.makePatCandidates = cms.Sequence(process.makePatElectrons
                                         # +process.makePatMETs
                                         )
process.patMyDefaultSequence = cms.Sequence(process.makePatCandidates)


process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")

# this is needed to compute islation (temporary?) 
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")




##
## setting up filters
#

process.ElectronID60 = cms.EDFilter("CandViewSelector",
                                    src = cms.InputTag("patElectrons"),
                                    cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta) >' + str(20) + ') && ( electronID("simpleEleId60cIso") == 7 )')
                                    )
process.ElectronID60Cut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronID60"),
                                       minNumber = cms.uint32(2)
                                       )

process.ElectronID80 = cms.EDFilter("CandViewSelector",
                                    src = cms.InputTag("patElectrons"),
                                    cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta) >' + str(20) + ') && ( electronID("simpleEleId80cIso") == 7 )')
                                    )
process.ElectronID80Cut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronID80"),
                                       minNumber = cms.uint32(2)
                                       )

process.ElectronID90 = cms.EDFilter("CandViewSelector",
                                    src = cms.InputTag("patElectrons"),
                                    cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta) >' + str(20) + ') && ( electronID("simpleEleId90cIso") == 7 )')
                                    )
process.ElectronID90Cut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronID90"),
                                       minNumber = cms.uint32(2)
                                )

process.ElectronID95 = cms.EDFilter("CandViewSelector",
                                    src = cms.InputTag("patElectrons"),
                                    cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta) >' + str(20) + ') && ( electronID("simpleEleId95cIso") == 7 )')
                                    )
process.ElectronID95Cut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronID95"),
                                       minNumber = cms.uint32(2)
                                )



# this is the ntuple producer
process.load("CalibCalorimetry.EcalTiming.ecalTimeEleTree_cfi")
process.ecalTimeEleTree.fileName = 'EcalTimeTree'
process.ecalTimeEleTree.muonCollection = cms.InputTag("muons")
process.ecalTimeEleTree.runNum = 999999

# this is where Il 
#process.load("ECALTime.SCwithPUData.scwithpudata_cfi")
from Hgg.ClusteringWithPU.scwithpudata_cfi import *
process.base          = SCwithPUData.clone()
process.base.myName   = cms.string('base')
process.base.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base.xi       = cms.double(0.01)

process.base_0          = SCwithPUData.clone()
process.base_0.myName   = cms.string('base-xi0.0')
process.base_0.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_0.xi       = cms.double(0.0)

process.base_0005          = SCwithPUData.clone()
process.base_0005.myName   = cms.string('base-xi0.005')
process.base_0005.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_0005.xi       = cms.double(0.005)

process.base_001          = SCwithPUData.clone()
process.base_001.myName   = cms.string('base-xi0.01')
process.base_001.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_001.xi       = cms.double(0.01)

process.base_0015          = SCwithPUData.clone()
process.base_0015.myName   = cms.string('base-xi0.015')
process.base_0015.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_0015.xi       = cms.double(0.015)

process.base_002          = SCwithPUData.clone()
process.base_002.myName   = cms.string('base-xi0.02')
process.base_002.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_002.xi       = cms.double(0.02)

process.base_0025          = SCwithPUData.clone()
process.base_0025.myName   = cms.string('base-xi0.02')
process.base_0025.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_0025.xi       = cms.double(0.025)

process.base_003          = SCwithPUData.clone()
process.base_003.myName   = cms.string('base-xi0.03')
process.base_003.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_003.xi       = cms.double(0.03)

process.base_0035          = SCwithPUData.clone()
process.base_0035.myName   = cms.string('base-xi0.03')
process.base_0035.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_0035.xi       = cms.double(0.035)

process.base_004          = SCwithPUData.clone()
process.base_004.myName   = cms.string('base-xi0.04')
process.base_004.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_004.xi       = cms.double(0.04)

process.base_005          = SCwithPUData.clone()
process.base_005.myName   = cms.string('base-xi0.05')
process.base_005.eleWorkingPoint   = cms.string('simpleEleId95relIso')
process.base_005.xi       = cms.double(0.05)


process.plotsWP90         = SCwithPUData.clone()
process.plotsWP90.myName  = cms.string('WP90')
process.plotsWP90.eleWorkingPoint   = cms.string('simpleEleId90relIso')


process.plotsWP80             = SCwithPUData.clone()
process.plotsWP80.myName      = cms.string('plosWP80')
process.plotsWP80.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80.xi          = cms.double(0.00)

process.plotsWP80_0005          = SCwithPUData.clone()
process.plotsWP80_0005.myName   = cms.string('plotsWP80-xi0.005')
process.plotsWP80_0005.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_0005.xi       = cms.double(0.005)

process.plotsWP80_001          = SCwithPUData.clone()
process.plotsWP80_001.myName   = cms.string('plotsWP80-xi0.01')
process.plotsWP80_001.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_001.xi       = cms.double(0.01)

process.plotsWP80_0015          = SCwithPUData.clone()
process.plotsWP80_0015.myName   = cms.string('plotsWP80-xi0.015')
process.plotsWP80_0015.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_0015.xi       = cms.double(0.015)

process.plotsWP80_002          = SCwithPUData.clone()
process.plotsWP80_002.myName   = cms.string('plotsWP80-xi0.02')
process.plotsWP80_002.eleWorkingPoint = cms.string('simpleEleId80relIso')
process.plotsWP80_002.xi       = cms.double(0.02)

process.plotsWP80_0025          = SCwithPUData.clone()
process.plotsWP80_0025.myName   = cms.string('plotsWP80-xi0.025')
process.plotsWP80_0025.eleWorkingPoint = cms.string('simpleEleId80relIso')
process.plotsWP80_0025.xi       = cms.double(0.025)

process.plotsWP80_003          = SCwithPUData.clone()
process.plotsWP80_003.myName   = cms.string('plotsWP80-xi0.03')
process.plotsWP80_003.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_003.xi       = cms.double(0.03)

process.plotsWP80_0035          = SCwithPUData.clone()
process.plotsWP80_0035.myName   = cms.string('plotsWP80-xi0.035')
process.plotsWP80_0035.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_0035.xi       = cms.double(0.035)

process.plotsWP80_004          = SCwithPUData.clone()
process.plotsWP80_004.myName   = cms.string('plotsWP80-xi0.04')
process.plotsWP80_004.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_004.xi       = cms.double(0.04)

process.plotsWP80_005          = SCwithPUData.clone()
process.plotsWP80_005.myName   = cms.string('plotsWP80-xi0.05')
process.plotsWP80_005.eleWorkingPoint   = cms.string('simpleEleId80relIso')
process.plotsWP80_005.xi       = cms.double(0.05)


process.plotsWP60         = SCwithPUData.clone()
process.plotsWP60.myName  = cms.string('WP60')
process.plotsWP60.eleWorkingPoint   = cms.string('simpleEleId60relIso')



# this is to check event content
process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")


# sequence for electrons above a certain Pt and passing a certain iso
# things standing, working point needs to be re-checked inside the analyzer too
process.WP60Seq = cms.Sequence(process.ElectronID60 * process.ElectronID60Cut * process.plotsWP60)
process.WP80Seq = cms.Sequence(process.ElectronID80 * process.ElectronID80Cut
                             * process.plotsWP80 *  process.plotsWP80_0005 *
                               process.plotsWP80_001 * process.plotsWP80_0015 * process.plotsWP80_002 * process.plotsWP80_0025
                             * process.plotsWP80_003 * process.plotsWP80_0035 * process.plotsWP80_004 * process.plotsWP80_005)
process.WP90Seq = cms.Sequence(process.ElectronID90 * process.ElectronID90Cut * process.plotsWP90)
process.baseSeq  = cms.Sequence(process.ElectronID95 * process.ElectronID95Cut * process.base * process.base_0
                                * process.base_0005  * process.base_001  * process.base_0015 * process.base_002 * process.base_0025
                                * process.base_003   * process.base_0035 * process.base_004
                                * process.base_005)

process.PrepareSeq =  cms.Sequence(
    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *
    getattr(process,"patPF2PATSequence"+postfix)  *
    process.patMyDefaultSequence 
    #process.makePatMETs
    )


#################################################### this is what runs the show... ########
process.eleZforPU_WP60Seq = cms.Path( process.PrepareSeq * process.WP60Seq )
process.eleZforPU_WP80Seq = cms.Path( process.PrepareSeq * process.WP80Seq )
process.eleZforPU_WP90Seq = cms.Path( process.PrepareSeq * process.WP90Seq )
process.eleZforPU_baseSeq  = cms.Path( process.PrepareSeq * process.baseSeq )
