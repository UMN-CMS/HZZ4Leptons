# Auto generated configuration file
# using: 
# Revision: 1.303.2.3 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: ZEE_7TeV_cfi.py -s GEN,FASTSIM,HLT:GRun --pileup=HighLumiPileUp --geometry DB --conditions=auto:mc --eventcontent=RECOSIM --datatier GEN-SIM-DIGI-RECO -n 10 --no_exec

# Mon Jul 11 22:59:19 CEST 2011
# command: cmsDriver.py ZEE_7TeV_cfi.py -s GEN,FASTSIM,HLT:GRun --pileup=HighLumiPileUp --geometry DB --conditions=auto:mc --eventcontent=RECOSIM --datatier GEN-SIM-DIGI-RECO -n 10 --no_exec 

# this python is meant to generate Zee from fast sim
# run the 'historical' pu analysis which uses truth
# run also the other analysis which is geared towards real data for a cross check

import FWCore.ParameterSet.Config as cms

#process = cms.Process('HLT')
process = cms.Process('gf-fastSimForPAT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_HighLumiPileUp_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('FastSimulation.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# number of pile up events
process.famosPileUp.PileUpSimulator.averageNumber = 0

# import of standard&useful configurations
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(4000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.3 $'),
    annotation = cms.untracked.string('ZEE_7TeV_cfi.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('ZEE_7TeV_cfi_py_GEN_FASTSIM_HLT_PU.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('GEN-SIM-DIGI-RECO')
#    ),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('generation_step')
#    )
#)

# tfile service, where both analyzes send their results :-)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('PU_scluster.root')
)
# Output
process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring(),
    fileName = cms.untracked.string('file:pippo.root'),
    )


# Additional output definition

# Other statements
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.simulation = cms.Sequence(process.simulationWithFamos)
process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)
process.Realistic7TeV2011CollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic7TeV2011CollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic7TeV2011CollisionVtxSmearingParameters
process.GlobalTag.globaltag = 'MC_42_V12::All'

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1      !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MSEL         = 11 ', 
            'MDME( 174,1) = 0    !Z decay into d dbar', 
            'MDME( 175,1) = 0    !Z decay into u ubar', 
            'MDME( 176,1) = 0    !Z decay into s sbar', 
            'MDME( 177,1) = 0    !Z decay into c cbar', 
            'MDME( 178,1) = 0    !Z decay into b bbar', 
            'MDME( 179,1) = 0    !Z decay into t tbar', 
            'MDME( 182,1) = 1    !Z decay into e- e+', 
            'MDME( 183,1) = 0    !Z decay into nu_e nu_ebar', 
            'MDME( 184,1) = 0    !Z decay into mu- mu+', 
            'MDME( 185,1) = 0    !Z decay into nu_mu nu_mubar', 
            'MDME( 186,1) = 0    !Z decay into tau- tau+', 
            'MDME( 187,1) = 0    !Z decay into nu_tau nu_taubar', 
            'CKIN( 1)     = 40.  !(D=2. GeV)', 
            'CKIN( 2)     = -1.  !(D=-1. GeV)'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


################################# take care of running pat here ############################
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

#------------------
#Load PAT sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.tools.pfTools")
#postfix = "PFlow"
#
## THis is NOT MC => remove matching
removeMCMatching(process, ['All'])
#
#
## bugfix for DATA Run2011 (begin)
removeSpecificPATObjects( process, ['Taus'] )
process.patDefaultSequence.remove( process.patTaus )


### modify the final pat sequence: keep only electrons + METS (muons are needed for met corrections)
#process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#
###
process.patElectrons.isoDeposits = cms.PSet()
#
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
###
#process.patElectrons.addGenMatch = cms.bool(False)
#process.patElectrons.embedGenMatch = cms.bool(False)
###
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs *
					process.patElectrons)
process.makePatCandidates = cms.Sequence( process.makePatElectrons   )
process.patMyDefaultSequence = cms.Sequence(process.makePatCandidates)
#
# stuff already loaded
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")


###################################################################################################
#                          my analysis modules
###################################################################################################
# this is where Il 
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

process.baseSeq  = cms.Sequence(  process.base_0 * process.base
                                * process.base_0005  * process.base_001  * process.base_0015 * process.base_002 * process.base_0025
                                * process.base_003   * process.base_0035 * process.base_004
                                * process.base_005)



process.scwithpuanalyzer = cms.EDAnalyzer('SCwithPUAnalysis')
process.scwithpuanalyzer.useRawEnergy = cms.bool(False)

process.analysisPath = cms.Sequence(
	#process.base  *
	process.baseSeq *
	process.scwithpuanalyzer
	)

# process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_genonly)
process.reconstruction = cms.Path(process.reconstructionWithFamos
				  * process.patMyDefaultSequence
				  )
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary
					    * process.analysisPath
					    )
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.reconstruction
			 #,process.RECOSIMoutput_step
			 ])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
