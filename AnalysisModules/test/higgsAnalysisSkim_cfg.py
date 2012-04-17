import FWCore.ParameterSet.Config as cms

import os

#--- Data/MC switch ---#
isMC=True
isData=not isMC

#--- Special flag for 44x/Fall11 ---#
is44x=True

#--- Flags for data taking era ---#
isRun2011A = True
dataEra    = 20111 
pileupEra  = 20111 
if is44x:
    pileupEra = 20113
#--- Placeholder for 2011B variables
if not isRun2011A:
    dataEra   = 20112
    pileupEra = 20112
    if is44x:
        pileupEra = 20114

#--- Flags for nominal studies ---#
runAnalysis = True

process = cms.Process("HZZ");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/Fall11/GluGluToHToZZTo4L_M-120_7TeV-powheg-pythia6/AODSIM/PU_S6_START44_V9B-v1/0000/28D38082-FE2E-E111-8BE3-00215E221EC6.root')
    # fileNames=cms.untracked.vstring('input.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Global Tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    if (is44x): 
        process.GlobalTag.globaltag=cms.string('START44_V12::All')
    else: 
        process.GlobalTag.globaltag=cms.string('START42_V13::All')
else:
    print "===============> Running on DATA <===================="
    if (is44x):
        process.GlobalTag.globaltag = cms.string('GR_R_44_V13::All')
    else:
        process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

#--- Output module: 
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('candevents.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands = cms.untracked.vstring("keep *")
                               )
    
#------------------------------#
#--- Include Generic Tracks ---#
#------------------------------#
#--- Generic PAT tracks modules stolen from ElectroWeakAnalysis/Skimming/python ---#
process.patAODTrackCandsUnfiltered = cms.EDProducer("ConcreteChargedCandidateProducer",
    src          = cms.InputTag("generalTracks"),
    particleType = cms.string('mu+')   # to fix mass hypothesis
)
process.patAODTrackCands = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patAODTrackCandsUnfiltered"),
    cut = cms.string('pt > 10')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.allPatTracks = patGenericParticles.clone(
    src = cms.InputTag("patAODTrackCands")
)
from PhysicsTools.PatAlgos.selectionLayer1.trackSelector_cfi import *
process.patTracksPt10 = selectedPatTracks.clone(
    cut = 'pt > 10.'
)
process.patTrackSequence = cms.Sequence( 
        process.patAODTrackCandsUnfiltered *
        process.patAODTrackCands *
        process.allPatTracks *
        process.patTracksPt10
)

## --------------------- ##
## Define the basic path ##
## --------------------- ##

#--- Electron identification ---#
process.load("RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi")
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
process.load("RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi")
from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *

process.eIDRobustLoose = eidCutBasedExt.clone()
process.eIDRobustLoose.electronIDType = 'robust'
process.eIDRobustLoose.electronQuality = 'loose'

process.eIDRobustLooseV00 = eidCutBasedExt.clone()
process.eIDRobustLooseV00.electronIDType = 'robust'
process.eIDRobustLooseV00.electronQuality = 'loose'
process.eIDRobustLooseV00.electronVersion = 'V00'

process.eIDRobustTight = eidCutBasedExt.clone()
process.eIDRobustTight.electronIDType = 'robust'
process.eIDRobustTight.electronQuality = 'tight'

process.eIDRobustHighEnergy = eidCutBasedExt.clone()
process.eIDRobustHighEnergy.electronIDType = 'robust'
process.eIDRobustHighEnergy.electronQuality = 'highenergy'

process.eIDLoose = eidCutBasedExt.clone()
process.eIDLoose.electronIDType = 'classbased'
process.eIDLoose.electronQuality = 'loose'

process.eIDTight = eidCutBasedExt.clone()
process.eIDTight.electronIDType = 'classbased'
process.eIDTight.electronQuality = 'tight'

process.eIDSequence = cms.Sequence(process.eIDRobustLoose+
                                   process.eIDRobustLooseV00+
                                   process.eIDRobustTight+
                                   process.eIDRobustHighEnergy+
                                   process.eIDLoose+
                                   process.eIDTight)

process.load("RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff") 

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

z1lepton1pt = 20.0 
z1lepton2pt = 10.0 
process.higgsAnalysis = cms.EDFilter( "Higgs", 
        DoLog            = cms.bool( False ), 
        muonTag          = cms.InputTag( "muons" ), 
        minZ1mu1pt       = cms.double( z1lepton1pt ), 
        minZ1mu2pt       = cms.double( z1lepton2pt ), 
        minZ2mu2pt       = cms.double( 5.0 ), 
        maxMuAbsEta      = cms.double( 2.4 ), 
        electronTag      = cms.InputTag( "gsfElectrons" ), 
        electronMap      = cms.InputTag( "eIDRobustLoose" ), 
        electronCutValue = cms.int32( 7 ), 
        minZ1e1pt        = cms.double( z1lepton1pt ), 
        minZ1e2pt        = cms.double( z1lepton2pt ), 
        minZ2e2pt        = cms.double( 7.0 ), 
        maxElecAbsEta    = cms.double( 2.5 ), 
        photonTag        = cms.InputTag( "photons" ), 
        fjPhotonRho      = cms.InputTag( "kt6PFJets25","rho" ), 
        minNoTrkpt       = cms.double( z1lepton2pt ), 
        minNoTrkAbsEta   = cms.double( 2.5 ), 
        maxNoTrkAbsEta   = cms.double( 3.0 ), 
        hfTag            = cms.InputTag( "hfRecoEcalCandidate" ), 
        minHFpt          = cms.double( z1lepton2pt ), 
        minHFElecAbsEta  = cms.double( 3.0 ), 
        maxHFElecAbsEta  = cms.double( 5.0 ), 
        minZ1Mass        = cms.double( 50.0 ),
        maxZ1Mass        = cms.double( 120.0 ),
        minZ2Mass        = cms.double( 12.0 ),
        maxZ2Mass        = cms.double( 120.0 ),
        min4objMass      = cms.double( 100.0 ),
        electronRelIsoLimit = cms.double( 100.0 ),
        muonRelIsoLimit     = cms.double( 100.0 ), 

        pileupEra = cms.int32( 20111 ),
        systPileupShift = cms.double( 0.4 )
)


# process.AnalysisIntroSequence = cms.Sequence( process.hfEMClusteringSequence ) 
process.AnalysisIntroSequence = cms.Sequence( process.eIDSequence + process.fjSequence25 + process.hfRecoEcalCandidate ) 

process.higgsGen = cms.EDFilter ( "HiggsGenAnalysis",
    filterFarElectronsOnly = cms.bool( False ) 
)

process.p = cms.Path(
    process.AnalysisIntroSequence + 
    process.higgsGen + 
    process.higgsAnalysis
)

 
#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("GluGlu_120GeV_Higgs_MC_Reco_Analysis_Cutlevel_0.root"),
)

