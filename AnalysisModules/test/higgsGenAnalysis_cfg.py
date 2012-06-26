import FWCore.ParameterSet.Config as cms

import os

needGenParts=True

process = cms.Process("AGEN");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/Fall11/GluGluToHToZZTo4L_M-115_7TeV-powheg-pythia6/AODSIM/PU_S6_START44_V9B-v1/0000/0A49A5CA-712E-E111-BCE7-00215E93E7AC.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("genanalysis8.root"),
)

#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("HeavyNu.AnalysisModules.heavynugenlevel_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.higgsGen = cms.EDFilter ( "HiggsGenAnalysis"
)

process.p = cms.Path( process.higgsGen )

process.higgsGen = cms.EDFilter ( "HiggsGenAnalysis",
    filterFarElectronsOnly = cms.bool( False ) 
)

