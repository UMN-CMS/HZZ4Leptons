import FWCore.ParameterSet.Config as cms

process = cms.Process("SCwithPUAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TimerService = cms.Service("TimerService")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(3000)
        )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('SCwithPUAnalysis_Hgg_PU20_v2.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
#        'file:/data2/kubota/CMSSW_3_9_4/src/H120GGgluonfusion_PU20_1.root',
        'file:/data/franzoni/data/superClusterPU/H120GGgluonfusion_PU20_9.root'
#        'rfio:/castor/cern.ch/user/f/futyand/SCwithPU/H120GGgluonfusion_PU20_1.root'
#        'rfio:/castor/cern.ch/user/f/futyand/SCwithPU/H120GGgluonfusion_PU20_2.root',
    )
)

process.scwithpuanalyzer = cms.EDAnalyzer('SCwithPUAnalysis'
)


process.p = cms.Path(process.scwithpuanalyzer)
