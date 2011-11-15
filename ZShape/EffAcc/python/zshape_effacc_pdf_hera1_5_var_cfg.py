import FWCore.ParameterSet.Config as cms

process = cms.Process("Zefficiency")
process.TimerService = cms.Service("TimerService")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound') 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring("file:///local/cms/user/pastika/powhegZ/base_PDF-10901/PythiaUEZ2/Zee_Powheg_PythiaUEZ2_1.root" )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_10M_partBUILDINGTTEST.root')
)

process.load("ZShape.EffAcc.SmearedElectronsProducer_cfi")
import ZShape.EffAcc.ZEfficiencyPDF_cfi 


process.pdf0 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf0.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf0.pdfReweightTargetId = 0

process.pdf1 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf1.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf1.pdfReweightTargetId = 1

process.pdf2 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf2.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf2.pdfReweightTargetId = 2

process.pdf3 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf3.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf3.pdfReweightTargetId = 3

process.pdf4 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf4.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf4.pdfReweightTargetId = 4

process.pdf5 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf5.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf5.pdfReweightTargetId = 5

process.pdf6 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf6.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf6.pdfReweightTargetId = 6

process.pdf7 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf7.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf7.pdfReweightTargetId = 7

process.pdf8 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf8.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf8.pdfReweightTargetId = 8

process.pdf9 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf9.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf9.pdfReweightTargetId = 9

process.pdf10 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf10.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf10.pdfReweightTargetId = 10

process.pdf11 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf11.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf11.pdfReweightTargetId = 11

process.pdf12 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf12.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf12.pdfReweightTargetId = 12

process.pdf13 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf13.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf13.pdfReweightTargetId = 13

process.pdf14 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf14.pdfReweightTargetName = 'HERAPDF1.5_VAR.LHgrid'
process.pdf14.pdfReweightTargetId = 14

process.p=cms.Path(process.SmearedElectronsProducer
     + process.pdf0
     + process.pdf1
     + process.pdf2
     + process.pdf3
     + process.pdf4
     + process.pdf5
     + process.pdf6
     + process.pdf7
     + process.pdf8
     + process.pdf9
     + process.pdf10
     + process.pdf11
     + process.pdf12
     + process.pdf13
     + process.pdf14
  )

