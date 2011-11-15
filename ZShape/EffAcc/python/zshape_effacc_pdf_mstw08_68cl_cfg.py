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
process.pdf0.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf0.pdfReweightTargetId = 0

process.pdf1 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf1.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf1.pdfReweightTargetId = 1

process.pdf2 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf2.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf2.pdfReweightTargetId = 2

process.pdf3 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf3.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf3.pdfReweightTargetId = 3

process.pdf4 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf4.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf4.pdfReweightTargetId = 4

process.pdf5 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf5.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf5.pdfReweightTargetId = 5

process.pdf6 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf6.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf6.pdfReweightTargetId = 6

process.pdf7 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf7.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf7.pdfReweightTargetId = 7

process.pdf8 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf8.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf8.pdfReweightTargetId = 8

process.pdf9 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf9.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf9.pdfReweightTargetId = 9

process.pdf10 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf10.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf10.pdfReweightTargetId = 10

process.pdf11 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf11.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf11.pdfReweightTargetId = 11

process.pdf12 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf12.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf12.pdfReweightTargetId = 12

process.pdf13 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf13.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf13.pdfReweightTargetId = 13

process.pdf14 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf14.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf14.pdfReweightTargetId = 14

process.pdf15 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf15.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf15.pdfReweightTargetId = 15

process.pdf16 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf16.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf16.pdfReweightTargetId = 16

process.pdf17 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf17.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf17.pdfReweightTargetId = 17

process.pdf18 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf18.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf18.pdfReweightTargetId = 18

process.pdf19 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf19.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf19.pdfReweightTargetId = 19

process.pdf20 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf20.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf20.pdfReweightTargetId = 20

process.pdf21 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf21.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf21.pdfReweightTargetId = 21

process.pdf22 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf22.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf22.pdfReweightTargetId = 22

process.pdf23 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf23.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf23.pdfReweightTargetId = 23

process.pdf24 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf24.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf24.pdfReweightTargetId = 24

process.pdf25 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf25.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf25.pdfReweightTargetId = 25

process.pdf26 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf26.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf26.pdfReweightTargetId = 26

process.pdf27 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf27.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf27.pdfReweightTargetId = 27

process.pdf28 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf28.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf28.pdfReweightTargetId = 28

process.pdf29 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf29.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf29.pdfReweightTargetId = 29

process.pdf30 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf30.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf30.pdfReweightTargetId = 30

process.pdf31 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf31.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf31.pdfReweightTargetId = 31

process.pdf32 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf32.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf32.pdfReweightTargetId = 32

process.pdf33 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf33.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf33.pdfReweightTargetId = 33

process.pdf34 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf34.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf34.pdfReweightTargetId = 34

process.pdf35 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf35.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf35.pdfReweightTargetId = 35

process.pdf36 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf36.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf36.pdfReweightTargetId = 36

process.pdf37 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf37.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf37.pdfReweightTargetId = 37

process.pdf38 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf38.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf38.pdfReweightTargetId = 38

process.pdf39 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf39.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf39.pdfReweightTargetId = 39

process.pdf40 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf40.pdfReweightTargetName = 'MSTW2008nlo68cl.LHgrid'
process.pdf40.pdfReweightTargetId = 40

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
     + process.pdf15
     + process.pdf16
     + process.pdf17
     + process.pdf18
     + process.pdf19
     + process.pdf20
     + process.pdf21
     + process.pdf22
     + process.pdf23
     + process.pdf24
     + process.pdf25
     + process.pdf26
     + process.pdf27
     + process.pdf28
     + process.pdf29
     + process.pdf30
     + process.pdf31
     + process.pdf32
     + process.pdf33
     + process.pdf34
     + process.pdf35
     + process.pdf36
     + process.pdf37
     + process.pdf38
     + process.pdf39
     + process.pdf40
  )

