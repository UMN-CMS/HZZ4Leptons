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
process.pdf0.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf0.pdfReweightTargetId = 0

process.pdf1 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf1.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf1.pdfReweightTargetId = 1

process.pdf2 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf2.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf2.pdfReweightTargetId = 2

process.pdf3 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf3.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf3.pdfReweightTargetId = 3

process.pdf4 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf4.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf4.pdfReweightTargetId = 4

process.pdf5 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf5.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf5.pdfReweightTargetId = 5

process.pdf6 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf6.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf6.pdfReweightTargetId = 6

process.pdf7 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf7.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf7.pdfReweightTargetId = 7

process.pdf8 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf8.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf8.pdfReweightTargetId = 8

process.pdf9 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf9.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf9.pdfReweightTargetId = 9

process.pdf10 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf10.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf10.pdfReweightTargetId = 10

process.pdf11 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf11.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf11.pdfReweightTargetId = 11

process.pdf12 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf12.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf12.pdfReweightTargetId = 12

process.pdf13 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf13.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf13.pdfReweightTargetId = 13

process.pdf14 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf14.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf14.pdfReweightTargetId = 14

process.pdf15 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf15.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf15.pdfReweightTargetId = 15

process.pdf16 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf16.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf16.pdfReweightTargetId = 16

process.pdf17 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf17.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf17.pdfReweightTargetId = 17

process.pdf18 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf18.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf18.pdfReweightTargetId = 18

process.pdf19 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf19.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf19.pdfReweightTargetId = 19

process.pdf20 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf20.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf20.pdfReweightTargetId = 20

process.pdf21 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf21.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf21.pdfReweightTargetId = 21

process.pdf22 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf22.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf22.pdfReweightTargetId = 22

process.pdf23 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf23.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf23.pdfReweightTargetId = 23

process.pdf24 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf24.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf24.pdfReweightTargetId = 24

process.pdf25 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf25.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf25.pdfReweightTargetId = 25

process.pdf26 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf26.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf26.pdfReweightTargetId = 26

process.pdf27 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf27.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf27.pdfReweightTargetId = 27

process.pdf28 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf28.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf28.pdfReweightTargetId = 28

process.pdf29 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf29.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf29.pdfReweightTargetId = 29

process.pdf30 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf30.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf30.pdfReweightTargetId = 30

process.pdf31 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf31.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf31.pdfReweightTargetId = 31

process.pdf32 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf32.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf32.pdfReweightTargetId = 32

process.pdf33 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf33.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf33.pdfReweightTargetId = 33

process.pdf34 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf34.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf34.pdfReweightTargetId = 34

process.pdf35 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf35.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf35.pdfReweightTargetId = 35

process.pdf36 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf36.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf36.pdfReweightTargetId = 36

process.pdf37 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf37.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf37.pdfReweightTargetId = 37

process.pdf38 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf38.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf38.pdfReweightTargetId = 38

process.pdf39 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf39.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf39.pdfReweightTargetId = 39

process.pdf40 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf40.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf40.pdfReweightTargetId = 40

process.pdf41 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf41.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf41.pdfReweightTargetId = 41

process.pdf42 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf42.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf42.pdfReweightTargetId = 42

process.pdf43 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf43.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf43.pdfReweightTargetId = 43

process.pdf44 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf44.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf44.pdfReweightTargetId = 44

process.pdf45 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf45.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf45.pdfReweightTargetId = 45

process.pdf46 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf46.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf46.pdfReweightTargetId = 46

process.pdf47 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf47.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf47.pdfReweightTargetId = 47

process.pdf48 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf48.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf48.pdfReweightTargetId = 48

process.pdf49 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf49.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf49.pdfReweightTargetId = 49

process.pdf50 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf50.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf50.pdfReweightTargetId = 50

process.pdf51 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf51.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf51.pdfReweightTargetId = 51

process.pdf52 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf52.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf52.pdfReweightTargetId = 52

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
     + process.pdf41
     + process.pdf42
     + process.pdf43
     + process.pdf44
     + process.pdf45
     + process.pdf46
     + process.pdf47
     + process.pdf48
     + process.pdf49
     + process.pdf50
     + process.pdf51
     + process.pdf52
  )

