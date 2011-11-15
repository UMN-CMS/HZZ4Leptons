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
process.pdf0.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf0.pdfReweightTargetId = 0

process.pdf1 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf1.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf1.pdfReweightTargetId = 1

process.pdf2 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf2.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf2.pdfReweightTargetId = 2

process.pdf3 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf3.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf3.pdfReweightTargetId = 3

process.pdf4 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf4.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf4.pdfReweightTargetId = 4

process.pdf5 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf5.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf5.pdfReweightTargetId = 5

process.pdf6 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf6.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf6.pdfReweightTargetId = 6

process.pdf7 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf7.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf7.pdfReweightTargetId = 7

process.pdf8 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf8.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf8.pdfReweightTargetId = 8

process.pdf9 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf9.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf9.pdfReweightTargetId = 9

process.pdf10 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf10.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf10.pdfReweightTargetId = 10

process.pdf11 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf11.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf11.pdfReweightTargetId = 11

process.pdf12 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf12.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf12.pdfReweightTargetId = 12

process.pdf13 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf13.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf13.pdfReweightTargetId = 13

process.pdf14 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf14.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf14.pdfReweightTargetId = 14

process.pdf15 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf15.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf15.pdfReweightTargetId = 15

process.pdf16 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf16.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf16.pdfReweightTargetId = 16

process.pdf17 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf17.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf17.pdfReweightTargetId = 17

process.pdf18 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf18.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf18.pdfReweightTargetId = 18

process.pdf19 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf19.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf19.pdfReweightTargetId = 19

process.pdf20 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf20.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf20.pdfReweightTargetId = 20

process.pdf21 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf21.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf21.pdfReweightTargetId = 21

process.pdf22 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf22.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf22.pdfReweightTargetId = 22

process.pdf23 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf23.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf23.pdfReweightTargetId = 23

process.pdf24 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf24.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf24.pdfReweightTargetId = 24

process.pdf25 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf25.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf25.pdfReweightTargetId = 25

process.pdf26 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf26.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf26.pdfReweightTargetId = 26

process.pdf27 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf27.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf27.pdfReweightTargetId = 27

process.pdf28 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf28.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf28.pdfReweightTargetId = 28

process.pdf29 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf29.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf29.pdfReweightTargetId = 29

process.pdf30 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf30.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf30.pdfReweightTargetId = 30

process.pdf31 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf31.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf31.pdfReweightTargetId = 31

process.pdf32 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf32.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf32.pdfReweightTargetId = 32

process.pdf33 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf33.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf33.pdfReweightTargetId = 33

process.pdf34 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf34.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf34.pdfReweightTargetId = 34

process.pdf35 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf35.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf35.pdfReweightTargetId = 35

process.pdf36 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf36.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf36.pdfReweightTargetId = 36

process.pdf37 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf37.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf37.pdfReweightTargetId = 37

process.pdf38 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf38.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf38.pdfReweightTargetId = 38

process.pdf39 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf39.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf39.pdfReweightTargetId = 39

process.pdf40 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf40.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf40.pdfReweightTargetId = 40

process.pdf41 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf41.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf41.pdfReweightTargetId = 41

process.pdf42 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf42.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf42.pdfReweightTargetId = 42

process.pdf43 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf43.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf43.pdfReweightTargetId = 43

process.pdf44 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf44.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf44.pdfReweightTargetId = 44

process.pdf45 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf45.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf45.pdfReweightTargetId = 45

process.pdf46 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf46.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf46.pdfReweightTargetId = 46

process.pdf47 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf47.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf47.pdfReweightTargetId = 47

process.pdf48 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf48.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf48.pdfReweightTargetId = 48

process.pdf49 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf49.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf49.pdfReweightTargetId = 49

process.pdf50 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf50.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf50.pdfReweightTargetId = 50

process.pdf51 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf51.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf51.pdfReweightTargetId = 51

process.pdf52 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf52.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf52.pdfReweightTargetId = 52

process.pdf53 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf53.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf53.pdfReweightTargetId = 53

process.pdf54 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf54.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf54.pdfReweightTargetId = 54

process.pdf55 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf55.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf55.pdfReweightTargetId = 55

process.pdf56 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf56.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf56.pdfReweightTargetId = 56

process.pdf57 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf57.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf57.pdfReweightTargetId = 57

process.pdf58 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf58.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf58.pdfReweightTargetId = 58

process.pdf59 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf59.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf59.pdfReweightTargetId = 59

process.pdf60 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf60.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf60.pdfReweightTargetId = 60

process.pdf61 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf61.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf61.pdfReweightTargetId = 61

process.pdf62 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf62.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf62.pdfReweightTargetId = 62

process.pdf63 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf63.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf63.pdfReweightTargetId = 63

process.pdf64 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf64.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf64.pdfReweightTargetId = 64

process.pdf65 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf65.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf65.pdfReweightTargetId = 65

process.pdf66 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf66.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf66.pdfReweightTargetId = 66

process.pdf67 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf67.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf67.pdfReweightTargetId = 67

process.pdf68 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf68.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf68.pdfReweightTargetId = 68

process.pdf69 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf69.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf69.pdfReweightTargetId = 69

process.pdf70 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf70.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf70.pdfReweightTargetId = 70

process.pdf71 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf71.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf71.pdfReweightTargetId = 71

process.pdf72 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf72.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf72.pdfReweightTargetId = 72

process.pdf73 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf73.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf73.pdfReweightTargetId = 73

process.pdf74 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf74.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf74.pdfReweightTargetId = 74

process.pdf75 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf75.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf75.pdfReweightTargetId = 75

process.pdf76 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf76.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf76.pdfReweightTargetId = 76

process.pdf77 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf77.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf77.pdfReweightTargetId = 77

process.pdf78 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf78.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf78.pdfReweightTargetId = 78

process.pdf79 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf79.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf79.pdfReweightTargetId = 79

process.pdf80 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf80.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf80.pdfReweightTargetId = 80

process.pdf81 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf81.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf81.pdfReweightTargetId = 81

process.pdf82 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf82.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf82.pdfReweightTargetId = 82

process.pdf83 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf83.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf83.pdfReweightTargetId = 83

process.pdf84 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf84.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf84.pdfReweightTargetId = 84

process.pdf85 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf85.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf85.pdfReweightTargetId = 85

process.pdf86 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf86.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf86.pdfReweightTargetId = 86

process.pdf87 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf87.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf87.pdfReweightTargetId = 87

process.pdf88 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf88.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf88.pdfReweightTargetId = 88

process.pdf89 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf89.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf89.pdfReweightTargetId = 89

process.pdf90 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf90.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf90.pdfReweightTargetId = 90

process.pdf91 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf91.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf91.pdfReweightTargetId = 91

process.pdf92 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf92.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf92.pdfReweightTargetId = 92

process.pdf93 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf93.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf93.pdfReweightTargetId = 93

process.pdf94 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf94.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf94.pdfReweightTargetId = 94

process.pdf95 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf95.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf95.pdfReweightTargetId = 95

process.pdf96 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf96.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf96.pdfReweightTargetId = 96

process.pdf97 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf97.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf97.pdfReweightTargetId = 97

process.pdf98 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf98.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf98.pdfReweightTargetId = 98

process.pdf99 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf99.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf99.pdfReweightTargetId = 99

process.pdf100 =ZShape.EffAcc.ZEfficiencyPDF_cfi.mcEff.clone()
process.pdf100.pdfReweightTargetName = 'NNPDF21_100.LHgrid'
process.pdf100.pdfReweightTargetId = 100

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
     + process.pdf53
     + process.pdf54
     + process.pdf55
     + process.pdf56
     + process.pdf57
     + process.pdf58
     + process.pdf59
     + process.pdf60
     + process.pdf61
     + process.pdf62
     + process.pdf63
     + process.pdf64
     + process.pdf65
     + process.pdf66
     + process.pdf67
     + process.pdf68
     + process.pdf69
     + process.pdf70
     + process.pdf71
     + process.pdf72
     + process.pdf73
     + process.pdf74
     + process.pdf75
     + process.pdf76
     + process.pdf77
     + process.pdf78
     + process.pdf79
     + process.pdf80
     + process.pdf81
     + process.pdf82
     + process.pdf83
     + process.pdf84
     + process.pdf85
     + process.pdf86
     + process.pdf87
     + process.pdf88
     + process.pdf89
     + process.pdf90
     + process.pdf91
     + process.pdf92
     + process.pdf93
     + process.pdf94
     + process.pdf95
     + process.pdf96
     + process.pdf97
     + process.pdf98
     + process.pdf99
     + process.pdf100
  )

