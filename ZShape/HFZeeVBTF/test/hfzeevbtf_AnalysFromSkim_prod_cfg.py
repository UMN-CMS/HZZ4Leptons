import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
                            fileMode = cms.untracked.string('NOMERGE')
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(    )
)
process.load("RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('')
)

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

process.Id-Iso-Rej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter
process.Id-Iso = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which I carry along to allow the keeping also of electrons that have not passed conversion rejection
    acceptedElectronIDs = cms.vint32( 3, 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())


)

process.Id-Rej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 5 , 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)

process.Iso-Rej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 6 , 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)






process.Id = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 1 , 3, 5, 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


process.Iso = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


process.Rej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
#      ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 4, 5, 6 , 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)






process.hfRecoEcalCandidate.intercept2DCut=0.3

# 
# import HLTrigger.HLTfilters.hltHighLevelDev_cfi
# process.EG_1e28 = HLTrigger.HLTfilters.hltHighLevelDev_cfi.hltHighLevelDev.clone(andOr = True)
# process.EG_1e28.HLTPaths = (
# #     "HLT_Mu9",
# #     "HLT_Mu5_Track0_Jpsi"
#     "HLT_Photon10_L1R",
#     "HLT_Photon15_L1R",
#     "HLT_Photon15_LooseEcalIso_L1R",
#     "HLT_Photon20_L1R",
#     "HLT_Photon30_L1R_8E29",
#     "HLT_DoublePhoton4_Jpsi_L1R",
#     "HLT_DoublePhoton4_Upsilon_L1R",
#     "HLT_DoublePhoton4_eeRes_L1R",
#     "HLT_DoublePhoton5_eeRes_L1R", #added to match the /cdaq/physics/firstCollisions10/v2.0/HLT_7TeV/V5 table
#     "HLT_DoublePhoton5_Jpsi_L1R",
#     "HLT_DoublePhoton5_Upsilon_L1R",
#     "HLT_DoublePhoton5_L1R",
#     "HLT_DoublePhoton10_L1R",
#     "HLT_DoubleEle5_SW_L1R",
#     "HLT_Ele20_LW_L1R",
#     "HLT_Ele15_SiStrip_L1R",
#     "HLT_Ele15_SC10_LW_L1R",
#     "HLT_Ele15_LW_L1R",
#     "HLT_Ele10_LW_EleId_L1R",
#     "HLT_Ele10_LW_L1R",
#     "HLT_Photon15_TrackIso_L1R"
#     )
# process.EG_1e28.HLTPathsPrescales  = cms.vuint32(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# process.EG_1e28.HLTOverallPrescale = cms.uint32(1)
# process.EG_1e28.throw = False
# process.EG_1e28.andOr = True


process.p = cms.Path(
    # process.EG_1e28 *
    process.hfRecoEcalCandidate
    *process.demo
    *process.demoLoose
    )
