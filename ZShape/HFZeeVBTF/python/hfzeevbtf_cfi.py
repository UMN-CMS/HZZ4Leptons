import FWCore.ParameterSet.Config as cms

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi


hfzeevbtf = cms.EDFilter('HFZeeVBTF',
                    acceptedElectronIDs = cms.vint32( 7 ),
                    ECALid = cms.string('simpleEleId80cIso'),
                    minEtECAL = cms.double(20),
                    minEtHF   = cms.double(20),
                    Zmass     = cms.vdouble(40,130),
                    FILTERid  = cms.untracked.string(''),
                    # the following params are, respectively: e9e25_loose, e9e25_tight,  var2d_loose, var2d_tight,  eCOREe9_loose, eCOREe9_tight,  eSeL_loose, eSeL_tight;
                    # set to -9999/+9999 if you want to neglect a cut, resp eCOREe9/eSeL
                    hFselParams =  cms.vdouble(0.90, 0.94,      0.32, 0.45,    -9999, -9999,     9999, 9999),
                    #hFselParams =  cms.vdouble(0.90, 0.94,      0.32, 0.45,    1, 2,     3, 4),
                    #hFElectrons = cms.InputTag("hfRecoEcalCandidateLoose","",""),
                    hFElectrons = cms.InputTag("hfRecoEcalCandidate","",""),
                    myName = cms.string('HFZeeVBTF-IdIsoRej'),
                    DoLog = cms.bool(False),
                    skipTrigger = cms.bool(False),
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
