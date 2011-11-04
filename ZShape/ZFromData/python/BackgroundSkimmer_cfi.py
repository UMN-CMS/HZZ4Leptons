import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *

#from PhysicsTools.TagAndProbe.tag_probe_electron_cfi import *

from RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff import *

from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import *

from Geometry.CaloEventSetup.CaloGeometry_cff import *

from Configuration.EventContent.EventContent_cff import *

hfRecoEcalCandidate.Correct = True
##hfRecoEcalCandidate.e9e25Cut = 0
hfRecoEcalCandidate.intercept2DCut = 0.25

hfSuperClusterCandidate = hfRecoEcalCandidate.clone()
hfSuperClusterCandidate.e9e25Cut = 0
hfSuperClusterCandidate.intercept2DCut = -99


    ############################################# 
    ############  GsfElectrons  ################
    #############################################
 

electrons = cms.EDFilter("ElectronDuplicateRemover",
    src = cms.untracked.string('pixelMatchGsfElectrons'),
    ptMin = cms.untracked.double(15.0),
    EndcapMinEta = cms.untracked.double(1.56),
    ptMax = cms.untracked.double(1000.0),
    BarrelMaxEta = cms.untracked.double(1.4442),
    EndcapMaxEta = cms.untracked.double(3.0)
)

theGsfElectrons = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("electrons"),
    cut = cms.string('pt > 15.0')
)

HFElectronID = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfRecoEcalCandidate"),
    cut = cms.string('et > 10.0')
)

theGsfHf = cms.EDFilter("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("HFElectronID"))
)


tpMapEnd = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(200.0),
    TagCollection = cms.InputTag("theGsfElectrons"),
    MassMinCut = cms.untracked.double(50.0),
    ProbeCollection = cms.InputTag("theGsfHf")
)

#electron_sequence = cms.Sequence(electrons * theGsfElectrons * theGsfHf * theIsolation * eidRobust * theId * theHLT * HFElectronID )

electron_sequence = cms.Sequence(electrons * theGsfElectrons *  HFElectronID * theGsfHf * tpMapEnd )




#lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence * truthMatch_sequence)
#lepton_cands = cms.Sequence(genParticles *  ZIntoElectronsEventProducer * SmearedElectronsProducer * hfEMClusteringSequence * sc_sequence * electrons * theGsfElectrons)
lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * electron_sequence )


   

