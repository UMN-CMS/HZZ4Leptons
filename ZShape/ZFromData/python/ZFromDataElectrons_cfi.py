import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *

#from PhysicsTools.TagAndProbe.tag_probe_electron_cfi import *

from RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff import *

from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import *

from Geometry.CaloEventSetup.CaloGeometry_cff import *

from Configuration.EventContent.EventContent_cff import *

hfRecoEcalCandidate.Correct = True
##hfRecoEcalCandidate.e9e25Cut = 0
##hfRecoEcalCandidate.intercept2DCut = -99

hfSuperClusterCandidate = hfRecoEcalCandidate.clone()
hfSuperClusterCandidate.e9e25Cut = 0
hfSuperClusterCandidate.intercept2DCut = -99

from RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff import *

#  Calculate efficiency for *SuperClusters passing as GsfElectron* 
#
#  Tag           =  isolated GsfElectron with Robust ID, passing HLT, and  
#                    within the fiducial volume of ECAL
#  Probe --> Passing Probe   =  
#  SC --> GsfElectron --> isolation --> id --> Trigger
#############################################################




    ############################################# 
    ############  SuperClusters  ################
    ############################################# 


HybridSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedHybridSuperClusters"),
    particleType = cms.string('gamma')
)
EBSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("HybridSuperClusters"),
    cut = cms.string('abs( eta ) < 1.4442')
)



EndcapSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    particleType = cms.string('gamma')
)
EESuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("EndcapSuperClusters"),
    cut = cms.string('abs( eta ) > 1.560 & abs( eta ) < 3.0')
)



#HFSuperClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
#    src = cms.InputTag("hfEMClusters"),
#    particleType = cms.string('gamma')
#)
#HFSuperClusters = cms.EDFilter("CandViewSelector",
#    src = cms.InputTag("HFSuperClusterCands"),
#    cut = cms.string('abs( eta ) > 3.0')
#)
#theHFSuperClusters = cms.EDFilter("CandViewSelector",
#    src = cms.InputTag("HFSuperClusters"),
#    cut = cms.string('et > 10.0')
#)

theHFSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfSuperClusterCandidate"),
    cut = cms.string('et > 10.0')
)


allSuperClusters = cms.EDFilter("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"),cms.InputTag("theHFSuperClusters"))
)

thepreSuperClusters = cms.EDFilter("PdgIdAndStatusCandViewSelector",
    src = cms.InputTag("genParticles"),
    pdgId = cms.vint32(11),
    status = cms.vint32(3)
)

# My final selection of superCluster candidates 
theSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("allSuperClusters"),
    cut = cms.string('et  > 5.0')
)

theEEHFGapSuperClusters =  cms.EDFilter("CandViewSelector",
    src = cms.InputTag("EESuperClusters"),
    cut = cms.string('(abs( eta ) > 2.5 & abs( eta ) < 3.0)')
)

sc_sequence = cms.Sequence( (HybridSuperClusters * EBSuperClusters + EndcapSuperClusters * EESuperClusters + hfSuperClusterCandidate * theHFSuperClusters) *allSuperClusters * thepreSuperClusters * (theEEHFGapSuperClusters + theSuperClusters))




    ############################################# 
    ############  GsfElectrons  ################
    #############################################
 

electrons = cms.EDFilter("ElectronDuplicateRemover",
    src = cms.untracked.string('pixelMatchGsfElectrons'),
    ptMin = cms.untracked.double(20.0),
    EndcapMinEta = cms.untracked.double(1.56),
    ptMax = cms.untracked.double(1000.0),
    BarrelMaxEta = cms.untracked.double(1.4442),
    EndcapMaxEta = cms.untracked.double(3.0)
)

theGsfElectrons = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("electrons"),
    cut = cms.string('et > 10.0 & pt < 105.0')
)

thePt20GsfElectrons = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("electrons"),
    cut = cms.string('pt > 20.0 & (( caloEnergy * sin( caloPosition.theta ) )  > 10.0)')
)

HFElectronID = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hfRecoEcalCandidate"),
    cut = cms.string('et > 10.0')
)

theGsfHf = cms.EDFilter("CandViewMerger",
    #src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"), cms.InputTag("theEEHFGapSuperClusters"))
    src = cms.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theHFSuperClusters"))
)

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
eleIsoDepositTk.src = cms.InputTag('theGsfElectrons')
eleIsoDepositEcalFromHits.src = cms.InputTag('theTrackerIsolation')
eleIsoDepositHcalFromTowers.src = cms.InputTag('theEcalIsolation')


theTrackerIsolation = cms.EDProducer("IsolatedElectronCandProducer",
	electronProducer = cms.InputTag('theGsfElectrons'),
        #electronProducer = cms.InputTag('pixelMatchGsfElectrons'),
	IsolationProducer = cms.InputTag('eleIsoFromDepsTk'),
	isoCutEB = cms.double(7.2),
	isoCutEE = cms.double(5.1),
	EBEEBound = cms.double(1.5)
)

theEcalIsolation = cms.EDProducer("IsolatedElectronCandProducer",
	electronProducer = cms.InputTag('theTrackerIsolation'),
        #electronProducer = cms.InputTag('electrons'),
        #electronProducer = cms.InputTag('pixelMatchGsfElectrons'),
	IsolationProducer = cms.InputTag('eleIsoFromDepsEcalFromHits'),
	isoCutEB = cms.double(5.7),
	isoCutEE = cms.double(5.0),
	EBEEBound = cms.double(1.5)
)

theIsolation = cms.EDProducer("IsolatedElectronCandProducer",
	electronProducer = cms.InputTag('theEcalIsolation'),
        #electronProducer = cms.InputTag('pixelMatchGsfElectrons'),
	IsolationProducer = cms.InputTag('eleIsoFromDepsHcalFromTowers'),
	isoCutEB = cms.double(8.1),
	isoCutEE = cms.double(3.4),
	EBEEBound = cms.double(1.5)
)


#  isolation  ################
#theIsolation = cms.EDProducer("IsolatedElectronCandProducer",
#    absolut = cms.bool(False),
#    trackProducer = cms.InputTag('generalTracks'),
#    isoCut = cms.double(0.2),
#    intRadius = cms.double(0.02),
#    electronProducer = cms.InputTag('theGsfElectrons'),
#    extRadius = cms.double(0.2),
#    ptMin = cms.double(1.5),
#    maxVtxDist = cms.double(0.1)
#)

# Cut-based Robust electron ID  ######
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobust = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobust.src = cms.InputTag('theIsolation')
eidRobust.robustEleIDCuts = cms.PSet(
            barrel = cms.vdouble(999.9, 0.010, 999.9, 0.0071),
            endcap = cms.vdouble(999.9, 0.028, 999.9, 0.0066)
                )


theId = cms.EDProducer("eidCandProducer",
    electronCollection = cms.untracked.InputTag('theIsolation'),  
    electronLabelLoose = cms.InputTag('eidRobust')
)

# Trigger  ##################
theHLT = cms.EDProducer("eTriggerCandProducer",
    InputProducer = cms.InputTag('theId'),              
    #hltTag = cms.untracked.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt10TrackIsolFilter","","HLT")
    hltTag = cms.untracked.InputTag("hltL1NonIsoHLTLooseIsoSingleElectronLWEt15TrackIsolFilter","","HLT") 
)

#electron_sequence = cms.Sequence(electrons * theGsfElectrons * theGsfHf * theIsolation * eidRobust * theId * theHLT * HFElectronID )
electron_sequence = cms.Sequence(electrons * theGsfElectrons * theGsfHf * thePt20GsfElectrons * eleIsoDepositTk * eleIsoFromDepsTk * theTrackerIsolation * eleIsoDepositEcalFromHits * eleIsoFromDepsEcalFromHits * theEcalIsolation * eleIsoDepositHcalFromTowers * eleIsoFromDepsHcalFromTowers * theIsolation  * eidRobust * theId * theHLT * HFElectronID )




 



    ################################################################ 
    ############  All Tag / Probe Association Maps  ###############
    ################################################################ 
    #// Remember that tag will always be "theHLT" collection.
    #//
    #// Probe can be one of the following collections: 
    #// "theSuperClusters", "theGsfElectrons", "theIsolation", "theId".
    #// 
    #// Passing Probe can be one of the following collections: 
    #// "theGsfElectrons", "theIsolation", "theId", "theHLT".
    #//




tpMapSuperClusters = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    #TagCollection = cms.InputTag("theHLT"),
    TagCollection = cms.InputTag("theSuperClusters"),
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theSuperClusters"),
    RequireOS = cms.untracked.bool(False)
)

tpMappreSuperClusters = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(130.0),
    TagCollection = cms.InputTag("theHLT"),
    MassMinCut = cms.untracked.double(50.0),
    ProbeCollection = cms.InputTag("thepreSuperClusters"),
    RequireOS = cms.untracked.bool(False)
)

tpMapGsfElectrons = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    #TagCollection = cms.InputTag("theHLT"),
    TagCollection = cms.InputTag("theGsfElectrons"),
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theGsfElectrons"),
    RequireOS = cms.untracked.bool(False)
)

tpMapIsolation = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    TagCollection = cms.InputTag("theHLT"),
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theIsolation"),
    RequireOS = cms.untracked.bool(False)
)

tpMapId = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    TagCollection = cms.InputTag("theHLT"),
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theId"),
    RequireOS = cms.untracked.bool(False)
)

tpMapHFSuperClusters = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    TagCollection = cms.InputTag("theHLT"),          #theID??
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theHFSuperClusters"),
    RequireOS = cms.untracked.bool(False)
)

tpMapGsfAndHF = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    #TagCollection = cms.InputTag("theGsfElectrons"),          #theID??
    TagCollection = cms.InputTag("theHLT"),
    MassMinCut = cms.untracked.double(60.0),
    ProbeCollection = cms.InputTag("theGsfHf"),
    RequireOS = cms.untracked.bool(False)
)

tpMap_sequence = cms.Sequence(tpMappreSuperClusters + tpMapSuperClusters + tpMapGsfElectrons + tpMapIsolation + tpMapId + tpMapHFSuperClusters + tpMapGsfAndHF)




    ################################################################ 
    ############  All Truth-matched collections  ###################
    ################################################################ 
    #// find generator particles matching by DeltaR

preSuperClustersMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("thepreSuperClusters"),
    distMin = cms.double(0.1),
    matched = cms.InputTag("genParticles")
)

SuperClustersMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theSuperClusters"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

GsfElectronsMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theGsfElectrons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

IsolationMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theIsolation"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

IdMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theId"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HLTMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theHLT"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HFSCMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("theHFSuperClusters"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

HFIDMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(11),
    src = cms.InputTag("HFElectronID"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)


truthMatch_sequence = cms.Sequence(preSuperClustersMatch + SuperClustersMatch + GsfElectronsMatch + IsolationMatch + IdMatch + HLTMatch + HFSCMatch + HFIDMatch)

lepton_cands = cms.Sequence(genParticles * hfEMClusteringSequence * sc_sequence * electron_sequence * tpMap_sequence * truthMatch_sequence)




   

