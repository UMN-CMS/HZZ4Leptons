import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")


process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring(     )
                  
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
process.GlobalTag.globaltag = cms.string('GR_R_35X_V8B::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
##
#
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## MET creation     <=== WARNING: YOU MAY WANT TO MODIFY THIS PART OF THE CODE       %%%%%%%%%%%%%
##                                specify the names of the MET collections that you need here %%%%
##                                                                                             #%%
## if you don't specify anything the default MET is the raw Calo MET                           #%%
## if you don't specify anything the default MET is the raw Calo MET                           #%%
process.caloMET = process.patMETs.clone(                                                       #%%
        metSource = cms.InputTag("met","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
process.tcMET = process.patMETs.clone(                                                         #%%
        metSource = cms.InputTag("tcMet","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
process.pfMET = process.patMETs.clone(                                                         #%%
        metSource = cms.InputTag("pfMet","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
## specify here what you want to have on the plots! <===== MET THAT YOU WANT ON THE PLOTS  %%%%%%%
myMetCollection   = 'caloMET'
myPfMetCollection =   'pfMET'
myTcMetCollection =   'tcMET'

## specify here what you want to have on the plots! <===== MET THAT YOU WANT ON THE PLOTS  %%%%%%%
## myDesiredMetCollection = 'layer1RawCaloMETs'
## modify the sequence of the MET creation:                                                    #%%
process.makePatMETs = cms.Sequence(process.caloMET*process.tcMET*process.pfMET)

## GF changed here: more than one pat done here

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## modify the final pat sequence: keep only electrons + METS (muons are needed for met corrections)
process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#process.patElectronIsolation = cms.Sequence(process.egammaIsolationSequence)

process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),    
    )
##
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)
##
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons)
# process.makePatMuons may be needed depending on how you calculate the MET
process.makePatCandidates = cms.Sequence(process.makePatElectrons+process.makePatMETs)
process.patDefaultSequence = cms.Sequence(process.makePatCandidates)


##
##  ################################################################################
##
##  generator level filter
##
##  ################################################################################
##  conditions on different particles are OR-ed
process.genFilter = cms.EDFilter("MCSingleParticleFilter",
                                 Status = cms.untracked.vint32(1, 1),
                                 MinEta = cms.untracked.vdouble(-2.4, -2.4),
                                 MaxEta = cms.untracked.vdouble(2.4, 2.4),
                                 MinPt = cms.untracked.vdouble(15.0, 15.0),
                                 ParticleID = cms.untracked.vint32(11, -11) ## these are electron and positron
                                 )





##
##  ################################################################################
##
##  the filter to select the candidates from the data samples
##
##
## WARNING: you may want to modify this item:
HLT_process_name = "HLT"   # 
# trigger path selection
# HLT_path_name    = "HLT_Ele15_LW_L1R"
HLT_path_name     = "HLT_Photon10_L1R"
# trigger filter name
HLT_filter_name  =  "hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter"
# HLT_filter_name  =  "hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter"
#
HLT_path_name_extra   = "HLT_Photon15_L1R" #= "HLT_Ele15_LW_L1R" #
HLT_filter_name_extra = "hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter"
process.z1legFilter=cms.EDFilter('WenuCandidateFilter',
                                  ### the input collections needed:
                                  electronCollectionTag = cms.untracked.InputTag("patElectrons","","PAT"),
                                  metCollectionTag = cms.untracked.InputTag(myMetCollection,"","PAT"),
                                  pfMetCollectionTag = cms.untracked.InputTag(myPfMetCollection,"","PAT"),
                                  tcMetCollectionTag = cms.untracked.InputTag(myTcMetCollection,"","PAT"),
                                  triggerCollectionTag = cms.untracked.InputTag("TriggerResults","",HLT_process_name),
                                  triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLT_process_name),
                                  hltpath = cms.untracked.string(HLT_path_name),
                                  hltpathFilter = cms.untracked.InputTag(HLT_filter_name,"",HLT_process_name),
                                  ebRecHits = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                                  eeRecHits = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                                  PrimaryVerticesCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
                                  ### here the preselection is applied
                                  # fiducial cuts:
                                  BarrelMaxEta = cms.untracked.double(1.4442),
                                  EndCapMinEta = cms.untracked.double(1.566),
                                  EndCapMaxEta = cms.untracked.double(2.5),
                                  # demand ecal driven electron:
                                  useEcalDrivenElectrons = cms.untracked.bool(True),
                                  # demand offline spike cleaning with the Swiss Cross criterion:
                                  useSpikeRejection = cms.untracked.bool(True),
                                  spikeCleaningSwissCrossCut = cms.untracked.double(0.95),
                                  # demand geometrically matched to an HLT object with ET>15GeV
                                  useTriggerInfo = cms.untracked.bool(False), # GF 
                                  electronMatched2HLT = cms.untracked.bool(False), # GF 
                                  electronMatched2HLT_DR = cms.untracked.double(0.2), # GF 
                                  useHLTObjectETCut = cms.untracked.bool(False), # GF 
                                  hltObjectETCut = cms.untracked.double(15.),
                                  useExtraTrigger = cms.untracked.bool(False), # GF 
                                  hltpathExtra = cms.untracked.string(HLT_path_name_extra),
                                  hltpathFilterExtra = cms.untracked.InputTag(HLT_filter_name_extra,"",HLT_process_name),
                                  # ET Cut in the SC
                                  ETCut = cms.untracked.double(10.), # GF : mistakenly set to 20 in the first pass
                                  METCut = cms.untracked.double(0.),
                                  # reject events with a 2nd electron with ET > 20 that passes the WP95%
                                  vetoSecondElectronEvents = cms.untracked.bool(False),
                                  storeSecondElectron = cms.untracked.bool(False),
                                  ETCut2ndEle = cms.untracked.double(20.),
                                  vetoSecondElectronIDType = cms.untracked.string("simpleEleId80relIso"), # GF
                                  vetoSecondElectronIDSign = cms.untracked.string("="),
                                  vetoSecondElectronIDValue = cms.untracked.double(7.),
                                  # Other parameters of the code - leave them as they are
                                  useValidFirstPXBHit = cms.untracked.bool(False),
                                  useConversionRejection = cms.untracked.bool(False),
                                  useExpectedMissingHits = cms.untracked.bool(False),
                                  maxNumberOfExpectedMissingHits = cms.untracked.int32(1),
                                  # calculate some new cuts
                                  calculateValidFirstPXBHit = cms.untracked.bool(True),
                                  calculateConversionRejection = cms.untracked.bool(True),
                                  calculateExpectedMissingHits = cms.untracked.bool(True),
                                  # we are dealing with DATA
                                  dataMagneticFieldSetUp = cms.untracked.bool(False),
                                  dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
                                  )

# = cms.EDFilter('WenuCandidateFilter',
#                                   # cfg for data
#                                   dataMagneticFieldSetUp = cms.untracked.bool(False), # set to False in MC!
#                                   dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
#                                   # cuts
#                                   ETCut = cms.untracked.double(10.),
#                                   METCut = cms.untracked.double(0.),
#                                   # 2nd electron in W events
#                                   vetoSecondElectronEvents = cms.untracked.bool(False),
#                                   ETCut2ndEle = cms.untracked.double(20.),
#                                   vetoSecondElectronIDType = cms.untracked.string("simpleEleId80relIso"),
#                                   vetoSecondElectronIDSign = cms.untracked.string("="),
#                                   vetoSecondElectronIDValue = cms.untracked.double(7.),
#                                   # trigger 
#                                   useTriggerInfo = cms.untracked.bool(False),
#                                   triggerCollectionTag = cms.untracked.InputTag("TriggerResults","",HLT_process_name),
#                                   triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLT_process_name),
#                                   hltpath = cms.untracked.string(HLT_path_name), 
#                                   hltpathFilter = cms.untracked.InputTag(HLT_filter_name,"",HLT_process_name),
#                                   electronMatched2HLT = cms.untracked.bool(False),
#                                   electronMatched2HLT_DR = cms.untracked.double(0.2),
#                                   # additional preselection cuts
#                                   useValidFirstPXBHit = cms.untracked.bool(False),
#                                   useConversionRejection = cms.untracked.bool(False),
#                                   useExpectedMissingHits = cms.untracked.bool(False),
#                                   maxNumberOfExpectedMissingHits = cms.untracked.int32(1),
#                                   # calculate some new cuts
#                                   calculateValidFirstPXBHit = cms.untracked.bool(True),
#                                   calculateConversionRejection = cms.untracked.bool(True),
#                                   calculateExpectedMissingHits = cms.untracked.bool(True),
#                                   # electrons and MET
#                                   electronCollectionTag = cms.untracked.InputTag("patElectrons","","PAT"),
#                                   metCollectionTag = cms.untracked.InputTag(myDesiredMetCollection,"","PAT"),
# 
#                                   )
####################################################################################







# to access values of EldId cuts as defined by the EWK group
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter
process.demo = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 7 ),
                              
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


# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter
process.demoLoose = cms.EDAnalyzer('HFZeeVBTF',
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


# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run BEFORE the Wenu EWK filter, so to have normalization of number of MC events processed
process.demoBefCuts = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90cIso'),
    DoLog = cms.bool(True),
#   this instance of the analysis code is just to count events, so it does not matter how tight eleID might be
    acceptedElectronIDs = cms.vint32( 7 ),
                                     
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



# THIS IS PROD
process.TFileService = cms.Service("TFileService",
       fileName = cms.string(""),
)



# this is the filter requiring the same HLT bits which compose the EG dataset
# see: https://twiki.cern.ch/twiki/bin/view/CMS/PhysicsSecondaryDatasets  AND
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Configuration/Skimming/test/SDmaker_6SD_2CS_PDMinBias_1e28_cfg.py?revision=1.1&view=markup
import HLTrigger.HLTfilters.hltHighLevelDev_cfi
process.EG_1e28 = HLTrigger.HLTfilters.hltHighLevelDev_cfi.hltHighLevelDev.clone(andOr = True)
process.EG_1e28.HLTPaths = (
#     "HLT_Mu9",
#     "HLT_Mu5_Track0_Jpsi"
    "HLT_Photon10_L1R",
    "HLT_Photon15_L1R",
    "HLT_Photon15_LooseEcalIso_L1R",
    "HLT_Photon20_L1R",
    "HLT_Photon30_L1R_8E29",
    "HLT_DoublePhoton4_Jpsi_L1R",
    "HLT_DoublePhoton4_Upsilon_L1R",
    "HLT_DoublePhoton4_eeRes_L1R",
    "HLT_DoublePhoton5_eeRes_L1R", #added to match the /cdaq/physics/firstCollisions10/v2.0/HLT_7TeV/V5 table
    "HLT_DoublePhoton5_Jpsi_L1R",
    "HLT_DoublePhoton5_Upsilon_L1R",
    "HLT_DoublePhoton5_L1R",
    "HLT_DoublePhoton10_L1R",
    "HLT_DoubleEle5_SW_L1R",
    "HLT_Ele20_LW_L1R",
    "HLT_Ele15_SiStrip_L1R",
    "HLT_Ele15_SC10_LW_L1R",
    "HLT_Ele15_LW_L1R",
    "HLT_Ele10_LW_EleId_L1R",
    "HLT_Ele10_LW_L1R",
    "HLT_Photon15_TrackIso_L1R"
    )
process.EG_1e28.HLTPathsPrescales  = cms.vuint32(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
process.EG_1e28.HLTOverallPrescale = cms.uint32(1)
process.EG_1e28.throw = False
process.EG_1e28.andOr = True


process.load("RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi")
process.hfRecoEcalCandidate.intercept2DCut=0.3


process.z1lPath = cms.Path(
    process.patDefaultSequence *
    process.demoBefCuts *
    #  process.genFilter *
    process.EG_1e28 *
    process.z1legFilter *
    process.hfRecoEcalCandidate *
    process.demo *
    process.demoLoose
    )
