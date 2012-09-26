import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")


process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring(     )
                  
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )




# this is the equivalent of the JSON file 
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(


# from: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-140399_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt

## # # # # # # # # # # # # # # # # # # # # # # # #  
# * runrange 132440 - 137436 : Jun14th reprocessing

#  	'132440:157-132440:378',
#  	'132596:382-132596:382',
#  	'132596:447-132596:447',
#  	'132598:174-132598:175',
#  	'132599:1-132599:379',
#  	'132599:381-132599:437',
#  	'132601:1-132601:207',
#  	'132601:209-132601:259',
#  	'132601:261-132601:1107',
#  	'132602:1-132602:70',
#  	'132605:1-132605:444',
#  	'132605:446-132605:522',
#  	'132605:526-132605:622',
#  	'132605:624-132605:814',
#  	'132605:816-132605:829',
#  	'132605:831-132605:867',
#  	'132605:896-132605:942',
#  	'132606:1-132606:26',
#  	'132656:1-132656:111',
#  	'132658:1-132658:51',
#  	'132658:56-132658:120',
#  	'132658:127-132658:148',
#  	'132659:1-132659:75',
#  	'132661:1-132661:115',
#  	'132662:1-132662:9',
#  	'132662:25-132662:74',
#  	'132716:220-132716:436',
#  	'132716:440-132716:487',
#  	'132716:491-132716:586',
#  	'132959:326-132959:334',
#  	'132960:1-132960:124',
#  	'133029:101-133029:115',
#  	'133029:129-133029:331',
#  	'133031:1-133031:18',
#  	'133035:1-133035:63',
#  	'133035:67-133035:301',
#  	'133036:1-133036:222',
#  	'133046:1-133046:43',
#  	'133046:45-133046:210',
#  	'133046:213-133046:227',
#  	'133046:229-133046:323',
#  	'133082:1-133082:336',
#  	'133082:523-133082:592',
#  	'133082:595-133082:608',
#  	'133158:65-133158:786',
#  	'133874:166-133874:297',
#  	'133874:299-133874:721',
#  	'133874:724-133874:814',
#  	'133874:817-133874:864',
#  	'133875:1-133875:20',
#  	'133875:22-133875:37',
#  	'133876:1-133876:315',
#  	'133877:1-133877:77',
#  	'133877:82-133877:104',
#  	'133877:113-133877:231',
#  	'133877:236-133877:294',
#  	'133877:297-133877:437',
#  	'133877:439-133877:622',
#  	'133877:624-133877:853',
#  	'133877:857-133877:1472',
#  	'133877:1474-133877:1640',
#  	'133877:1643-133877:1930',
#  	'133881:1-133881:71',
#  	'133881:74-133881:223',
#  	'133881:225-133881:551',
#  	'133885:1-133885:132',
#  	'133885:134-133885:710',
#  	'133927:1-133927:44',
#  	'133928:1-133928:645',
#  	'135059:59-135059:67',
#  	'135149:297-135149:337',
#  	'135149:339-135149:754',
#  	'135149:756-135149:932',
#  	'135149:934-135149:937',
#  	'135149:942-135149:993',
#  	'135149:995-135149:1031',
#  	'135149:1033-135149:1098',
#  	'135149:1102-135149:1808',
#  	'135149:1811-135149:2269',
#  	'135149:2274-135149:2524',
#  	'135149:2528-135149:2713',
#  	'135149:2715-135149:3098',
#  	'135149:3100-135149:3102',
#  	'135149:3105-135149:3179',
#  	'135149:3182-135149:3303',
#  	'135149:3305-135149:3380',
#  	'135175:55-135175:545',
#  	'135175:548-135175:561',
#  	'135175:563-135175:790',
#  	'135175:792-135175:1046',
#  	'135521:60-135521:108',
#  	'135521:110-135521:359',
#  	'135521:361-135521:440',
#  	'135521:442-135521:488',
#  	'135523:1-135523:64',
#  	'135523:66-135523:109',
#  	'135523:113-135523:124',
#  	'135523:126-135523:211',
#  	'135525:1-135525:3',
#  	'135525:6-135525:119',
#  	'135525:121-135525:143',
#  	'135525:145-135525:275',
#  	'135525:293-135525:381',
#  	'135525:384-135525:435',
#  	'135525:437-135525:451',
#  	'135528:1-135528:91',
#  	'135528:94-135528:95',
#  	'135528:98-135528:142',
#  	'135528:145-135528:147',
#  	'135528:149-135528:308',
#  	'135528:310-135528:454',
#  	'135528:456-135528:606',
#  	'135528:608-135528:609',
#  	'135528:611-135528:770',
#  	'135528:773-135528:776',
#  	'135528:779-135528:813',
#  	'135528:816-135528:912',
#  	'135528:915-135528:924',
#  	'135528:926-135528:1082',
#  	'135528:1084-135528:1213',
#  	'135528:1215-135528:1436',
#  	'135535:1-135535:42',
#  	'135535:75-135535:167',
#  	'135535:169-135535:231',
#  	'135573:102-135573:110',
#  	'135573:113-135573:118',
#  	'135573:120-135573:154',
#  	'135575:2-135575:210',
#  	'135575:213-135575:241',
#  	'135575:243-135575:264',
#  	'135575:266-135575:381',
#  	'135575:384-135575:638',
#  	'135575:645-135575:1161',
#  	'135575:1163-135575:1253',
#  	'135735:31-135735:42',
#  	'135735:44-135735:149',
#  	'135735:151-135735:234',
#  	'135735:236-135735:320',
#  	'136066:181-136066:297',
#  	'136066:299-136066:319',
#  	'136066:321-136066:321',
#  	'136066:323-136066:336',
#  	'136066:338-136066:348',
#  	'136066:350-136066:366',
#  	'136066:368-136066:529',
#  	'136066:532-136066:544',
#  	'136066:546-136066:595',
#  	'136066:597-136066:719',
#  	'136066:721-136066:1004',
#  	'136066:1006-136066:1184',
#  	'136080:249-136080:250',
#  	'136080:252-136080:255',
#  	'136082:1-136082:173',
#  	'136082:175-136082:422',
#  	'136082:477-136082:504',
#  	'136082:506-136082:506',
#  	'136087:250-136087:265',
#  	'136087:267-136087:277',
#  	'136087:279-136087:310',
#  	'136087:312-136087:315',
#  	'136088:1-136088:3',
#  	'136088:6-136088:17',
#  	'136088:19-136088:39',
#  	'136088:41-136088:47',
#  	'136088:50-136088:50',
#  	'136088:52-136088:52',
#  	'136088:54-136088:63',
#  	'136088:65-136088:65',
#  	'136088:67-136088:89',
#  	'136088:91-136088:113',
#  	'136088:115-136088:131',
#  	'136088:133-136088:137',
#  	'136088:139-136088:184',
#  	'136088:187-136088:191',
#  	'136088:193-136088:194',
#  	'136088:196-136088:206',
#  	'136088:208-136088:208',
#  	'136088:210-136088:262',
#  	'136097:1-136097:1',
#  	'136097:3-136097:3',
#  	'136097:5-136097:5',
#  	'136097:7-136097:7',
#  	'136097:9-136097:9',
#  	'136097:11-136097:12',
#  	'136097:14-136097:78',
#  	'136098:1-136098:10',
#  	'136100:1-136100:94',
#  	'136119:1-136119:36',
#  	'137027:98-137027:146',
#  	'137027:149-137027:162',
#  	'137027:165-137027:187',
#  	'137028:1-137028:107',
#  	'137028:111-137028:157',
#  	'137028:163-137028:169',
#  	'137028:171-137028:180',
#  	'137028:184-137028:200',
#  	'137028:203-137028:209',
#  	'137028:211-137028:232',
#  	'137028:234-137028:274',
#  	'137028:276-137028:338',
#  	'137028:341-137028:475',
#  	'137028:478-137028:484'
# 


## # # # # # # # # # # # # # # # # # # # # # # # #  
#  * runrange 137437 - 139558 : prompt reco
# 
#   	'138564:1-138564:4',
#   	'138571:1-138571:13',
#   	'138572:1-138572:213',
#   	'138737:1-138737:66',
#   	'138738:1-138738:10',
#   	'138739:1-138739:13',
#   	'138742:1-138742:20',
#   	'138742:22-138742:38',
#   	'138744:1-138744:24',
#   	'138745:1-138745:10',
#   	'138746:1-138746:129',
#   	'138747:1-138747:71',
#   	'138747:73-138747:131',
#   	'138750:1-138750:46',
#   	'138750:49-138750:208',
#   	'138750:210-138750:623',
#   	'138750:626-138750:709',
#   	'138751:1-138751:110',
#   	'138751:112-138751:147',
#   	'138919:62-138919:150',
#   	'138920:1-138920:55',
#   	'138921:1-138921:165',
#   	'138921:168-138921:181',
#   	'138921:183-138921:191',
#   	'138923:1-138923:2',
#   	'138924:1-138924:62',
#   	'138937:1-138937:26',
#   	'138939:1-138939:26',
#   	'139020:227-139020:316',
#   	'139020:319-139020:617',
#   	'139096:193-139096:210',
#   	'139096:212-139096:233',
#   	'139096:235-139096:236',
#   	'139096:239-139096:239',
#   	'139096:242-139096:267',
#   	'139098:1-139098:92',
#   	'139098:94-139098:121',
#   	'139098:123-139098:160',
#   	'139098:162-139098:179',
#   	'139100:1-139100:13',
#   	'139100:15-139100:15',
#   	'139100:17-139100:19',
#   	'139100:21-139100:102',
#   	'139100:104-139100:104',
#   	'139100:106-139100:205',
#   	'139100:207-139100:307',
#   	'139102:1-139102:12',
#   	'139102:14-139102:55',
#   	'139103:7-139103:154',
#   	'139103:157-139103:414',
#   	'139103:416-139103:416',
#   	'139103:418-139103:449',
#   	'139195:9-139195:14',
#   	'139195:16-139195:29',
#   	'139195:31-139195:81',
#   	'139195:83-139195:92',
#   	'139239:164-139239:214',
#   	'139239:217-139239:217',
#   	'139239:219-139239:222',
#   	'139239:224-139239:225',
#   	'139239:227-139239:270',
#   	'139347:214-139347:357',
#   	'139347:359-139347:516',
#   	'139356:175-139356:186',
#   	'139356:188-139356:192',
#   	'139360:1-139360:30',
#   	'139364:1-139364:48',
#   	'139364:51-139364:72',
#   	'139365:1-139365:9',
#   	'139365:12-139365:14',
#   	'139365:16-139365:82',
#   	'139365:84-139365:84',
#   	'139365:86-139365:112',
#   	'139365:114-139365:166',
#   	'139365:168-139365:179',
#   	'139365:181-139365:250',
#   	'139368:1-139368:34',
#   	'139368:36-139368:100',
#   	'139368:102-139368:150',
#   	'139368:152-139368:166',
#   	'139370:1-139370:175',
#   	'139370:177-139370:254',
#   	'139370:257-139370:445',
#   	'139370:447-139370:468',
#   	'139370:472-139370:504',
#   	'139370:506-139370:562',
#   	'139370:564-139370:581',
#   	'139370:583-139370:619',
#   	'139370:622-139370:645',
#   	'139372:20-139372:46',
#   	'139372:48-139372:173',
#   	'139372:175-139372:194',
#   	'139375:1-139375:47',
#   	'139399:75-139399:89',
#   	'139399:94-139399:107',
#   	'139399:109-139399:109',
#   	'139399:111-139399:125',
#   	'139400:1-139400:22',
#   	'139407:1-139407:976',
#   	'139407:978-139407:1074',
#   	'139407:1076-139407:1251',
#   	'139411:4-139411:67',
#   	'139411:69-139411:135',
#   	'139457:1-139457:73',
#   	'139458:1-139458:103',
#   	'139458:107-139458:326',
#   	'139458:329-139458:395',
#   	'139459:1-139459:2',
#   	'139459:4-139459:55'
# 
# 
# ## # # # # # # # # # # # # # # # # # # # # # # # #  
# #  * runrange 139559 - 140174 : Jul16th reprocessing
# #  from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_139779-1340159_7TeV_July16thReReco_Collisions10_COMSSWConfig.txt
#           '139779:4-139779:52',
#                   '139780:4-139780:75',
#                   '139781:4-139781:200',
#                   '139783:4-139783:36',
#                   '139784:4-139784:71',
#                   '139786:4-139786:101',
#                   '139788:4-139788:35',
#                   '139789:4-139789:112',
#                   '139790:4-139790:21',
#                   '139965:209-139965:213',
#                   '139966:1-139966:51',
#                   '139967:1-139967:55',
#                   '139968:1-139968:5',
#                   '139969:1-139969:78',
#                   '139971:1-139971:304',
#                   '139972:1-139972:38',
#                   '139973:1-139973:118',
#                   '139974:1-139974:53',
#                   '139975:1-139975:15',
#                   '139975:17-139975:54',
#                   '140058:111-140058:224',
#                   '140059:1-140059:610',
#                   '140059:613-140059:1160',
#                   '140070:1-140070:11',
#                   '140076:1-140076:60',
#                   '140119:1-140119:14',
#                   '140123:1-140123:53',
#                   '140123:58-140123:62',
#                   '140123:67-140123:83',
#                   '140123:85-140123:86',
#                   '140123:88-140123:93',
#                   '140123:95-140123:102',
#                   '140124:1-140124:390',
#                   '140124:392-140124:938',
#                   '140124:940-140124:1733',
#                   '140126:1-140126:322',
#                   '140133:1-140133:65',
#                   '140158:5-140158:404',
#                   '140159:1-140159:379'
#   

# ## # # # # # # # # # # # # # # # # # # # # # # # #  
# #       * runrange 140175 - 999999 : prompt reco
 
   	'140160:1-140160:217',
   	'140180:1-140180:55',
   	'140181:1-140181:13',
   	'140331:222-140331:544',
   	'140331:546-140331:614',
   	'140331:626-140331:626',
   	'140331:635-140331:749',
   	'140359:1-140359:35',
   	'140359:40-140359:43',
   	'140361:1-140361:151',
   	'140362:1-140362:354',
   	'140379:123-140379:244',
   	'140381:16-140381:37',
   	'140382:1-140382:128',
   	'140382:131-140382:214',
   	'140383:1-140383:276',
   	'140383:281-140383:447',
   	'140383:453-140383:487',
   	'140385:1-140385:4',
   	'140385:26-140385:30',
   	'140385:50-140385:185',
   	'140386:1-140386:23',
   	'140387:7-140387:137'
   

)




## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
process.GlobalTag.globaltag = cms.string('GR_R_36X_V11A::All') ## GF changed here
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
 ## GF changed here: 
# you have to tell the ID that it is data
process.simpleEleId95relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId95cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60cIso.dataMagneticFieldSetUp = cms.bool(True)
#

process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons)
# process.makePatMuons may be needed depending on how you calculate the MET
process.makePatCandidates = cms.Sequence(process.makePatElectrons+process.makePatMETs)
process.patDefaultSequence = cms.Sequence(process.makePatCandidates)
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
#
HLT_path_name_extra   = "HLT_Photon15_L1R" #= "HLT_Ele15_LW_L1R" #
HLT_filter_name_extra = "hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter"

process.z1legFilter =cms.EDFilter('WenuCandidateFilter',
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
                                  dataMagneticFieldSetUp = cms.untracked.bool(True),
                                  dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
                                  )

####################################################################################
process.z1lPath = cms.Path(process.patDefaultSequence*process.z1legFilter)


process.z1lOutputModule = cms.OutputModule( "PoolOutputModule",
       fileName = cms.untracked.string(""),
                                             maxSize = cms.untracked.int32(300000),
                                             SelectEvents = cms.untracked.PSet(
                                             SelectEvents = cms.vstring('z1lPath',),
                                                                              )
                                             )

process.outpath = cms.EndPath(process.z1lOutputModule)
