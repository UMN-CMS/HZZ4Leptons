import FWCore.ParameterSet.Config as cms

#--- This configuration was updated to include the following ---#
#--- Run 2011A: prompt V6                                    ---#
#--- Run 2011B: prompt V1 up to run 178479                   ---#
#--- ONLY lumi sections with HLT_Mu40_eta2p1 unprescaled     ---#
#--- There are more prompt V6 events, but not here...        ---#
lumisToProcess = cms.untracked.VLuminosityBlockRange()
lumisToProcess.extend([
        #--- Prompt V6 with HLT_Mu40_eta2p1 ---#
	'173236:128-173236:231',
	'173240:1-173240:124',
	'173241:1-173241:85',
	'173241:92-173241:381',
	'173241:388-173241:759',
	'173243:1-173243:89',
	'173380:75-173380:204',
	'173380:209-173380:209',
	'173380:217-173380:217',
	'173381:1-173381:134',
	'173381:139-173381:221',
	'173381:226-173381:226',
	'173381:231-173381:294',
	'173389:18-173389:106',
	'173389:111-173389:169',
	'173389:174-173389:224',
	'173389:230-173389:262',
	'173389:276-173389:653',
	'173406:37-173406:286',
	'173430:72-173430:163',
	'173438:32-173438:56',
	'173439:1-173439:274',
	'173439:280-173439:313',
	'173439:319-173439:467',
	'173439:472-173439:756',
	'173657:59-173657:93',
	'173658:1-173658:110',
	'173659:1-173659:61',
	'173659:66-173659:317',
	'173660:1-173660:81',
	'173660:86-173660:164',
	'173660:169-173660:357',
	'173660:359-173660:362',
	'173661:1-173661:56',
	'173661:62-173661:63',
	'173661:65-173661:125',
	'173663:1-173663:71',
	'173663:73-173663:120',
	'173664:5-173664:16',
	'173692:50-173692:143',
	'173692:149-173692:468',
	'173692:474-173692:1123',
	'173692:1129-173692:1162',
	'173692:1168-173692:1321',
	'173692:1326-173692:1361',
	'173692:1366-173692:1383',
	'173692:1389-173692:1431',
	'173692:1436-173692:1818',
	'173692:1824-173692:2361',
	'173692:2366-173692:2710',
	'173692:2715-173692:2754',
        #--- Run 2011B Prompt V1 ---#
	'175860:1-175860:82',
	'175863:1-175863:75',
	'175865:1-175865:58',
	'175866:1-175866:353',
	'175866:356-175866:511',
	'175872:60-175872:90',
	'175873:1-175873:32',
	'175874:1-175874:33',
	'175874:36-175874:108',
	'175877:1-175877:69',
	'175881:1-175881:160',
	'175886:1-175886:280',
	'175887:1-175887:18',
	'175887:21-175887:138',
	'175888:1-175888:424',
	'175906:79-175906:191',
	'175910:1-175910:12',
	'175921:64-175921:428',
	'175973:126-175973:183',
	'175973:186-175973:268',
	'175974:1-175974:124',
	'175975:1-175975:90',
	'175975:93-175975:101',
	'175975:104-175975:390',
	'175976:1-175976:33',
	'175976:36-175976:51',
	'175976:54-175976:102',
	'175990:56-175990:164',
	'176023:44-176023:129',
	'176161:6-176161:7',
	'176161:13-176161:14',
	'176161:16-176161:23',
	'176163:1-176163:4',
	'176163:7-176163:63',
	'176165:1-176165:29',
	'176167:1-176167:5',
	'176169:1-176169:140',
	'176201:52-176201:52',
	'176201:54-176201:55',
	'176201:58-176201:58',
	'176201:61-176201:63',
	'176201:67-176201:263',
	'176201:265-176201:570',
	'176201:573-176201:656',
	'176202:1-176202:179',
	'176202:182-176202:264',
	'176206:1-176206:101',
	'176207:1-176207:230',
	'176286:55-176286:224',
	'176286:226-176286:247',
	'176286:250-176286:252',
	'176286:254-176286:486',
	'176289:1-176289:43',
	'176289:46-176289:75',
	'176304:69-176304:69',
	'176304:73-176304:73',
	'176304:77-176304:77',
	'176304:79-176304:132',
	'176304:135-176304:422',
	'176304:425-176304:516',
	'176308:1-176308:45',
	'176308:48-176308:130',
	'176308:133-176308:293',
	'176309:1-176309:84',
	'176309:87-176309:122',
	'176309:125-176309:354',
	'176309:356-176309:467',
	'176309:470-176309:512',
	'176309:514-176309:627',
	'176309:630-176309:1177',
	'176309:1180-176309:1641',
	'176467:1-176467:85',
	'176467:88-176467:164',
	'176468:1-176468:124',
	'176468:126-176468:210',
	'176469:1-176469:5',
	'176545:5-176545:5',
	'176547:1-176547:151',
	'176547:154-176547:172',
	'176547:175-176547:184',
	'176548:1-176548:72',
	'176548:75-176548:387',
	'176548:390-176548:487',
	'176548:489-176548:494',
	'176548:496-176548:575',
	'176548:577-176548:1024',
	'176697:50-176697:66',
	'176697:68-176697:344',
	'176701:1-176701:145',
	'176702:1-176702:4',
	'176702:6-176702:42',
	'176702:45-176702:103',
	'176702:105-176702:284',
	'176702:287-176702:291',
	'176702:294-176702:351',
	'176702:353-176702:599',
	'176702:622-176702:684',
	'176702:686-176702:723',
	'176702:726-176702:794',
	'176702:796-176702:799',
	'176702:801-176702:839',
	'176765:84-176765:141',
	'176771:57-176771:60',
	'176771:62-176771:156',
	'176795:38-176795:105',
	'176796:1-176796:24',
	'176796:31-176796:100',
	'176797:1-176797:26',
	'176797:33-176797:89',
	'176797:95-176797:131',
	'176797:136-176797:196',
	'176797:203-176797:237',
	'176797:244-176797:252',
	'176799:1-176799:27',
	'176799:30-176799:35',
	'176799:38-176799:238',
	'176801:1-176801:92',
	'176801:95-176801:173',
	'176801:176-176801:177',
	'176801:179-176801:249',
	'176805:1-176805:97',
	'176807:1-176807:74',
	'176807:76-176807:255',
	'176807:257-176807:461',
	'176841:105-176841:157',
	'176841:159-176841:182',
	'176841:184-176841:235',
	'176841:238-176841:241',
	'176844:1-176844:110',
	'176848:1-176848:25',
	'176850:1-176850:29',
	'176850:32-176850:301',
	'176860:1-176860:78',
	'176860:81-176860:271',
	'176860:273-176860:335',
	'176860:338-176860:450',
	'176868:1-176868:495',
	'176885:67-176885:152',
	'176885:164-176885:235',
	'176886:1-176886:335',
	'176886:337-176886:416',
	'176886:418-176886:425',
	'176886:427-176886:429',
	'176886:431-176886:499',
	'176886:502-176886:503',
	'176886:505-176886:505',
	'176886:507-176886:507',
	'176886:509-176886:514',
	'176886:516-176886:678',
	'176889:1-176889:81',
	'176889:83-176889:102',
	'176928:1-176928:99',
	'176928:102-176928:121',
	'176929:1-176929:59',
	'176929:62-176929:198',
	'176929:201-176929:218',
	'176929:221-176929:320',
	'176933:1-176933:71',
	'176933:73-176933:358',
	'176933:361-176933:415',
	'176933:418-176933:438',
	'176959:17-176959:39',
	'176982:1-176982:10',
	'177053:54-177053:72',
	'177053:74-177053:87',
	'177053:90-177053:208',
	'177053:211-177053:217',
	'177053:219-177053:258',
	'177053:260-177053:380',
	'177053:385-177053:427',
	'177053:429-177053:570',
	'177053:573-177053:709',
	'177053:712-177053:740',
	'177074:52-177074:120',
	'177074:123-177074:468',
	'177074:470-177074:554',
	'177074:556-177074:735',
	'177074:737-177074:769',
	'177088:57-177088:83',
	'177088:85-177088:95',
	'177095:49-177095:111',
	'177095:113-177095:204',
	'177096:1-177096:14',
	'177096:17-177096:190',
	'177131:49-177131:139',
	'177138:55-177138:57',
	'177138:59-177138:114',
	'177139:1-177139:206',
	'177139:208-177139:249',
	'177139:252-177139:335',
	'177139:337-177139:375',
	'177139:378-177139:627',
	'177140:1-177140:106',
	'177140:108-177140:492',
	'177140:495-177140:557',
	'177141:1-177141:57',
	'177141:59-177141:154',
	'177141:157-177141:288',
	'177141:290-177141:492',
	'177141:494-177141:628',
	'177183:4-177183:74',
	'177183:76-177183:199',
	'177183:202-177183:219',
	'177184:1-177184:12',
	'177184:15-177184:41',
	'177201:60-177201:65',
	'177201:68-177201:73',
	'177201:75-177201:212',
	'177201:215-177201:334',
	'177201:337-177201:516',
	'177201:519-177201:603',
	'177201:606-177201:621',
	'177201:623-177201:799',
	'177222:50-177222:342',
	'177317:1-177317:60',
	'177318:1-177318:81',
	'177319:1-177319:54',
	'177319:56-177319:224',
	'177319:227-177319:286',
	'177319:289-177319:351',
	'177319:354-177319:368',
	'177449:56-177449:80',
	'177449:83-177449:234',
	'177449:236-177449:392',
	'177449:394-177449:434',
	'177452:1-177452:156',
	'177452:158-177452:169',
	'177718:53-177718:62',
	'177718:64-177718:482',
	'177718:485-177718:656',
	'177718:659-177718:697',
	'177718:699-177718:701',
	'177718:703-177718:703',
	'177718:706-177718:1572',
	'177718:1575-177718:1635',
	'177719:1-177719:194',
	'177719:197-177719:694',
	'177719:696-177719:704',
	'177730:49-177730:88',
	'177730:91-177730:234',
	'177730:237-177730:448',
	'177730:451-177730:701',
	'177730:704-177730:773',
	'177730:776-177730:931',
	'177730:933-177730:1470',
	'177730:1473-177730:1763',
	'177730:1766-177730:1956',
	'177730:1959-177730:2266',
	'177730:2268-177730:2303',
	'177776:73-177776:108',
	'177782:72-177782:75',
	'177782:78-177782:226',
	'177783:1-177783:94',
	'177783:97-177783:176',
	'177783:178-177783:234',
	'177783:237-177783:363',
	'177788:1-177788:46',
	'177789:1-177789:38',
	'177789:41-177789:82',
	'177789:84-177789:86',
	'177790:1-177790:189',
	'177790:191-177790:715',
	'177791:1-177791:27',
	'177791:30-177791:109',
	'177791:137-177791:198',
	'177875:71-177875:214',
	'177875:217-177875:367',
	'177875:370-177875:448',
	'177875:451-177875:589',
	'177878:1-177878:138',
	'177878:143-177878:492',
	'177878:495-177878:837',
	'177878:839-177878:849',
	'178098:16-178098:120',
	'178098:122-178098:142',
	'178098:144-178098:152',
	'178098:155-178098:179',
	'178098:181-178098:249',
	'178098:251-178098:549',
	'178098:551-178098:576',
	'178098:579-178098:619',
	'178098:622-178098:677',
	'178098:679-178098:707',
	'178099:1-178099:192',
	'178100:1-178100:40',
	'178100:43-178100:163',
	'178100:166-178100:511',
	'178100:514-178100:529',
	'178100:532-178100:1114',
	'178100:1117-178100:1613',
	'178101:1-178101:61',
	'178101:64-178101:75',
	'178102:1-178102:13',
	'178110:51-178110:133',
	'178110:135-178110:197',
	'178110:199-178110:246',
	'178110:252-178110:396',
	'178110:401-178110:468',
	'178116:1-178116:120',
	'178116:122-178116:615',
	'178151:47-178151:122',
	'178160:51-178160:456',
	'178160:458-178160:534',
	'178162:1-178162:4',
	'178162:7-178162:55',
	'178162:58-178162:146',
	'178162:149-178162:175',
	'178365:43-178365:96',
	'178365:98-178365:238',
	'178365:241-178365:380',
	'178365:382-178365:558',
	'178365:560-178365:660',
	'178365:663-178365:1050',
	'178365:1053-178365:1059',
	'178367:22-178367:135',
	'178367:138-178367:161',
	'178367:163-178367:264',
	'178367:266-178367:316',
	'178367:318-178367:366',
	'178367:369-178367:392',
	'178367:395-178367:548',
	'178367:560-178367:711',
	'178367:713-178367:737',
	'178380:1-178380:81',
	'178380:84-178380:195',
	'178420:53-178420:61',
	'178420:63-178420:70',
	'178421:1-178421:3',
	'178421:47-178421:77',
	'178421:79-178421:135',
	'178421:137-178421:149',
	'178421:151-178421:265',
	'178421:268-178421:293',
	'178421:295-178421:320',
	'178421:323-178421:512',
	'178421:514-178421:586',
	'178421:588-178421:742',
	'178421:744-178421:756',
	'178421:759-178421:915',
	'178421:917-178421:1207',
	'178424:1-178424:536',
	'178424:538-178424:594',
	'178424:597-178424:837',
	'178479:132-178479:245',
	'178479:247-178479:354',
	'178479:356-178479:576',
	'178479:578-178479:985',
	'178479:988-178479:1004'
])