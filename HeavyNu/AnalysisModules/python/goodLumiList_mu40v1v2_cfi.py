import FWCore.ParameterSet.Config as cms

lumisToProcess = cms.untracked.VLuminosityBlockRange(
	'165088:107-165088:266',
	'165098:124-165098:187',
	'165098:190-165098:193',
	'165098:195-165098:248',
	'165098:250-165098:254',
	'165098:256-165098:331',
	'165098:333-165098:367',
	'165098:369-165098:415',
	'165099:1-165099:105',
	'165102:1-165102:185',
	'165103:1-165103:440',
	'165120:82-165120:97',
	'165121:1-165121:466',
	'165205:80-165205:248',
	'165208:1-165208:101',
	'165364:45-165364:111',
	'165364:114-165364:147',
	'165364:160-165364:807',
	'165364:809-165364:1220',
	'165364:1260-165364:1301',
	'165402:1-165402:28',
	'165415:58-165415:85',
	'165415:88-165415:640',
	'165415:643-165415:707',
	'165415:712-165415:777',
	'165415:780-165415:1356',
	'165415:1360-165415:1383',
	'165467:39-165467:708',
	'165472:1-165472:184',
	'165472:186-165472:882',
	'165486:37-165486:102',
	'165487:1-165487:151',
	'165506:54-165506:170',
	'165514:72-165514:244',
	'165514:246-165514:560',
	'165514:562-165514:567',
	'165548:1-165548:363',
	'165548:365-165548:381',
	'165548:384-165548:589',
	'165558:1-165558:62',
	'165567:54-165567:109',
	'165567:114-165567:309',
	'165567:315-165567:631',
	'165570:1-165570:2',
	'165570:5-165570:83',
	'165570:88-165570:942',
	'165570:944-165570:946',
	'165617:26-165617:52',
	'165617:54-165617:143',
	'165617:145-165617:288',
	'165620:14-165620:19',
	'165633:56-165633:62',
	'165633:64-165633:64',
	'165633:66-165633:317',
	'165633:319-165633:500',
	'165970:67-165970:329',
	'165970:331-165970:335',
	'165993:71-165993:873',
	'165993:879-165993:1660',
	'165993:1665-165993:1697',
	'166011:1-166011:81',
	'166011:83-166011:83',
	'166033:35-166033:53',
	'166033:59-166033:330',
	'166033:336-166033:355',
	'166033:360-166033:444',
	'166033:450-166033:606',
	'166033:613-166033:707',
	'166033:713-166033:1233',
	'166034:1-166034:109',
	'166034:115-166034:228',
	'166034:234-166034:307',
	'166049:53-166049:86',
	'166049:88-166049:236',
	'166049:242-166049:674',
	'166149:1-166149:2',
	'166150:1-166150:99',
	'166150:101-166150:116',
	'166161:38-166161:120',
	'166161:122-166161:123',
	'166161:126-166161:126',
	'166163:1-166163:12',
	'166163:14-166163:33',
	'166164:1-166164:32',
	'166164:34-166164:40',
	'166346:48-166346:210',
	'166346:212-166346:215',
	'166374:46-166374:64',
	'166374:66-166374:188',
	'166380:1-166380:367',
	'166380:373-166380:711',
	'166380:715-166380:1400',
	'166380:1406-166380:1809',
	'166408:67-166408:283',
	'166408:291-166408:947',
	'166408:953-166408:1236',
	'166429:33-166429:89',
	'166438:32-166438:85',
	'166438:87-166438:856',
	'166438:858-166438:866',
	'166462:78-166462:102',
	'166462:108-166462:317',
	'166462:323-166462:526',
	'166486:54-166486:75',
	'166486:80-166486:95',
	'166486:97-166486:174',
	'166502:43-166502:78',
	'166502:83-166502:109',
	'166512:42-166512:430',
	'166512:432-166512:487',
	'166512:491-166512:605',
	'166512:611-166512:1279',
	'166512:1281-166512:1818',
	'166512:1821-166512:1862',
	'166512:1868-166512:1868',
	'166512:1870-166512:1871',
	'166512:1873-166512:1874',
	'166514:1-166514:455',
	'166514:460-166514:464',
	'166530:43-166530:63',
	'166554:46-166554:218',
	'166554:224-166554:287',
	'166554:290-166554:317',
	'166554:320-166554:595',
	'166554:597-166554:730',
	'166554:732-166554:734',
	'166554:736-166554:736',
	'166563:1-166563:276',
	'166563:492-166563:748',
	'166565:1-166565:147',
	'166565:153-166565:312',
	'166565:316-166565:467',
	'166565:469-166565:898',
	'166699:55-166699:234',
	'166699:240-166699:415',
	'166699:421-166699:477',
	'166699:483-166699:677',
	'166699:681-166699:912',
	'166701:1-166701:13',
	'166701:16-166701:319',
	'166701:324-166701:506',
	'166701:513-166701:551',
	'166701:557-166701:672',
	'166701:681-166701:705',
	'166701:712-166701:724',
	'166701:731-166701:757',
	'166701:764-166701:777',
	'166701:783-166701:792',
	'166763:46-166763:168',
	'166763:174-166763:650',
	'166781:41-166781:111',
	'166781:115-166781:115',
	'166781:117-166781:233',
	'166781:236-166781:253',
	'166781:255-166781:382',
	'166782:1-166782:569',
	'166784:1-166784:114',
	'166784:119-166784:276',
	'166784:281-166784:365',
	'166787:1-166787:55',
	'166787:60-166787:127',
	'166787:132-166787:364',
	'166839:43-166839:173',
	'166839:178-166839:297',
	'166839:299-166839:302',
	'166841:1-166841:845',
	'166841:851-166841:876',
	'166841:882-166841:977',
	'166841:984-166841:984',
	'166841:988-166841:992',
	'166841:998-166841:1015',
	'166842:1-166842:170',
	'166859:62-166859:418',
	'166859:421-166859:421',
	'166859:423-166859:423',
	'166860:1-166860:21',
	'166861:1-166861:6',
	'166861:8-166861:13',
	'166864:1-166864:29',
	'166864:31-166864:77',
	'166864:79-166864:99',
	'166864:102-166864:119',
	'166864:125-166864:247',
	'166864:249-166864:307',
	'166864:311-166864:365',
	'166864:367-166864:374',
	'166864:378-166864:454',
	'166864:478-166864:537',
	'166888:56-166888:90',
	'166888:93-166888:154',
	'166888:156-166888:394',
	'166888:398-166888:470',
	'166889:1-166889:73',
	'166889:79-166889:228',
	'166890:1-166890:441',
	'166894:1-166894:190',
	'166895:1-166895:66',
	'166895:72-166895:597',
	'166895:599-166895:603',
	'166911:58-166911:76',
	'166911:81-166911:104',
	'166922:1-166922:39',
	'166922:41-166922:105',
	'166922:110-166922:340',
	'166922:345-166922:418',
	'166922:423-166922:747',
	'166922:752-166922:769',
	'166922:773-166922:773',
	'166923:1-166923:382',
	'166923:389-166923:470',
	'166946:41-166946:72',
	'166946:75-166946:201',
	'166950:1-166950:1',
	'166950:8-166950:31',
	'166950:36-166950:210',
	'166950:216-166950:877',
	'166950:883-166950:950',
	'166950:956-166950:1012',
	'166950:1018-166950:1321',
	'166950:1327-166950:1345',
	'166950:1347-166950:1438',
	'166960:1-166960:137',
	'166960:143-166960:166',
	'166966:1-166966:238',
	'166967:1-166967:220',
	'167039:20-167039:92',
	'167039:98-167039:228',
	'167041:1-167041:336',
	'167041:339-167041:391',
	'167041:396-167041:462',
	'167041:467-167041:663',
	'167043:1-167043:125',
	'167043:130-167043:235',
)