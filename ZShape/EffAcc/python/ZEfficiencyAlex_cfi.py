import FWCore.ParameterSet.Config as cms


mcEff = cms.EDAnalyzer("ZEfficiencyCalculator",
        Effs = cms.untracked.VPSet(
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/HLT80-PtEta.txt'),
                name = cms.untracked.string('HLT-EtaDet')
                ), 
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/HLT80-PtEta.txt'),
                name = cms.untracked.string('HLT-GSF')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId-Pt.txt'),
                name = cms.untracked.string('ElectronId-EtaDet')
                ), 
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/Iso-PtEta.txt'),
                name = cms.untracked.string('Iso-Pt')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId-Pt.txt'),
                name = cms.untracked.string('EID95')
                ), 
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/Iso-PtEta.txt'),
                name = cms.untracked.string('ISO95')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/ElectronId80-Pt.txt'),
                name = cms.untracked.string('EID80')
                ), 
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/Iso80-PtEta.txt'),
                name = cms.untracked.string('ISO80')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/GsfTrack-PtEta.txt'),
                name = cms.untracked.string('GsfTrack-EtaDet')
                ), 
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/Supercluster-Eta.txt'),
                name = cms.untracked.string('Supercluster-Eta')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/HFEId-PtEta.txt'), 
                name = cms.untracked.string('HFElectronId-EtaDet')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/HFEIdT-Pt.txt'), 
                name = cms.untracked.string('HFTightElectronId-EtaDet')
                ),
            cms.PSet(
                effFile = cms.FileInPath('ZShape/EffAcc/data/HFSC-Et_910.txt'), 
                name = cms.untracked.string('HFSuperCluster-Et')
                ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/NTEId-EtaDet.txt'),
                    name = cms.untracked.string('NTLooseElectronId-EtaDet')
                    ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/NTEIdT-EtaDet.txt'),
                    name = cms.untracked.string('NTTightElectronId-EtaDet')
                    ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/WP95-PtEta.txt'),
                    name = cms.untracked.string('WP95')
                    ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/WP90-PtEta.txt'),
                    name = cms.untracked.string('WP90')
                    ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/WP85-PtEta.txt'),
                    name = cms.untracked.string('WP85')
                    ),
            cms.PSet(
                    effFile = cms.FileInPath('ZShape/EffAcc/data/WP80-PtEta.txt'),
                    name = cms.untracked.string('WP80')
                    )                                    
            ),
    writeHistoBeforeEndJob = cms.untracked.bool(False),
    zsrc = cms.untracked.InputTag("SmearedElectronsProducer","ZEventParticles"),
    zElectronsCollection = cms.untracked.InputTag("SmearedElectronsProducer","ZEventParticles"),
    zTreeLevelElectronsCollection = cms.untracked.InputTag("ZIntoElectronsEventProducer","ZEventEle3"), quiet = cms.untracked.bool(True),
    zElectronsCone = cms.double(0.00),
    acceptanceMassWindow = cms.vdouble(60.0,120.0),
    # options include : ECALScale+,ECALScale-,HFScale+,HFScale-
    systematic = cms.untracked.string(""),                  
    ZDefs = cms.untracked.VPSet(
            cms.PSet(
                name = cms.untracked.string('ECAL_Efficiency_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","Supercluster-Eta", "GsfTrack-EtaDet","WP80"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('HF_Efficiency_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(HF)","PT20","HFSuperCluster-Et", "HFElectronId-EtaDet"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('Trackless_Efficiency_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","Supercluster-Eta","NTTightElectronId-EtaDet"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('ECAL_Background_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","!WP95"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('HF_Background_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(HF)","PT20","!HFElectronId-EtaDet"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('Trackless_Background_Chain'),
                e1 = cms.untracked.vstring( "WP80","WP80","WP80"),
                e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","!NTLooseElectronId-EtaDet"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                name = cms.untracked.string('ECAL_Background_90_Chain'),
                e1 = cms.untracked.vstring( "WP90","WP90","WP90","WP90"),
                e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT20","!GsfTrack-EtaDet","!WP95"),
                Z = cms.untracked.vstring('m(60,120)'),
                ptorder = cms.untracked.bool (True)
                ),
            cms.PSet(
                    name = cms.untracked.string('HF_Background_90_Chain'),
                    e1 = cms.untracked.vstring( "WP90","WP90","WP90"),
                    e2 = cms.untracked.vstring("ACC(HF)","PT20","!HFElectronId-EtaDet"),
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),
            cms.PSet(
                    name = cms.untracked.string('Trackless_Background_90_Chain'),
                    e1 = cms.untracked.vstring( "WP90","WP90","WP90"),
                    e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT20","!NTLooseElectronId-EtaDet"),
                    Z = cms.untracked.vstring('m(60,120)'),
                    ptorder = cms.untracked.bool (True)
                    ),


            )
)

       					
## 	cms.PSet(
##             name = cms.untracked.string('ECAL95-ECAL95'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet( 
##             name = cms.untracked.string('ECAL80-HF'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFElectronId-EtaDet"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (False)
## 	),
## 	cms.PSet(
##             name = cms.untracked.string('ECAL80-HFTight'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (False)
##         ),
## 	cms.PSet(
##             name = cms.untracked.string('ECAL80-NTLoose'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTLooseElectronId-EtaDet"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (False)
## 	),
##         cms.PSet(
##             name = cms.untracked.string('ECAL80-ECAL95'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##           ),
##         cms.PSet(
##             name = cms.untracked.string('ECAL80-ECAL95-MUO'),
##             e1 = cms.untracked.vstring("ACC(ECAL+MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL+MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##           ),
##         cms.PSet(
##             name = cms.untracked.string('ECAL80-ECAL95-MUO-HIGHPT'),
##             e1 = cms.untracked.vstring("ACC(ECAL+MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL+MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)','q(90,600)'),
##             ptorder = cms.untracked.bool (True)
##           ),
##         cms.PSet(
##             name = cms.untracked.string('MUON-MUON'),
##             e1 = cms.untracked.vstring("ACC(MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(MUON)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##           ),
##         cms.PSet(
##             name = cms.untracked.string('MUON-MUON-G20'),
##             e1 = cms.untracked.vstring("ACC(MUON)","GPT20","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(MUON)","GPT20","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP95"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##           ),
## 	cms.PSet(
##             name = cms.untracked.string('ECAL80-NTTight'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL-TRK)","PT10","Supercluster-Eta","PT20","NTTightElectronId-EtaDet"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (False)
## 	),
## 	cms.PSet(
##             name = cms.untracked.string('ECAL80-ECAL80'),
##             e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##             e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                        "GsfTrack-EtaDet","PT20","WP80"),		
##             Z = cms.untracked.vstring('m(60,120)'),
##             ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                  name = cms.untracked.string('ECAL80-HFTight-20'),
##                  e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                             "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
##                  e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
##                  Z = cms.untracked.vstring('m(60,120)'),
##                  ptorder = cms.untracked.bool (False)
##         ),
##         cms.PSet(
##                  name = cms.untracked.string('ECAL80-HFTight-25'),
##                  e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
##                                             "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
##                  e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
##                  Z = cms.untracked.vstring('m(60,120)'),
##                  ptorder = cms.untracked.bool (False)
##         ),
##         cms.PSet(
##                  name = cms.untracked.string('ECAL80-HFTight-30'),
##                  e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
##                                             "GsfTrack-EtaDet","ISO80","EID80", "HLT-GSF"),
##                  e2 = cms.untracked.vstring("ACC(HF)","PT10","HFSuperCluster-Et","PT20","HFTightElectronId-EtaDet"),
##                  Z = cms.untracked.vstring('m(60,120)'),
##                  ptorder = cms.untracked.bool (False)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95-ECAL95-20'),
##                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","ISO95","EID95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95-ECAL95-25'),
##                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
##                                            "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT25",
##                                            "GsfTrack-EtaDet","ISO95","EID95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95-ECAL95-30'),
##                 e1 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
##                                            "GsfTrack-EtaDet","ISO95","EID95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(ECAL+TRK)","PT10","Supercluster-Eta","PT30",
##                                            "GsfTrack-EtaDet","ISO95","EID95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),           
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95EE-ECAL95EE'),
##                 e1 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95EB-ECAL95EB'),
##                 e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL95EB-ECAL95EE'),
##                 e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (False)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL80EE-ECAL80EE'),
##                 e1 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL80EB-ECAL80EB'),
##                 e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL80EB-ECAL80EE'),
##                 e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP80"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (False)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL80EE-ECAL95EE'),
##                 e1 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EE+TRK)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),
##         cms.PSet(
##                 name = cms.untracked.string('ECAL80EB-ECAL95EB'),
##                 e1 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","WP95","WP80", "HLT-GSF"),
##                 e2 = cms.untracked.vstring("ACC(EB)","PT10","Supercluster-Eta","PT20",
##                                            "GsfTrack-EtaDet","PT20","WP95"),
##                 Z = cms.untracked.vstring('m(60,120)'),
##                 ptorder = cms.untracked.bool (True)
##         ),                         
       
