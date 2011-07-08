import FWCore.ParameterSet.Config as cms

SCwithPUData = cms.EDAnalyzer('SCwithPUData',
                              DoLog             = cms.bool(False),
                              myName            = cms.string('nameOfSCwithPUData'),
                              eleWorkingPoint   = cms.string('simpleEleId80relIso'),
                              vertexCollection  = cms.InputTag("offlinePrimaryVertices",""),
                              patElectrons      = cms.InputTag("patElectrons",""),
                              ETCut             = cms.double(20),
                              etaMin            = cms.double(0),       # limit to EB clusters for now
                              etaMax            = cms.double(1.4442),
                              #etaMax            = cms.double(9999),
                              acceptedElectronIDs = cms.vint32( 7 ),
                              xi                = cms.double( 0.0 )
                              )
