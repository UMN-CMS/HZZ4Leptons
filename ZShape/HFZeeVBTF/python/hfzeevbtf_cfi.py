import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HFZeeVBTF',
# source
  acceptedElectronIDs = cms.vint32( 7 )
)
