#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "ZShape/ZFromData/interface/ZFromData.h"
#include "ZShape/ZFromData/interface/MCAnalyTrial.h"


DEFINE_SEAL_MODULE();

DEFINE_ANOTHER_FWK_MODULE(ZFromData);
DEFINE_ANOTHER_FWK_MODULE(MCAnalyTrial);
