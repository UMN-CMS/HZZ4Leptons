#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "ZShape/ZFromData/interface/ZFromData.h"
#include "ZShape/ZFromData/interface/MCAnalyTrial.h"
#include "ZShape/ZFromData/interface/InclusiveCheck.h"
#include "ZShape/ZFromData/interface/BackgroundSkimmer.h"


DEFINE_SEAL_MODULE();

DEFINE_ANOTHER_FWK_MODULE(ZFromData);
DEFINE_ANOTHER_FWK_MODULE(MCAnalyTrial);
DEFINE_ANOTHER_FWK_MODULE(InclusiveCheck);
DEFINE_ANOTHER_FWK_MODULE(BackgroundSkimmer);

