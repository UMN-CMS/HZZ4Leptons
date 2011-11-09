#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "ZShape/ZFromData/interface/ZFromData.h"
#include "ZShape/ZFromData/interface/MCAnalyTrial.h"
#include "ZShape/ZFromData/interface/InclusiveCheck.h"
#include "ZShape/ZFromData/interface/BackgroundSkimmer.h"
#include "ZShape/ZFromData/interface/IsoIdFilter.h"


DEFINE_FWK_MODULE(ZFromData);
DEFINE_FWK_MODULE(MCAnalyTrial);
DEFINE_FWK_MODULE(InclusiveCheck);
DEFINE_FWK_MODULE(BackgroundSkimmer);
DEFINE_FWK_MODULE(IsoIdFilter);


// Where are the seals!?
