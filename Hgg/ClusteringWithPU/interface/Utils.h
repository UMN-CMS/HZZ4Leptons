#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#define PI    3.14159
#define TWOPI 6.28318539

#define R_ECAL           136.5  // radius of maximum containement
#define Z_Endcap         328.0
#define etaBarrelEndcap  1.479
#define eeCrySize        2.5



float signum(float x);
float getPhiSize(reco::SuperClusterCollection::const_iterator scIt , float isBarrel );
float getPhiDistance(DetId theSeed, DetId theCry );
float getEtaDistance(DetId theSeed, DetId theCry );
