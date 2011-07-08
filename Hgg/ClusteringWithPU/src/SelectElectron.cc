#include <algorithm>
#include "Hgg/ClusteringWithPU/interface/SelectElectron.h"

void SelectElectron::add(int v)
{ allowed.push_back(v);}

bool SelectElectron::doesElePass(int v)
{ 
  if ( find(allowed.begin(),allowed.end(),v) != allowed.end() ) return true;
  else return false;
}
