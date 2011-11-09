#include "ZShape/Base/interface/ZShapeElectron.h"
#include <iostream>

void ZShapeElectron::dump(std::ostream& s) const {
  s << "(eta,phi,pt : deteta) " << p4_.eta() << ", " << p4_.phi() << ", " << p4_.Pt() << " : " << detEta_ << std::endl;
  for (std::map<std::string,bool>::const_iterator j=cutMap_.begin(); j!=cutMap_.end(); ++j) {
    s << "  \"" << j->first << "\" : " << ((j->second)?("pass"):("fail")) << std::endl;
  }
}

void ZShapeElectron::clear() {
  cutMap_.clear();
}

bool ZShapeElectron::passed(const std::string& cutName) const { // false if not present!
  std::map<std::string,bool>::const_iterator i=cutMap_.find(cutName);
  return ((i==cutMap_.end())?(false):(i->second));
}

void ZShapeElectron::cutResult(const std::string& cutName, bool didpass) {
  cutMap_[cutName]=didpass;
  cutMap_["!"+cutName]=!didpass;
}

double ZShapeElectron::detectorEta(const math::XYZPoint& vtx) const {
  const double EBr=1290+8.9*6; // front face + 6 Xo
  const double EEz=3144+8.9*6; // "front face" + 6 Xo 
  const double cutover=1.51;
  const double cutover2=3.0;
  // current version ignores any X/Y displacement of the virtex
  double de=p4_.eta();

  double aeta=fabs(de);
  
  if (aeta>cutover2) { // do nothing (HF)
  } else if (aeta>cutover) { // EE correction
    double EEz_here=(de<0)?(-EEz):(EEz);  
    de=asinh(sinh(de)*1.0/(1.0-vtx.z()/EEz_here));
  } else { // EB correction
    de=asinh(sinh(de)+vtx.z()/EBr);
  }
  return de;
}
