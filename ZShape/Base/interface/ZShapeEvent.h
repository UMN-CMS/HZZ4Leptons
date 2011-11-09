#ifndef ZSHAPE_BASE_ZSHAPEEVENT_H
#define ZSHAPE_BASE_ZSHAPEEVENT_H 1

#include "DataFormats/Math/interface/Vector3D.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <map>
#include <string>
#include <ostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"


/** \class ZShapeEvent
  *  
  * $Date: 2011/02/28 14:17:01 $
  * $Revision: 1.8 $
  * \author J. Mans - Minnesota
  */
class ZShapeEvent {
public:
  ZShapeEvent();
  // primary vertex
  ::math::XYZPoint vtx_;
  int n_elec;
  int n_TLelec;
  int n_gsf20;

  double m() const { return mass_; }
  double Y() const { return rap_; }
  double qT() const { return qT_; }
  float EvtPVz() const { return mPVz_; }
  float EvtBSz() const { return mBSz_; }
  float EvtMET() const { return mMET_; }
  float EvtTCMET() const { return mtcMET_; }
  float EvtPFMET() const { return mpfMET_; }
  uint32_t nvert() const { return mNPV_; }
  int npu() const {return npu_;}


  ZShapeElectron& elec(int i) { return ((i==0)?(*e1_):(*e2_)); }
  const ZShapeElectron& elec(int i) const { return ((i==0)?(*e1_):(*e2_)); }
  reco::GenParticle& elecTreeLevel(int i) { return ((i==0)?(*eTL1_):(*eTL2_)); }
  const reco::GenParticle& elecTreeLevel(int i) const { return ((i==0)?(*eTL1_):(*eTL2_)); }
  
  //This loads the Extra event information
  void initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  virtual ~ZShapeEvent();
  virtual void clear();
  virtual void dump(std::ostream& s) const;
  virtual void dump() const;
  void afterLoad();
private:
  virtual ZShapeElectron* allocateElectron();
  virtual reco::GenParticle* allocateTreeLevelElectron();
  ZShapeElectron* e1_;
  ZShapeElectron* e2_;
  reco::GenParticle * eTL1_;
  reco::GenParticle * eTL2_;
  double mass_, rap_, qT_;
  //These EXTRA Event Variables are mutable as they are filled via a const method. Sign=Significance
  uint32_t run_, lumi_, event_, mNPV_;
  float mPVx_,mPVy_,mPVz_,mBSx_,mBSy_,mBSz_; 
  float mMET_,mSumET_,mMETSign_,mtcMET_,mtcSumET_,
	        mtcMETSign_,mpfMET_,mpfSumET_,mpfMETSign_;
  int npu_;

};

#endif
