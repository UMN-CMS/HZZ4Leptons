#include "ZShape/Base/interface/ZShapeEvent.h"
#include <DataFormats/Math/interface/LorentzVector.h>
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include <iostream>

ZShapeElectron* ZShapeEvent::allocateElectron() {
  return new ZShapeElectron(); // vanilla
}

reco::GenParticle* ZShapeEvent::allocateTreeLevelElectron() {
  return new reco::GenParticle(); // vanilla(?)
}

void ZShapeEvent::dump() const {
  dump(std::cout);
}

void ZShapeEvent::dump(std::ostream& s) const {
  s << "Vertex: " << vtx_.x() << ", " << vtx_.y() << ", " << vtx_.z() << std::endl;
  for (int i=0; i<2; i++) {
    s << "Electron " << i+1;
    elec(i).dump(s);
  }
}

void ZShapeEvent::clear() {
  if (e1_==0) e1_=allocateElectron();
  if (e2_==0) e2_=allocateElectron();
  elec(0).clear();
  elec(1).clear();
  n_elec=0;
  n_gsf20=0;

  if(eTL1_==0) eTL1_=allocateTreeLevelElectron();
  if(eTL2_==0) eTL2_=allocateTreeLevelElectron();
  n_TLelec=0;

  // initializing
  eTL1_->setP4(reco::Particle::LorentzVector(0,0,0,0));
  eTL1_->setCharge(0);
  eTL1_->setVertex(math::XYZPoint(0,0,0));
  eTL1_->setPdgId(0);
  eTL1_->setStatus(0);
  eTL2_->setP4(reco::Particle::LorentzVector(0,0,0,0));
  eTL2_->setCharge(0);
  eTL2_->setVertex(math::XYZPoint(0,0,0));
  eTL2_->setPdgId(0);
  eTL2_->setStatus(0);

  free_i = 0;

}

ZShapeEvent::ZShapeEvent() : e1_(0),e2_(0),eTL1_(0),eTL2_(0) {
  npu_ = 0;
  mNPV_ = 0;
}

ZShapeEvent::~ZShapeEvent() {
  if (e1_!=0) delete e1_;
  if (e2_!=0) delete e2_;
  if (eTL1_!=0) delete eTL1_;
  if (eTL2_!=0) delete eTL2_;
}

void ZShapeEvent::afterLoad() {
  math::PtEtaPhiMLorentzVector ZP4 = elec(0).p4_ + elec(1).p4_;
  mass_=ZP4.M();
  rap_=ZP4.Rapidity();
  qT_=ZP4.Pt();
}

void ZShapeEvent::initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  /// *********** store some event variables: MET, SumET ******
  //////////// Primary vertex //////////////
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  mNPV_ = 0;
  mPVx_ =  100.0;
  mPVy_ =  100.0;
  mPVz_ =  100.0;

  for(unsigned int ind=0;ind<recVtxs->size();ind++) {
      if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>4) 
              && (fabs((*recVtxs)[ind].z())<=24.0) &&  
              ((*recVtxs)[ind].position().Rho()<=2.0) ) {
          mNPV_++;
          if(mNPV_==1) { // store the first good primary vertex
              mPVx_ = (*recVtxs)[ind].x();
              mPVy_ = (*recVtxs)[ind].y();
              mPVz_ = (*recVtxs)[ind].z();
          }
      }
  }

  // Pile up counting
  if (!iEvent.isRealData()){
      edm::Handle<std::vector<PileupSummaryInfo> > pPU;
      iEvent.getByLabel("addPileupInfo", pPU);
      if (pPU.isValid()){
          npu_ = pPU->size();
      } else {
          npu_ = -1; 
      }
  }

  //////////// Beam spot //////////////
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot);
  mBSx_ = beamSpot->position().X();
  mBSy_ = beamSpot->position().Y();
  mBSz_ = beamSpot->position().Z();


  ////////////// CaloMET //////
  edm::Handle<reco::CaloMETCollection> met;
  iEvent.getByLabel("met",met);
  if (met->size() == 0) {
    mMET_   = -1;
    mSumET_ = -1;
    mMETSign_ = -1;
  }
  else {
    mMET_   = (*met)[0].et();
    mSumET_ = (*met)[0].sumEt();
    mMETSign_ = (*met)[0].significance();
  }

  /////// TcMET information /////
  edm::Handle<reco::METCollection> tcmet;
  iEvent.getByLabel("tcMet", tcmet);
  if (tcmet->size() == 0) {
    mtcMET_   = -1;
    mtcSumET_ = -1;
    mtcMETSign_ = -1;
  }
  else {
    mtcMET_   = (*tcmet)[0].et();
    mtcSumET_ = (*tcmet)[0].sumEt();
    mtcMETSign_ = (*tcmet)[0].significance();
  }

  /////// PfMET information /////
  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByLabel("pfMet", pfmet);
  if (pfmet->size() == 0) {
    mpfMET_   = -1;
    mpfSumET_ = -1;
    mpfMETSign_ = -1;
  }
  else {
    mpfMET_   = (*pfmet)[0].et();
    mpfSumET_ = (*pfmet)[0].sumEt();
    mpfMETSign_ = (*pfmet)[0].significance();
  }

}

