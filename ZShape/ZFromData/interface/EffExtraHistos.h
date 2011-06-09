#ifndef EFFEXHS_INC
#define EFFEXHS_INC

#include <TH1.h>
#include <TH2.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "ZShape/EffAcc/interface/WgtProducer.h"
// holds histograms for Z and e+/- quantities
// filled to carachterize candidates' pupulation at every step along selection path

class EffExtraHistos {

 public:
  void Book(TFileDirectory& tfd);
  void Fill(const ZShapeElectron& e1, const ZShapeElectron& e2,const ZShapeElectron& em1, const ZShapeElectron& em2,double wgt, bool doMC);

 private:

  TH1 *sc_mZ_, *sc_Y_, *sc_pt_, *sc_ptMon_;
  
  TH1 *DmZ_,*DYZ_, *DptZ_;
  TH1 *De1eta_,*De2eta_,*De1phi_,*De2phi_,*De1pt_,*De2pt_;

  TH1 *MCmZ_,*MCYZ_, *MCptZ_,*atZ_,*alZ_ ,*ptZ_,*MCatZ_,*MCalZ_;

  TH1 *MCe1eta_,*MCe2eta_,*MCe1phi_,*MCe2phi_,*MCe1pt_,*MCe2pt_,*MCeeta_,*MCephi_;
  TH2 *wgt_;
  
}; 

#endif
