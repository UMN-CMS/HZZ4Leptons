#include "Hgg/ClusteringWithPU/interface/Utils.h"


#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


float signum(float x) {    return (x>0)?1:((x<0)?-1:0); }




// PhiSize is in unit of: crystals for EB 
//                        radiants for EE (change this to crystals too)             
float getPhiSize(reco::SuperClusterCollection::const_iterator scIt , float isBarrel ){
  // scIt is an iterator pointing to a supercluster
  // isBarrel has to be passed from outside - resolve there the clusters across EB and EE
  
  std::vector< std::pair<DetId, float> >  theHitsAndFractions =  (*scIt).hitsAndFractions();  
  
  float phiLow =0;
  float phiHigh=0;
  float phiFirstEB =-9999;  // a sort of seeding is done with the first crystal
  float phiFirstEE =-9999;  // in order to handle the boundaries where phi is discontinuos
  
  for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
      idsIt != theHitsAndFractions.end(); ++idsIt) 
    {
      
      float thePhi=-99999;
      
      if      ( isBarrel && (*idsIt).first.subdetId()==EcalBarrel ) { 

	thePhi = EBDetId(  (*idsIt).first.rawId()  ).iphi() ; 
	if (phiFirstEB<-9998) phiFirstEB = thePhi;
	
	if (fabs( thePhi - phiFirstEB ) > 180 ) // this means cluster is crossing the iphi==1 <--> iphi==360 boundary
	  {
	    //	    std::cout << "====== EB thePhi was: " << thePhi << " and now has become: " << ( thePhi + 360 * signum( thePhi - phiFirstEB ) * -1 )<< std::endl;
	    thePhi -= 360 * signum( thePhi - phiFirstEB );
	  }
      }

      else if ( (!isBarrel) &&  (*idsIt).first.subdetId()==EcalEndcap ) {
	double ix = EEDetId(  (*idsIt).first.rawId()  ).ix() -50;
	double iy = EEDetId(  (*idsIt).first.rawId()  ).iy() -50;
	thePhi = atan2( ix,  iy );

	if (phiFirstEE<-9998) phiFirstEE = thePhi;

	if (fabs( thePhi - phiFirstEE ) > PI ) // this means cluster is crossing the phi==-pi <--> phi==pi boundary
	  {
	    //	    std::cout << "====== EE thePhi was: " << thePhi << " and now has become: " << ( thePhi + TWOPI * signum( thePhi - phiFirstEE ) * -1 )<< std::endl;
	    thePhi -= TWOPI * signum( thePhi - phiFirstEE );
	  }
	
      }
      
      if (phiLow==0 && phiHigh==0)
	{ phiLow  = thePhi;
	  phiHigh = thePhi;	 }
      if       (thePhi < phiLow)  phiLow  = thePhi;
      else if  (thePhi > phiHigh) phiHigh = thePhi;

    }

  float thePhiSize = phiHigh-phiLow;

  // PhiSize is in unit of crystals for EB 
  //                       radiants for EE (change this to crystals too)             
  return thePhiSize;
}




// returns distance between two crystals, in units of crystal size (both  for EB and EE)
float getPhiDistance(DetId theSeed, DetId theCry ){

  if ( theSeed.subdetId()==EcalBarrel && theCry.subdetId()==EcalBarrel)     // the barrel case
    {
      float theSeedPhi = EBDetId ( theSeed.rawId() ).iphi();
      float theCryPhi  = EBDetId ( theCry.rawId() ).iphi();
      if (fabs( theCryPhi - theSeedPhi ) > 180 )                            // this means cluster is crossing the iphi==1 <--> iphi==360 boundary
	{
	  theCryPhi -= 360 * signum( theCryPhi - theSeedPhi );
	}
      return (theCryPhi-theSeedPhi);    // this is in units of EB crystals
    }
  else if ( theSeed.subdetId()==EcalEndcap && theCry.subdetId()==EcalEndcap) // the endcap case
    {
      float theSeedX = EEDetId ( theSeed.rawId() ).ix() -50;
      float theSeedY = EEDetId ( theSeed.rawId() ).iy() -50;
      float theSeedPhi = atan2( theSeedX , theSeedY); 
      float theCryX  = EEDetId ( theCry.rawId() ).ix()  -50;
      float theCryY  = EEDetId ( theCry.rawId() ).iy()  -50;
      float theCryPhi  = atan2( theCryX , theCryY); 
      
      if ( fabs( theCryPhi - theSeedPhi ) > PI ) // this means cluster is crossing the phi==-pi <--> phi==pi boundary
	{
	  theCryPhi -= TWOPI * signum ( theCryPhi - theSeedPhi );
	}
      //        angular span              radius in units of EE crystal size 
      return ( (theCryPhi - theSeedPhi) * pow( theSeedX*theSeedX + theSeedY*theSeedY, 0.5) ) ;
    }
  else // unsupported mixed case
    {
      return -9999;
    }
  
}





// returns distance between two crystals, in units of crystal size (both  for EB and EE)
float getEtaDistance(DetId theSeed, DetId theCry ){

  if ( theSeed.subdetId()==EcalBarrel && theCry.subdetId()==EcalBarrel)     // the barrel case
    {
      float theSeedEta = EBDetId ( theSeed.rawId() ).ieta();
      float theCryEta  = EBDetId ( theCry.rawId() ).ieta();
      
      if ( theSeedEta * theCryEta > 0 )     // crystal and seed are on the same side of EB
	{ 
	  return (theCryEta - theSeedEta);
	}
      else	                           // crystal and seed are on opposite sides of EB
	{
	  if (theSeedEta<0) return (theCryEta - theSeedEta -1);
	  else              return (theCryEta - theSeedEta +1);
	}
    }
  else if ( theSeed.subdetId()==EcalEndcap && theCry.subdetId()==EcalEndcap) // the endcap case
    {
      // to be implemented w/ geometry...
      return -9999;
    }

  return -9999;

}
