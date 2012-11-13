#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <functional>

#include "ECALTime/EcalTimePi0/interface/EcalTimePi0TreeContent.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

// initial authors P. Govoni et al
// authors: S. Cooper and G. Franzoni (UMN)

//! main program
int main (int argc, char** argv)
{

  // default output file
  std::string outputRootName = "outputHistos.root";
  std::string stringGenericOption    = "--";
  std::string stringHelp             = "--help";
  std::string stringInputFileName    = "--i";
  std::string stringOutFileName      = "--o";
  std::string stringNumEvents        = "--n";

  std::vector<std::string> listOfFiles;
  int numEvents=-1;

  //gf: support development
  //std::cout << "\nargc:       " << argc << std::endl;
  //for (int v=0; v<argc; v++ ){      std::cout << "argument: " << v << " argv: " << argv[v] << std::endl;    }

  
  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, e.g. --help\n" << std::endl ;
    exit (1) ;  
  }

  // loop over input options
  for (int v=1; v<argc; v++ )
    {
      //std::cout << "argv number " << v << " is: " << argv[v] << std::endl;
      
      if (argv[v] == stringHelp) { // help message
	std::cout << " --help : display help" << std::endl ;
	std::cout << " -o : set name of output root file name (e.g. histograms.root)" << std::endl ;
	std::cout << " name of input file : list name of input files ntuples" << std::endl ;     
	exit(1);      }

      
      if (argv[v] == stringNumEvents) { // set number of events
	std::cout << "argument " << v << " is " << argv[v]  << std::endl ;
	numEvents=atoi(argv[v+1]);
	v++;
      }

      else if (argv[v] == stringOutFileName) { // set output file
	outputRootName = argv[v+1];
	v++;
      }

      // handle here the case of multiple arguments for input files
      else if (argv[v] == stringInputFileName){// && v<(argc-1) ) {

	for (int u=v+1; u<argc; u++) {
	  
	  if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
	    if ( 0==listOfFiles.size())  {std::cout << "no input files listed" << std::cout;}
	    //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
	    break;
	  }

	  else {  listOfFiles.push_back(argv[u]);
	    v++;
	  }

	}// loop on arguments following --i

	continue;

      }//end 'if input files'

      
      else
	{std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program


  
  if (listOfFiles.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }
  




  // Tree construction
  TChain * chain = new TChain ("EcalTimePi0Analysis") ;
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
  int nEntries = chain->GetEntries () ;
  if (numEvents==-1) numEvents = nEntries;
  std::cout << "\n\tFOUND " << nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: " <<  numEvents << " events" << std::endl;
  std::cout << "\tOutput file: " <<  outputRootName << std::endl;

           
  EcalTimePi0TreeContent treeVars ; 
  setBranchAddresses (chain, treeVars) ;



  //loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {
      chain->GetEntry (entry) ;

      bool speak=false;
      if (entry<10 || entry%100==0  || true) speak=true;

      if (speak) std::cout << "------> reading entry " << entry << " <------\n" ; 


      // loop on calorimetric quantities

      if (speak)  std::cout << "  found " << treeVars.nSuperClusters << " superclusters" << std::endl ;
      if (speak)  std::cout << "  found " << treeVars.nClusters << " basic clusters" << std::endl ;


        /////////////////////////////////////////////////////
	//loop on basic cluster

	for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	  {

	    float eBC=0; // calculate energy of BC for validation
	    for (int cryInBC=0; cryInBC < treeVars.nXtalsInCluster[bCluster]; cryInBC++){
	      eBC+= treeVars.xtalInBCEnergy[bCluster][cryInBC];}


	    if (speak)  std::cout << "\tbCluster: num"               << bCluster 
				  << "\t eBC: "                      << treeVars.clusterEnergy[bCluster]
				  << "\t eBC_predicted: "            << eBC
				  << "\n\t et: "                     << treeVars.clusterTransverseEnergy[bCluster]
				  << "\t predicted et: "             << treeVars.clusterEnergy[bCluster]*sin(2*atan(exp(-1* treeVars.clusterEta[bCluster] )) )
				  << " eta: "                        << treeVars.clusterEta[bCluster]
				  << "\n\t num crystals: "           << treeVars.nXtalsInCluster[bCluster]
				  << "\n\t\tfirst crystal:  \tieta " << treeVars.xtalInBCIEta[bCluster][0] 
				  << "\teta "                        << treeVars.xtalInBCEta[bCluster][0] 
				  << " \t energy "                   << treeVars.xtalInBCEnergy[bCluster][0] 
				  << " \t ADC "                      << treeVars.xtalInBCAmplitudeADC[bCluster][0] 
				  << " \t time "                     << treeVars.xtalInBCTime[bCluster][0] 
				  << std::endl;
	  }






      /////////////////////////////////////////////////////
      //loop on superclusters
      for (int SCindex = 0 ; SCindex < treeVars.nSuperClusters ; ++SCindex)
//        {
//	  if (speak) std::cout << " SC n  = " << SCindex << std::endl;
//	  if (speak) std::cout << " --->  Sup X = " << treeVars.superClusterX[SCindex] << " Sup Y = " << treeVars.superClusterY[SCindex] << " Sup Z = " << treeVars.superClusterZ[SCindex] << "\n" ;
//	  std::cout << "    found " << treeVars.nClustersInSuperCluster[SCindex] 
//                    << " clusters in supercluster\n" ;    
//          //PG loop over clusters in supercluster
//          for (int BCindex = treeVars.clusterIndexInSuperCluster[SCindex] ; 
//               BCindex < treeVars.clusterIndexInSuperCluster[SCindex] + 
//		 treeVars.nClustersInSuperCluster[SCindex] ; 
//               ++BCindex)
//            { 
//                
//	      std::cout << "      found " << treeVars.nXtalsInCluster[SCindex] 
//                        << " crystals in cluster\n" ;    
//
//              //PG loop over crystals in cluster
//              for (int XTLindex = treeVars.xtalIndexInCluster[BCindex] ; 
//                   XTLindex < treeVars.xtalIndexInCluster[BCindex] + 
//		     treeVars.nXtalsInCluster[BCindex] ; 
//                   ++XTLindex)
//                {
//                  if (!EBDetId::validHashIndex (treeVars.xtalHashedIndex[XTLindex]))
//                    {
//		      std::cerr << "ERROR crystal " 
//                                << treeVars.xtalHashedIndex[XTLindex] 
//                                << " has invalid DetId" << std::endl ;
//                      continue ;
//                    }          
//                  EBDetId dummy = EBDetId::unhashIndex (treeVars.xtalHashedIndex[XTLindex]) ;   
//		  std::cout << "        found crystal " 
//                            << treeVars.xtalHashedIndex[XTLindex]
//                            << " at (" << dummy.ieta () << "," << dummy.iphi ()
//                            << ") with energy " 
//                            << treeVars.xtalEnergy[XTLindex] 
//                            << " GeV \n" ;       
//                } //PG loop over crystals in cluster
//
//            } // loop over clusters in supercluster
//
//	  std::cout << "    found " << treeVars.nXtalsInSuperCluster[SCindex] 
//                    << " crystals in supercluster\n" ;    
//          //PG loop over crystals in SUPERcluster
//          for (int XTLindex = treeVars.xtalIndexInSuperCluster[SCindex] ; 
//               XTLindex < treeVars.xtalIndexInSuperCluster[SCindex] + 
//		 treeVars.nXtalsInSuperCluster[SCindex] ; 
//               ++XTLindex)
//            {
//	      std::cout << "      found crystal " 
//                        << treeVars.xtalHashedIndex[XTLindex]
//                        << " with energy " 
//                        << treeVars.xtalEnergy[XTLindex] 
//                        << " GeV \n" ;    
//            
//            } //loop over crystals in SUPERcluster
//        } //loop on superclusters

      std::cout << "  found " << treeVars.nXtals << " crystals\n" ;    
      //PG loop over crystals
      for (int XTLindex = 0 ; XTLindex < treeVars.nXtals ; ++XTLindex)
        {
          EBDetId dummy = EBDetId::unhashIndex (treeVars.xtalHashedIndex[XTLindex]) ;   
        } //PG loop over crystals

      

    } //PG loop over entries


  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;
          
  saving.Close () ;

  delete chain ;
  return 0 ;

}
