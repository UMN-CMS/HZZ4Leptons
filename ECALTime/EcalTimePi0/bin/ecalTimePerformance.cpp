#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <math.h> 
#include <assert.h>

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"
#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"
#include "ECALTime/EcalTimePi0/interface/EcalObjectTime.h"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"

typedef std::set<std::pair<int,int> > SetOfIntPairs;

// authors: S. Cooper and G. Franzoni (UMN)

#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins     35
#define numAoSigmaBins  25

#define NSlices 54   // # of slices in log(A)
#define LogStep 0.05 // size of Slices log(A)



// -------- Globals ----------------------------------------
EcalTimeTreeContent treeVars_; 
TFile* saving_;
std::vector<std::string> listOfFiles_;
bool speak_=false;
char buffer_ [75];
std::string bufferTitle_; 
// default settings
std::string outputRootName_ = "outputHistos.root";
int   numEvents_      = -1;
unsigned int  minRun_ = 0;
unsigned int  maxRun_ = 9999999;
unsigned int  minLS_ = 0;
unsigned int  maxLS_ = 9999999;
float eTGammaMinEB_   = 0.2;
float s4s9GammaMinEB_ = 0.85;
float eTPi0MinEB_     = 0.65;
float eTGammaMinEE_   = 0.250;
float s4s9GammaMinEE_ = 0.85;
float eTPi0MinEE_     = 0.800;
float swissCrossMaxEB_ = 0.95; // 1-E4/E1
float swissCrossMaxEE_ = 0.95; // 1-E4/E1
// based on range and bins, the bin width is 50 ps
float rangeTDistro_ = 3; // 1-E4/E1
int   binsTDistro_  = 120; // 1-E4/E1
std::vector<std::vector<double> > trigIncludeVector;
std::vector<std::vector<double> > trigExcludeVector;
std::vector<std::vector<double> > ttrigIncludeVector;
std::vector<std::vector<double> > ttrigExcludeVector;



int  minEntriesForFit_ = 7;
int  flagOneVertex_ = 0;
bool limitFit_(true); 
//std::string fitOption_(""); // default: use chi2 method
std::string fitOption_("L"); // use likelihood method


// -------- Histograms -------------------------------------
TH1 * nVertices_;
TH1F* mass_;
TH1F* dZvertices_;
TH1F* Zvertices_;


// ---------------------------------------------------------------------------------------
// - Function to decide to include/exclude event based on the vectors passed for triggers 
bool includeEvent(int* triggers,
    int numTriggers,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
 bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=includeVector.size();++i){
      if(includeVector[i].size()==1 && triggers[ti]==includeVector[i][0]) keepEvent=true;
      else if(includeVector[i].size()==2 && (triggers[ti]>=includeVector[i][0] && triggers[ti]<=includeVector[i][1])) keepEvent=true;
    }
  }
  if(!keepEvent)
    return false;

  keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=excludeVector.size();++i){
      if(excludeVector[i].size()==1 && triggers[ti]==excludeVector[i][0]) keepEvent=false;
      else if(excludeVector[i].size()==2 && (triggers[ti]>=excludeVector[i][0] && triggers[ti]<=excludeVector[i][1])) keepEvent=false;
    }
  }

  return keepEvent;
}

// ---------------------------------------------------------------------------------------
// ------- Function to decide to include/exclude event based on the vectors passed -------
bool includeEvent(double eventParameter,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
  bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for(uint i=0; i!=includeVector.size();++i){
    if(includeVector[i].size()==1 && eventParameter==includeVector[i][0])
      keepEvent=true;
    else if(includeVector[i].size()==2 && (eventParameter>=includeVector[i][0] && eventParameter<=includeVector[i][1]))
      keepEvent=true;
  }
  if(!keepEvent) // if it's not in our include list, skip it
    return false;

  keepEvent = true;
  for(uint i=0; i!=excludeVector.size();++i){
    if(excludeVector[i].size()==1 && eventParameter==excludeVector[i][0])
      keepEvent=false;
    else if(excludeVector[i].size()==2 && (eventParameter>=excludeVector[i][0] && eventParameter<=excludeVector[i][1]))
      keepEvent=false;
  }

  return keepEvent; // if someone includes and excludes, exseedion will overrule

}



// ---------------------------------------------------------------------------------------
// ------------------ Function to split arg input lists by comma -------------------------
std::vector<std::string> split(std::string msg, std::string separator)
{
  boost::char_separator<char> sep(separator.c_str());
  boost::tokenizer<boost::char_separator<char> > tok(msg, sep );
  std::vector<std::string> token;
  for ( boost::tokenizer<boost::char_separator<char> >::const_iterator i = tok.begin(); i != tok.end(); ++i ) {
    token.push_back(std::string(*i));
  }
  return token;
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to generate include/exclude vectors -----------------------
void genIncludeExcludeVectors(std::string optionString,
			      std::vector<std::vector<double> >& includeVector,
			      std::vector<std::vector<double> >& excludeVector)
{
  std::vector<std::string> rangeStringVector;
  std::vector<double> rangeIntVector;

  if(optionString != "-1"){
    std::vector<std::string> stringVector = split(optionString,",") ;

    for (uint i=0 ; i<stringVector.size() ; i++) {
      bool exclude = false;

      if(stringVector[i].at(0)=='x'){
        exclude = true;
        stringVector[i].erase(0,1);
      }
      rangeStringVector = split(stringVector[i],"-") ;

      rangeIntVector.clear();
      for(uint j=0; j<rangeStringVector.size();j++) {
        rangeIntVector.push_back(atof(rangeStringVector[j].c_str()));
      }
      if(exclude) excludeVector.push_back(rangeIntVector);
      else includeVector.push_back(rangeIntVector);

    }
  }
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to parse the command-line arguments------------------------
void parseArguments(int argc, char** argv)
{
  std::string stringGenericOption    = "--";
  std::string stringHelp             = "--help";
  std::string stringInputFileName    = "--i";
  std::string stringOutFileName      = "--o";
  std::string stringETGammaMinEB     = "--eTGammaMinEB";
  std::string strings4s9GammaMinEB   = "--s4s9GammaMinEB";
  std::string stringeTPi0MinEB       = "--eTPi0MinEB";
  std::string stringETGammaMinEE     = "--eTGammaMinEE";
  std::string strings4s9GammaMinEE   = "--s4s9GammaMinEE";
  std::string stringeTPi0MinEE       = "--eTPi0MinEE";
  std::string stringminAmpliOverSigma= "--minAOverSigma";
  std::string stringNumEvents        = "--n";
  std::string stringMinRun           = "--minRun";
  std::string stringMaxRun           = "--maxRun";
  std::string stringMinLS            = "--minLS";
  std::string stringMaxLS            = "--maxLS";
  std::string vertex                 = "--vertex";
  std::string stringTriggers         = "--trig";
  std::string stringTechTriggers     = "--techTrig";

  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, or use --help\n" << std::endl ;
    exit (1) ;  
  }

  // loop over input options
  for (int v=1; v<argc; v++ )
  {
    //std::cout << "argv number " << v << " is: " << argv[v] << std::endl;

    if (argv[v] == stringHelp) { // help message
      std::cout << " --help : display help" << std::endl ;
      std::cout << " --o : set name of output root file name (e.g. histograms.root)" << std::endl ;
      std::cout << " --n : number of events" << std::endl ;
      std::cout << " --i <list of strings> list of input files" << std::endl ;     
      std::cout << " --eTGammaMinEB: min eT for EB gammas" << std::endl;
      std::cout << " --s4s9GammaMinEB: min EB shower shape" << std::endl;
      std::cout << " --eTPi0MinEB min eT for EB pi0 candidate" << std::endl;
      std::cout << " --eTGammaMinEE: min eT for EE gammas" << std::endl;
      std::cout << " --s4s9GammaMinEE: min EE shower shape" << std::endl;
      std::cout << " --eTPi0MinEE: min eT for EE pi0 candidate" << std::endl;
      std::cout << " --minAOverSigma: min ampli considered for time" << std::endl;
      std::cout << " --minRun: lowest run number considered" << std::endl;
      std::cout << " --maxRun: highest run number considered" << std::endl;
      std::cout << " --minLS: lowest lumi section number considered" << std::endl;
      std::cout << " --maxLS: highest lumi section number considered" << std::endl;
      std::cout << " --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)" << std::endl;
      std::cout << " --trig: L1 triggers to include (exclude with x)" << std::endl;
      std::cout << " --techTrig: L1 technical triggers to include (exclude with x)" << std::endl;
      exit(1);      }


    else if (argv[v] == stringNumEvents) { // set number of events
      std::cout << "events number" << std::endl;
      numEvents_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxLS) { // set last LS of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinLS) { // set first LS of interval to be considered 
      std::cout << "min LS number" << std::endl;
      minLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxRun) { // set last run of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinRun) { // set first run of interval to be considered 
      std::cout << "min run number" << std::endl;
      minRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEB) { // choose et cut for EB single cluster
      eTGammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEB) { // choose cut for EB shower shape
      s4s9GammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEB) { // choose et cut for EB pi0 candidate
      eTPi0MinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEE) { // choose et cut for EE single cluster
      eTGammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEE) { // choose cut for EE shower shape
      s4s9GammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEE) { // choose et cut for EE pi0 candidate
      eTPi0MinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringOutFileName) { // set output file
      outputRootName_ = argv[v+1];
      v++;
    }
    else if (argv[v] == stringminAmpliOverSigma) { // set min amplitude considered for time measurement
      minAmpliOverSigma_  = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == vertex) { // collect requirement for one vertex only or not
      flagOneVertex_  = atof(argv[v+1]);
       if (flagOneVertex_!=0 && flagOneVertex_!=1 && flagOneVertex_!=2){
         std::cout << "Not a valid value for flagOneVertex_ (0,1,2). Returning." << std::endl;
	 exit (1);}
       v++;
    } 
    else if (argv[v] == stringTriggers) { // set L1 triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),trigIncludeVector,trigExcludeVector);
      v++;
    }
    else if (argv[v] == stringTechTriggers) { // set L1 technical triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),ttrigIncludeVector,ttrigExcludeVector);
      v++;
    }
    // handle here the case of multiple arguments for input files
    else if (argv[v] == stringInputFileName){// && v<(argc-1) ) 

      for (int u=v+1; u<argc; u++) {

        if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
          if ( 0==listOfFiles_.size())  {std::cout << "no input files listed" << std::cout;}
          //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
          break;
        }

        else {  listOfFiles_.push_back(argv[u]);
          v++;
        }

      }// loop on arguments following --i

      continue;

    }//end 'if input files'

    else
    {std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program
}

// ---------------------------------------------------------------------------------------
// -------------------- struct holding sed of  histograms  -------------------------------
struct HistSet{
  
  //book histogram set w/ common suffix inside the provided TFileDirectory
  //void book(edm::Service<TFileService>& td,const std::string&);
  void book(TFileDirectory subDir,const std::string&);
  
  // fill all histos of the set with the two electron candidates
  void fill(int sc1, int sc2, int cl1, int cl2);
  
  TH1 * nVertices_;
  TH1F* seedTime_;
  TH1F* secondTime_;
  TH1F* clusterTime_;

} theHists;


void HistSet::book(TFileDirectory subDir, const std::string& post) {

  nVertices_=subDir.make<TH1F>("num vertices","num vertices; num vertices",41,-0.5,40.5);
  mass_         =(TH1F*) subDir.make<TH1F>("mass","mass; m(ele,ele) [GeV]",80,50,130);

  seedTime_            =(TH1F*) subDir.make<TH1F>("seed time","seed time; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  secondTime_          =(TH1F*) subDir.make<TH1F>("second time","second time; t_{second} [ns]; num. secs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  clusterTime_         =(TH1F*) subDir.make<TH1F>("cluster time","cluster time; t_{cluster} [ns]; num. clusters/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);


}
  
void HistSet::fill(int sc1, int sc2, int bc1, int bc2 ){

  float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
  math::PtEtaPhiELorentzVectorD  el1(et1  ,
				     treeVars_.superClusterEta[sc1],
				     treeVars_.superClusterPhi[sc1],
				     treeVars_.superClusterRawEnergy[sc1] );  
  float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
  math::PtEtaPhiELorentzVectorD  el2(et2 ,
				     treeVars_.superClusterEta[sc2],
				     treeVars_.superClusterPhi[sc2],
				     treeVars_.superClusterRawEnergy[sc2] );
  math::PtEtaPhiELorentzVectorD diEle = el1;
  diEle += el2;
  
  // ////////////////////////
  mass_      ->Fill(diEle.M());
  float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
  dvertex       = sqrt(dvertex);

  ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
  ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);


  int vtxOfThisEle=-99;
  // look for the vertex which electrons  are attached to: 
  for(int u=0; u<treeVars_.nVertices; u++){
    // matching done with 1mm tolerance
    if( fabs(treeVars_.superClusterVertexZ[sc2]-treeVars_.vtxZ[u]) < 0.1) {       vtxOfThisEle=u;     }
    //std::cout << u << "\t" << treeVars_.superClusterVertexZ[sc2] << "\t" << treeVars_.vtxZ[u] << std::endl;
  }

  // take care of the seeds
  seedTime_            -> Fill(bcTime1.seedtime);  seedTime_->Fill(bcTime2.seedtime); 
  clusterTime_         -> Fill(bcTime1.time);              clusterTime_ ->Fill(bcTime2.time);

}
// end HistSet::fill


// ---------------------------------------------------------------------------------------
// ------------------ Function to initialize the histograms ------------------------------
void initializeHists(TFileDirectory subDir){

//  mass_         = subDir.make<TH1F>("mass global","mass (global); m(ele,ele) [GeV]",80,50,130);
//  dZvertices_   = subDir.make<TH1F>("dZvertices global","dZvertices (global); #DeltaZ(ele_{1},ele_{2}) [cm]",250,0,25);
//  Zvertices_    = subDir.make<TH1F>("Zvertices global","Zvertices (global); z vertex [cm]",250,-25,25);
  nVertices_=subDir.make<TH1F>("num vertices global","num vertices (global); num vertices",41,-0.5,40.5);

}//end initializeHists




// ---------------------------------------------------------------------------------------
//! main program
int main (int argc, char** argv)
{
  // First parse arguments
  parseArguments(argc, argv);

  if (listOfFiles_.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles_.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }

  // Tree construction
  // FIX should turn this string into a configurable 
  TChain * chain = new TChain ("EcalTimeAnalysis") ;  // ntuple producer in CMSSW CVS
  //TChain * chain = new TChain ("EcalTimePi0Analysis") ;  // ntuple producer in UserCode/UMN space
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
  int nEntries = chain->GetEntries () ;
  if (numEvents_==-1) numEvents_ = nEntries;
  std::cout << "\n\tFOUND "         <<  listOfFiles_.size() << " input files" << std::endl ;    
  std::cout << "\n\tFOUND "         <<  nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: "    <<  numEvents_ << " events" << std::endl;
  std::cout << "\tOutput file: "    <<  outputRootName_ << std::endl;
  std::cout << "\tminAOverSigma: "  <<  minAmpliOverSigma_ << std::endl;
  std::cout << "\teTGammaMinEB: "   <<  eTGammaMinEB_ << std::endl;
  std::cout << "\ts4s9GammaMinEB: " <<  s4s9GammaMinEB_ << std::endl;
  std::cout << "\teTPi0MinEB: "     <<  eTPi0MinEB_ << std::endl;
  std::cout << "\teTGammaMinEE: "   <<  eTGammaMinEE_ << std::endl;
  std::cout << "\ts4s9GammaMinEE: " <<  s4s9GammaMinEE_ << std::endl;
  std::cout << "\teTPi0MinEE: "     <<  eTPi0MinEE_ << std::endl;
  std::cout << "\tminRun: "         <<  minRun_ << std::endl;
  std::cout << "\tmaxRun: "         <<  maxRun_ << std::endl;
  std::cout << "\tminLS: "          <<  minLS_ << std::endl;
  std::cout << "\tmaxLS: "          <<  maxLS_ << std::endl;
	
  setBranchAddresses (chain, treeVars_);
  
  // setting up the TFileService in the ServiceRegistry;
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());
  std::vector<edm::ParameterSet> psets;
  edm::ParameterSet pSet;
  pSet.addParameter("@service_type",std::string("TFileService"));
  pSet.addParameter("fileName",std::string("TimePerf-plots.root")); // this is the file TFileService will write into
  psets.push_back(pSet);
  static edm::ServiceToken services(edm::ServiceRegistry::createSet(psets));
  static edm::ServiceRegistry::Operate operate(services);
  edm::Service<TFileService> fs;

  TFileDirectory subDirEBEB=fs->mkdir("EBEB");  
  HistSet plotsEBEB;
  plotsEBEB.book(subDirEBEB,std::string("EBEB"));
  
  TFileDirectory subDirEEEE=fs->mkdir("EEEE");  
  HistSet plotsEEEE;
  plotsEEEE.book(subDirEEEE,std::string("EEEE"));
  
  timeCorrector theCorr;
  std::cout << "\ncreated object theCorr to be used for timeVsAmpliCorrections" << std::endl;
  std::cout << "\ninitializing theCorr" << std::endl;
  theCorr.initEB( "EB" );
  theCorr.initEE( "EE" );


  //Initialize output root file
  //saving_ = new TFile(outputRootName_.c_str (),"recreate");

  // Initialize the histograms
  TFileDirectory subDirGeneral=fs->mkdir("General");  
  initializeHists(subDirGeneral);

  int eventCounter = 0;
  /////////////////////////////////////////////////////
  // Main loop over entries
  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
  {
    chain->GetEntry (entry) ;
    // Keep the event?
    bool keepEvent = includeEvent(treeVars_.l1ActiveTriggers,
        treeVars_.l1NActiveTriggers,trigIncludeVector,trigExcludeVector)
            && includeEvent(treeVars_.l1ActiveTechTriggers,
                treeVars_.l1NActiveTechTriggers,ttrigIncludeVector,ttrigExcludeVector);
    if(!keepEvent)
      continue;

    
    // do analysis if the run is in the desired range  
    if( treeVars_.runId<minRun_  || maxRun_<treeVars_.runId) continue;
    
    // do analysis if the LS is in the desired range  
    if( treeVars_.lumiSection<minLS_  || maxLS_<treeVars_.lumiSection) continue;
    
    bool verticesAreOnlyNextToNominalIP;
    int  count=0;
    
    for(int v=0; v<treeVars_.nVertices; v++  )
	{        if (fabs(treeVars_.vtxZ[0])<15) count++; }
    
    if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
    else                                                        verticesAreOnlyNextToNominalIP = false;
    
    //    --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
    if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
    if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) )  continue;
    
    // if evet being actually processed, increment counter of analyzed events
    eventCounter++;
    
    speak_=false;
    if (entry<10 || entry%10000==0) speak_=true;

    if (speak_)  std::cout << "\n\n------> reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ; 
    if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;


    ///////////////////////////////////////////////////////////////////////
    // outer loop on supercluster
    for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){

      float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
      if (et1<20) continue;

      math::PtEtaPhiELorentzVectorD  el1(et1  ,
					 treeVars_.superClusterEta[sc1],
					 treeVars_.superClusterPhi[sc1],
					 treeVars_.superClusterRawEnergy[sc1] );

      ///////////////////////////////////////////////////////////////////////
      // inner loop on supercluster
      for (int sc2=(sc1+1); sc2<treeVars_.nSuperClusters; sc2++){

	float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
	if (et2<20) continue;
	
	math::PtEtaPhiELorentzVectorD  el2(et2 ,
					   treeVars_.superClusterEta[sc2],
					   treeVars_.superClusterPhi[sc2],
					   treeVars_.superClusterRawEnergy[sc2] );

	
	math::PtEtaPhiELorentzVectorD diEle = el1;
	diEle += el2;

	// ////////////////////////
	mass_      ->Fill(diEle.M());//GF
	nVertices_->Fill(treeVars_.nVertices);

	float tmpEne=-9999;
	// loop on BC and match to sc1  ===============
	int bc1=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( (pow(treeVars_.superClusterEta[sc1]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc1]-treeVars_.clusterPhi[bc],2) ) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc1=bc;
	  }// end - if good bc candidate
	}// end - loop over BC

	
	tmpEne=-9999;
	// loop on BC and match to sc2 ==============
	int bc2=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( pow(treeVars_.superClusterEta[sc2]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc2]-treeVars_.clusterPhi[bc],2) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc2=bc;
	  }// end - if good bc candidate
	}// end - loop over BC
	
	// protect in case of no matching
	if(bc1==-1 || bc2==-1) continue;
	if(0) {
	std::cout << "\n\nsc1 : " << treeVars_.superClusterEta[sc1] << " " << treeVars_.superClusterPhi[sc1] << " " << treeVars_.superClusterRawEnergy[sc1] << std::endl;
	std::cout << "bc1 : " << treeVars_.clusterEta[bc1] << " " << treeVars_.clusterPhi[bc1] << " " << treeVars_.clusterEnergy[bc1] << "\n"<< std::endl;
	std::cout << "sc2 : " << treeVars_.superClusterEta[sc2] << " " << treeVars_.superClusterPhi[sc2] << " " << treeVars_.superClusterRawEnergy[sc2] << std::endl;
	std::cout << "bc2 : " << treeVars_.clusterEta[bc2] << " " << treeVars_.clusterPhi[bc2] << " " << treeVars_.clusterEnergy[bc2] << std::endl;
	}
	
	ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
	ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);

	if(! (bcTime1.isvalid && bcTime2.isvalid) ) continue;

	// fill the structures which hold all the plots
	if      ( fabs(treeVars_.clusterEta[bc1])<1.4    &&  fabs(treeVars_.clusterEta[bc2])<1.4 ){
 	  plotsEBEB.fill(sc1,sc2, bc1,bc2);

	  float energyRatio1 = treeVars_.xtalInBCEnergy[bc1][bcTime1.seed];
	  if(bcTime1.second>-1) {energyRatio1 /= treeVars_.xtalInBCEnergy[bc1][bcTime1.second]; }
	  else { energyRatio1 /= 99999; }
	  float energyRatio2 = treeVars_.xtalInBCEnergy[bc2][bcTime2.seed];
	  if(bcTime2.second>-1) {energyRatio2 /= treeVars_.xtalInBCEnergy[bc2][bcTime2.second]; }
	  else { energyRatio2 /= 99999; }

	  
	}// if EBEB, and subcases
	else if ( fabs(treeVars_.clusterEta[bc1])>1.5    &&  fabs(treeVars_.clusterEta[bc2])>1.5 ) 	  plotsEEEE.fill(sc1,sc2, bc1,bc2);

	// if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
	break;	
	
      }// end loop sc2
    }// end loop sc1
    
  }   // end of loop over entries
  

  delete chain ;
  
  return 0 ;
}
