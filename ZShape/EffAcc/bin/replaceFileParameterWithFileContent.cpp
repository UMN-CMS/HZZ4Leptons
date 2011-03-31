#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <fstream>
#include <iomanip>

int main( int argc, char **argv ) {

  std::cout << "\nThis program converts templates to proper values." << std::endl;


  if(argc<5 ){
    std::cout << "\nlacking arguments. Usage: inFile outFile Parameter ParamaterFileInfo  \n" << std::endl;
    return -1;
  }
  else if (argc >5){
    std::cout << "too many arguments.   Usage: efficiencyFileHisto SigFile BkgFile SigTemplate BkgTemplate SigOut Bkg Out \n" << std::endl;
    return -1;
  }

  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  std::string parameter = argv[3];
  std::string parFileName = argv[4];

  std::ifstream ifn(inputFileName.c_str());
  std::ofstream theOutFile;
  theOutFile.open (outputFileName.c_str(),std::ios::out);
  
  char str [50000];
  while ( ifn.getline(str,50000)){
    std::string mystr = str; //me likes strings, easy to use predefined algos to manipulate!
    std::string::size_type pos = mystr.find(parameter);
    if (pos != std::string::npos) 
      {
	//I really should first output the mystr.substr(0,pos) and at the end the mystr.substr(pos+parameter.size()) at the end
	if (pos > 0 ) theOutFile << mystr.substr(0,pos);
	//Open the par file and output all lines
	std::ifstream pfn(parFileName.c_str());
	char pstr [50000];
	while ( pfn.getline(pstr,50000)){
	  theOutFile << pstr << std::endl;
	  
	}
	pfn.close();
	if (mystr.size() > pos+parameter.size()) theOutFile << mystr.substr(pos+parameter.size()) << std::endl;

      }
    else theOutFile << mystr << std::endl;


    
  }
  
  
  ifn.close();


}

