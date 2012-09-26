//------------------------------------------------------------
//-
//-
//- FILE BY JASON HAUPT FEB-27th-2008 12:47pm
//-
//- Purpose, to add the EE and EB txt files, and possibly
//-   to the HF when available.
////////////////////////////////////////////////////////////
#include <fstream>
#include <iomanip>

#include <string.h>
#include <iostream>



int main( int argc, char **argv ) {
  using namespace std;
  std::string theEffTitle      = argv[1];
  std::string textFileNameOut_ = argv[2];
  std::string textFileName1_   = argv[3];
  std::string textFileName2_   = argv[4];

  char str [1000]; char str1[1000];
  
  float binMin[1000];     float binMax[1000];  
  float values[1000];     float systErrP[1000];     float systErrM[1000];
  int   denominator[1000];

  int binNumber1[1000];
  float binMin1[1000];     float binMax1[1000];  
  float values1[1000];     float systErrP1[1000];     float systErrM1[1000];
  int   denominator1[1000];

  int binNumber2[1000];
  float binMin2[1000];     float binMax2[1000];  
  float values2[1000];     float systErrP2[1000];     float systErrM2[1000];
  int   denominator2[1000];
  //std::string theEffTitle;

  //This is the combined output file
  ofstream os;
  os.open(textFileNameOut_.c_str(),std::ios::out);
  
 
  ifstream is1(textFileName1_.c_str());
  // aiming at the table with numbers, inside the text file
  while ( !(strstr(str,"Nbin") &&           strstr(str,"binMin") &&
            strstr(str,"binMax") &&         strstr(str,"value") &&
            strstr(str,"systP") &&          strstr(str,"systM") &&
            strstr(str,"denom"))
          )
    {
      is1.getline(str,1000);
    }

  *str=*str1;
  ifstream is2(textFileName2_.c_str());
  // aiming at the table with numbers, inside the text file
  while ( !(strstr(str,"Nbin") &&           strstr(str,"binMin") &&
            strstr(str,"binMax") &&         strstr(str,"value") &&
            strstr(str,"systP") &&          strstr(str,"systM") &&
            strstr(str,"denom"))
          )
    {is2.getline(str,1000);}

  //Now I will output the header for the new file.
  os << "#put your comments after a '#'.\n";
  os << "\nefficiency name: " << theEffTitle <<"\n";
  os << "dimension: 1 \n\n\n";
  os << "Nbin"
               << "\t" << "binMin"
               << "\t"  << "binMax"
               << "\t"  <<  "value"
               << "\t"  << "systP"
               << "\t"  << "systM"
               << "\t" << "denom"
               << "\n";

  

  int binCounter1=0;
  while (       is1.getline(str,1000)){
    if ( strstr(str,"#") ) { continue;}
    //std::cout << str << std::endl; // just for debug
    sscanf (str,"%d %f %f %f %f %f %d", &binNumber1[binCounter1], &binMin1[binCounter1], &binMax1[binCounter1],
            &values1[binCounter1], &systErrP1[binCounter1], &systErrM1[binCounter1], &denominator1[binCounter1]);
    binCounter1++;
  }
  is1.close();

  int binCounter2=0;
  while (       is2.getline(str,1000)){
    if ( strstr(str,"#") ) { continue;}
    //std::cout << str << std::endl; // just for debug
    sscanf (str,"%d %f %f %f %f %f %d", &binNumber2[binCounter2], &binMin2[binCounter2], &binMax2[binCounter2],
            &values2[binCounter2], &systErrP2[binCounter2], &systErrM2[binCounter2], &denominator2[binCounter2]);
    binCounter2++;
  }
  is2.close();

  if (binCounter1 != binCounter2){std::cout << " Oh no, differet number of bins, I don't like this "<< binCounter1 << " " << binCounter2 << std::endl; return 2;}

  for (int r=0; r<binCounter1; ++r)
    {
      if (((values1[r] > 0) && (values2[r] > 0.)) || ((denominator1[r] > 0) && (denominator2[r] > 0)))
	{
	  binMin[r]=binMin1[r];
          binMax[r]=binMax1[r];
	  values[r]=0.;
	  systErrP[r]=1.;
	  systErrM[r]=-1.;
	  denominator[r]=0;
	}
      else if  (((values1[r] > 0) ||(denominator1[r] > 0) ))
	{
	  binMin[r]=binMin1[r];
          binMax[r]=binMax1[r];
	  values[r]=values1[r];
	  systErrP[r]=systErrP1[r];
	  systErrM[r]=systErrM1[r];
	  denominator[r]=denominator1[r];
	}
      else if  (((values2[r] > 0) ||(denominator2[r] > 0)))
	{
	  binMin[r]=binMin2[r];
          binMax[r]=binMax2[r];
	  values[r]=values2[r];
	  systErrP[r]=systErrP2[r];
	  systErrM[r]=systErrM2[r];
	  denominator[r]=denominator2[r];
	}
      else 
	{
	  binMin[r]=binMin1[r];
          binMax[r]=binMax1[r];
	  values[r]=0.;
	  systErrP[r]=1.;
	  systErrM[r]=-1.;
	  denominator[r]=0;
	}

      os           << (r+1)
		   << "\t" << std::setprecision(4) << binMin[r]
		   << "\t" << std::setprecision(4) << binMax[r]
		   << "\t" << std::setprecision(4) << values[r]
		   << "\t" << std::setprecision(4) << systErrP[r]
		   << "\t" << std::setprecision(4) << systErrM[r]
		   << "\t" << denominator[r]
		   << "\n";
    }
  os.close();

}




// End of file
