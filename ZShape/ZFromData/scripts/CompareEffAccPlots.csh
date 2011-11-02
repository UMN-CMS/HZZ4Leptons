#!/bin/tcsh -f
########################################
# Input and outpur directory arguments
#
#
########################################

set outdir=$1
set filetype=$2
set deffile="/home/grad/haupt/CMSSW/ZShape/CMSSW_2_1_12/src/ZShape/EffAcc/data"

root -l -b <<!
#include <map.h>
typedef std::vector<std::string> zdef;
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(10);
.L tdrstyle.C
setTDRStyle();

zdef zdefvec;

std::cout << " I'm Trying here " << std::endl;
std::map<std::string, zdef> zdefvecvec;
//zdefvec.push_back("Tight-EB-Loose-ECAL");
zdefvec.push_back("Tight-ECAL-Loose-ECAL"); 
zdefvec.push_back("Tight-ECAL-HF");
//zdefvec.push_back("Golden-EB-EB");


gStyle->SetTitleAlign(22);

.L compareEffAcc.C++
//81.17
//8.0593665
//16.118733
//6.88835 for M20.... (cross fingers)...
EffAccHistos *myEffAccHistos = new EffAccHistos("/data/whybee0c/user/haupt/Electrons/TNPTREE/ZeeFromDataJeremyTEMPLATE/histoJeremy09_only1hlt_TRIAL_Z20.root","/data/whybee0c/user/haupt/Electrons/TNPTREE/testProdZEfficiencySmeared_10MTEMPLATE/histo_10M_part.root",zdefvec, 6.88835);
myEffAccHistos->getFiles();
myEffAccHistos->getHistos();
myEffAccHistos->printSumHistos("${filetype}");
myEffAccHistos->printIndividualHistos("${filetype}");
myEffAccHistos->MakeResultFiles("${deffile}");

!

##echo "now we are making the webpages"
if (! -d $outdir) then
mkdir ${outdir}
endif

if (! -d $outdir/Compare ) then
mkdir ${outdir}/Compare
endif

if (! -d $outdir/FinalCompare ) then
mkdir ${outdir}/FinalCompare
endif

cp Z* ${outdir}/.

mv ${outdir}/Z_Compare* ${outdir}/Compare/.
./makeEffAccWebPage.csh ${outdir} ${filetype}


echo "That's all folks"
rm Z*${filetype}

##End of Efficiency*Acceptance Script
