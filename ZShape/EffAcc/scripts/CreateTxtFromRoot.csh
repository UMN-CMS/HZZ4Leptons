#!/bin/tcsh -f
###########################################################
#
# Simple Script to make the txt files from the root file
###########################################################

set rootfile=$1

foreach Eff (Iso GsfClass Sc ElecId HLT )
    foreach Side ( EE EB )
	set physvar="Eta"
	if ( $Eff == "Iso") then
	    set physvar="Pt"
	endif
	efficiencyRoot2Txt $rootfile ${Eff}Eff${Side}${physvar}.txt ${Eff}Eff ${physvar}${Side}
    end
    if ( $Eff != "Iso") then
        set myEff=${Eff}
        echo "THE EFF "${Eff}
	echo "my eff "${myEff}
	if ( $Eff == "ElecId") then
	    set myEff="ElectronId"
        else if ( $Eff == "Sc") then
	    set myEff="Supercluster"
        else if ( $Eff == "GsfClass") then
	    set myEff="GsfTrack"
     	endif
	efficiencyTxtCombiner ${myEff}-Eta ${myEff}-Eta.txt ${Eff}EffEBEta.txt ${Eff}EffEEEta.txt
    endif
end
echo "Done with the conversions"


#End of FILE
