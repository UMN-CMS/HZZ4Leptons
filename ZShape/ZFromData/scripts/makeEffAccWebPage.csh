#!/bin/tcsh -f
#################################################
# This super simple script makes a basic webpage
#
#################################################
##
## Argument 1 is the directory that contains the plots 
##
set mydir=`pwd`
set PlotDir=$1
set plottype=$2

cat > ${PlotDir}/index.html <<EOF
<HTML>

<HEAD><TITLE>ZShape Result WebPages for `date`</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
<h2><A name="EB"><FONT color="Black">ZShape Result Web Pages for `date`</FONT></A><BR></h2>
<a href="Compare">Link to MC to Data comparision histograms</a>
<p></p>
<a href="FinalCompare">Link to FINAL MC to Data comparision histograms</a>
<h3><A name="EB"><FONT color="Blue">Result Histograms</FONT></A><BR></h3>

<A HREF=ZResult_With_MC_Z0_Y.png> <img height="300" src="ZResult_With_MC_Z0_Y.png"> </A>
<A HREF=ZResult_Z0_Y.png> <img height="300" src="ZResult_Z0_Y.png"> </A>
<A HREF=ZResult_With_MC_With_FD_Z0_Y.png> <img height="300" src="ZResult_With_MC_With_FD_Z0_Y.png"> </A>
<p>These plots represent the Efficiency X Acceptance corrected Z rapdity distribution. 
This is also compared to the GENERATOR rapidity distribution in RED and the MC distribution in BLUE. </p>
<h3><A name="EB"><FONT color="Blue">Individual ZFrom Data Histograms</FONT></A><BR></h3>

EOF

foreach FromData (`/bin/ls ${PlotDir} |grep "ZFromData_"| grep "png"  `)
   echo '<A HREF='${FromData}'> <img height="300" src="'${FromData}'"> </A>' >> ${PlotDir}/index.html
end

echo '<p>These are the individual reconstructed Z Rapidity distributions of the different Z definitions</p>' >> ${PlotDir}/index.html

echo '<h3><A name="EB"><FONT color="Blue">Combined Data Z Rapidity</FONT></A><BR></h3>' >> ${PlotDir}/index.html
echo '<A HREF=ZEffRapTotal_Z0_Y.png> <img height="300" src="ZEffRapTotal_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<A HREF=ZRapTotalStacked_Z0_Y.png> <img height="300" src="ZRapTotalStacked_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<A HREF=ZRapTotalEach_Z0_Y.png> <img height="300" src="ZRapTotalEach_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<p>These are the combined reconstructed Z Rapidity distributions. This is shown in a combined form, a stacked form, and individual Z definition form.</p>' >> ${PlotDir}/index.html

echo '<h3><A name="EB"><FONT color="Blue">Efficiency X Acceptance For Single Z Definitions</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach EffAcc (`/bin/ls ${PlotDir} |grep "ZMC_EffAcc_" | grep "png" `)
   echo '<A HREF='${EffAcc}'> <img height="300" src="'${EffAcc}'"> </A>' >> ${PlotDir}/index.html
end

echo '<p>These are the Efficiency X Acceptance MC Z rapidity distributions for different Z definitions</p>' >> ${PlotDir}/index.html

echo '<h3><A name="EB"><FONT color="Blue">Efficiency X Acceptance </FONT></A><BR></h3>' >> ${PlotDir}/index.html
echo '<A HREF=ZEffAccTotal_Z0_Y.png> <img height="300" src="ZEffAccTotal_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<A HREF=ZEffAccTotalStacked_Z0_Y.png> <img height="300" src="ZEffAccTotalStacked_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<A HREF=ZEffAccTotalEach_Z0_Y.png> <img height="300" src="ZEffAccTotalEach_Z0_Y.png"> </A>' >> ${PlotDir}/index.html
echo '<p>These are the combined Efficiency X Acceptance MC Z Rapidity distributions. This is shown in a combined form, a stacked form, and individual Z definition form.</p>' >> ${PlotDir}/index.html

echo '<h3><A name="EB"><FONT color="Blue">Full MC Rapdity Distribution</FONT></A><BR></h3>' >> ${PlotDir}/index.html
echo '<A HREF=ZMCRAPFull_Z0_Y.png> <img height="300" src="ZMCRAPFull_Z0_Y.png"> </A>' >> ${PlotDir}/index.html



#NOW I MOVE ONTO THE COMPARISON SCRIPTS
set oldPlotDir=${PlotDir}
set PlotDir=${PlotDir}/Compare
cat > ${PlotDir}/index.html <<EOF
<HTML>

<HEAD><TITLE>ZShape MC to Data Comparison WebPages for `date`</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
<h2><A name="EB"><FONT color="Black">ZShape Web Pages for `date`</FONT></A><BR></h2>
<h3><A name="EB"><FONT color="Blue">Comparison Histograms</FONT></A><BR></h3>

<p>These plots represent the comparision to fully simulated and reconstructed MC to MC truth distributions. 
The comparision is only truely valid for cuts after the data selection (generally the GSF step).
They are also listed in the reverse order of the cuts. </p>

<p>The different ZDefs available are... </p>
EOF
foreach Ztype ( Tight-ECAL-Loose-ECAL Tight-ECAL-HF )
   if (`/bin/ls ${PlotDir} |grep "png" |grep -c $Ztype`) then
      echo '<A href="#'${Ztype}'"><FONT color="Black">Comparisons for Z Def '${Ztype}'</FONT></A><BR>' >> ${PlotDir}/index.html
   endif
end

foreach Ztype ( Tight-ECAL-Loose-ECAL Tight-ECAL-HF )
   if (`/bin/ls ${PlotDir} |grep "png" |grep -c $Ztype`) then 
      echo '<h2 id="'${Ztype}'"><A name="EB"><FONT color="Black">Comparisons for type '${Ztype}'</FONT></A><BR></h2>' >> ${PlotDir}/index.html
      foreach physplot ( Z0_Y Z0_Pt Z0_mass e1_P_t e1_eta e1_phi e2_P_t e2_eta e2_phi )
         echo '<A href="#'${physplot}${Ztype}'"><FONT color="Black">'${physplot}', </FONT></A>' >> ${PlotDir}/index.html
      end

      foreach physplot ( Z0_Y Z0_Pt Z0_mass e1_P_t e1_eta e1_phi e2_P_t e2_eta e2_phi )
         echo '<h4 id="'${physplot}${Ztype}'"><A name="EB"><FONT color="Blue">Comparisons for Variable '${physplot}'</FONT></A><BR></h4>' >> ${PlotDir}/index.html
         foreach cuts ( C09 C08 C07 C06 C05 C04 C03 )
            foreach FromData (`/bin/ls ${PlotDir} |grep "png" |grep $Ztype |grep $physplot | grep $cuts`)
            echo '<A HREF='${FromData}'> <img height="300" src="'${FromData}'"> </A>' >> ${PlotDir}/index.html
            end
         end
      end
   endif
end


if (`echo $plottype | grep -c "eps"` > 0 ) then
foreach FromData (`/bin/ls ${PlotDir} |grep ".eps"`)
   ~/eps2pdf.perl ${PlotDir}/$FromData
end

mv *pdf ${PlotDir}/.

foreach FromData (`/bin/ls ${PlotDir}/.. |grep ".eps"`)
   ~/eps2pdf.perl ${PlotDir}/../$FromData
end

mv *pdf ${PlotDir}/../.

endif



mv  ${PlotDir}/Z_CompareFINAL* ${PlotDir}/../FinalCompare/.


set PlotDir=${oldPlotDir}/FinalCompare
echo $PlotDir
cat > ${PlotDir}/index.html <<EOF
<HTML>

<HEAD><TITLE>ZShape MC to Data FINAL Comparison WebPages for `date`</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
<h2><A name="EB"><FONT color="Black">ZShape Web Pages for `date`</FONT></A><BR></h2>
<h3><A name="EB"><FONT color="Blue">Final Comparison Histograms</FONT></A><BR></h3>

EOF

foreach Ztype ( Tight-ECAL-Loose-ECAL Tight-ECAL-HF )
   echo '<h2 id="'${Ztype}'"><A name="EB"><FONT color="Black">Comparisons for type '${Ztype}'</FONT></A><BR></h2>' >> ${PlotDir}/index.html
   foreach FromData (`/bin/ls ${PlotDir} |grep ".png" | grep "Z_CompareFINAL" | grep $Ztype` )
      echo '<A HREF='${FromData}'> <img height="400" src="'${FromData}'"> </A>' >> ${PlotDir}/index.html
   end
end


#end of ZShape web page making script
