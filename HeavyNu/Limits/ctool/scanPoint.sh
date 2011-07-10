#!/bin/sh

cmswd=$1
jobwd=$2
inp=$3
mass=$4
outp=$5

method=$6
toys=$7
special=$8
seed=1789

source /local/cms/sw/cmsset_default.sh

cd ${cmswd}

cmsenv

cd ${jobwd}

if [ "${method}" == "MarkovChainMC" ]; then
   extra="-i 30000 --tries 50"
   extraObs="-i 30000 --tries 100"
fi

localOF=/tmp/sp_$$

rm ${outp}

echo "Observed Limit Pass" > ${outp}
echo "==============================================" >> ${outp}

combine -v0  -n WRmuObs -m ${mass} -M ${method} -H ProfileLikelihood -s ${seed} ${extraObs} ${inp} > ${localOF} 2>&1 

cat ${localOF} >> ${outp}

echo " " >> ${outp}

rm ${localOF}

echo "Expected Limit Pass" >> ${outp}
echo "==============================================" >> ${outp}

combine -v0 -t${toys} -n WRmu -m ${mass} -M ${method} -H ProfileLikelihood -s ${seed} ${extra} ${inp} > ${localOF} 2>&1 

cat ${localOF} >> ${outp}

rm ${localOF}




#echo combine -v0 -t${toys} -n WRmu -m ${mass} -M ${method} -H ProfileLikelihood -s ${seed} ${extra} ${inp} 

#-m 16001000 -s1002 -n WRmu -t100 -v0  
