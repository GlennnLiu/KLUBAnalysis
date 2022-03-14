#!/bin/bash
selections="s1b1jresolvedMcut s2b0jresolvedMcut sboostedLLMcut VBFloose"
masses="250 260 270 300 320 350 400 450 500 550 600 650 700 750 800 850 900 1000 1250 1500 1750 2000 2500 3000"
#masses="450 500 550 600 650 700 750 800 850"
mkdir submit_get_limits_res_combined1
cd submit_get_limits_res_combined1
for mass in $masses
do
f=get_limits_res_combined${mass}.sh
echo '#!/bin/bash' > $f
echo 'mass="'$mass'"' >>$f
echo 'var="HHKin_mass_raw"' >>$f
echo 'tag="cards_Legacy2018_11Oct2021"' >>$f
echo 'identifier=".test"' >>$f
#echo 'cd ..' >>$f
#echo 'mkdir $tag/combined_out' >>$f
echo 'echo $mass' >>$f
echo 'combine -M AsymptoticLimits /data_CMS/cms/liugeliang/HHbbtautau_limits/$tag/comb.Radion$mass.root -n $identifier --run blind --noFitAsimov -m $mass --freezeParameters SignalScale > /data_CMS/cms/liugeliang/HHbbtautau_limits/$tag/combined_out/comb.Radion$mass.log' >>$f
chmod +x $f
t3submit $f
done
