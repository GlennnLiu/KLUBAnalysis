#!/bin/bash
selections="s1b1jresolvedMcut s2b0jresolvedMcut sboostedLLMcut VBFloose"
masses="250 260 270 280 300 320 350 400 450 500 550 600 650 700 750 800 850 900 1000 1250 1500 1750 2000 2500 3000"
var="DNNoutSM_kl_1"


#tag="cards_TauTau2016_22Sep2020"
#tag="cards_MuTau2016_22Sep2020"
#tag="cards_ETau2016_22Sep2020"

#tag="cards_TauTau2017_22Sep2020"
#tag="cards_MuTau2017_22Sep2020"
#tag="cards_ETau2017_22Sep2020"

##tag="cards_TauTau2018_22Sep2020"
##tag="cards_MuTau2018_22Sep2020"
##tag="cards_ETau2018_22Sep2020"

#tag="cards_TauTau2018_res_17Feb2021"
#tag="cards_MuTau2018_res_17Feb2021"
#tag="cards_ETau2018_res_17Feb2021"

#tag="cards_2018_17Feb2021"
#tag="cards_CombChan2018_res_17Feb2021"

#tags="cards_TauTau2018_res_constbin_19Feb2021 cards_MuTau2018_res_constbin_19Feb2021 cards_ETau2018_res_constbin_19Feb2021"
tag="cards_CombChan2018_res_Legacy2018_29Mar2021_m"

#combine -M AsymptoticLimits cards_TauTau2018_res_15Feb2021/s1b1jresolvedMcutDNNoutSM_kl_1/hh_2018_2_C1_13TeV.m450.txt -n .test210215 --run blind --noFitAsimov

identifier=".test"

echo "TEST"

#for tag in $tags
#do
#    for sel in $selections
#    do 
#	mkdir $tag/$sel$var/combined_out
#	echo $sel
#	for mass in $masses
#	do
#	    echo $mass
#	    cd $tag #/$sel$var/combined_out
#	    #combine -M AsymptoticLimits ../comb.Radion$mass.txt -n $identifier --run blind --noFitAsimov -m $mass > comb.Radion$mass.log	
#	    combine -M AsymptoticLimits $sel$var/comb.Radion$mass.txt -n $identifier$tag$sel --run blind --noFitAsimov -m $mass > $sel$var/combined_out/comb.$sel.Radion$mass.log
#	    
#	    cd -
#	done
#    done
#done
cd $tag
for sel in $selections
do
echo $sel
mkdir $sel$var/combined_out
for mass in $masses
do
    echo $mass
    combine -M AsymptoticLimits $sel$var/comb.Radion$mass.txt -n $identifier --run blind --noFitAsimov --freezeParameters SignalScale -m $mass > $sel$var/combined_out/comb.Radion$mass.log
done
done
