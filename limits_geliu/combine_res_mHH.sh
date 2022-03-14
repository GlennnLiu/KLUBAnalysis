var="DNNoutSM_kl_1"
selections="s1b1jresolvedMcut s2b0jresolvedMcut sboostedLLMcut VBFloose"

mHHs="mHH250_335 mHH335_475 mHH475_725 mHH725_1100 mHH1100_3500 InvDR"

masses="250 260 270 280 300 320 350 400 450 500 550 600 650 700 750 800 850 900 1000 1250 1500 1750 2000 2500 3000"
#masses="300 320 400 600 850 1250"


tags="cards_MuTauLegacy2018_24May2021 cards_ETauLegacy2018_24May2021"
for tag in $tags
do
cd $tag
for sel in $selections
do
mkdir $sel
for mass in $masses
do
echo $mass
#combineCards.py -S ${sel}InvDR$var/hhres*.Radion$mass.txt >> ${sel}InvDR$var/comb.Radion$mass.txt
combineCards.py -S ${sel}*$var/hhres*.Radion$mass.txt > ${sel}/comb.Radion$mass.txt
echo 'SignalScale rateParam * Radion'$mass' 0.01' > add.txt
cat add.txt >> ${sel}/comb.Radion$mass.txt
text2workspace.py ${sel}/comb.Radion$mass.txt -o ${sel}/comb.Radion$mass.root
done
done
cd -
done


mHHs="mHH250_335 mHH335_475 mHH475_725 mHH725_1100 mHH1100_3500"

tags="cards_TauTauLegacy2018_24May2021"

for tag in $tags
do
cd $tag
for sel in $selections
do
mkdir $sel
for mass in $masses
do
echo $mass
#combineCards.py -S ${sel}InvDR$var/hhres*.Radion$mass.txt >> ${sel}InvDR$var/comb.Radion$mass.txt
combineCards.py -S ${sel}*$var/hhres*.Radion$mass.txt > ${sel}/comb.Radion$mass.txt
echo 'SignalScale rateParam * Radion'$mass' 0.01' > add.txt
cat add.txt >> ${sel}/comb.Radion$mass.txt
text2workspace.py ${sel}/comb.Radion$mass.txt -o ${sel}/comb.Radion$mass.root
done
done
cd -
done

