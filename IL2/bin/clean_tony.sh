set -x

cb=`tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/cb-IL2RA.csv | cut -f2 -d,`
dgap=`tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/dgap-IL2RA.csv | cut -f2 -d, `

#for x in $cb $dgap
for x in $dgap
do
    mkdir -p ~nikolas/dunwich/Projects/IL2/clean.Tony/$x
    fcs=`find /chiswick/data/store/facs/Tony/Tony -iname *$x*.fcs`
    for f in $fcs
    do
        f2=~nikolas/dunwich/Projects/IL2/clean.Tony/$x/`echo $f | rev | cut -d'/' -f 1,2 | rev | tr '/' '-'`
        cp $f $f2
    done
done


for f in ~nikolas/dunwich/Projects/IL2/clean.Tony/*/*.fcs
do
    python ~/bin/fcstools/getPanel.py $f
done

#find /chiswick/data/store/facs/Tony/Tony/21112013_normal*  -iname *.fcs -print0 | while read -d $'\0' f
#find /chiswick/data/store/facs/Tony/Tony/20130215_SEB_24h_Abatacept_25OHVitD3/CB01566Q_SEB* -iname *.fcs -print0 | while read -d $'\0' f

find /chiswick/data/store/facs/Tony/Tony/ -iname *.fcs -print0 | while read -d $'\0' f
do
    panel=`python ~/bin/fcstools/getPanel.py "$f" | tr '\t' '-' | sed 's/\s\+/-/g' `
    date=`python ~/bin/fcstools/getDate.py "$f"`
    echo "$f;$panel;$date"
done > ~nikolas/dunwich/Projects/IL2/tony-all-2014-01-07.csv



cb=`tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/cb-IL2RA.csv | cut -f2 -d,`
dgap=`tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/dgap-IL2RA.csv | cut -f2 -d, `


cb=`tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/cb-IL2RA.csv | cut -f1 -d,`


for x in `tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/cb-IL2RA.csv | cut -f1,2,3 -d,` `tail -n+2 ~nikolas/Projects/IL2/IL2RA-genotype/dgap-IL2RA.csv | cut -f1,2,3 -d,` 
do
    x1=`echo $x | cut -f1 -d,`
    grep $x1 ~nikolas/dunwich/Projects/IL2/tony-all-2014-01-07.csv | xargs -i echo {}";$x"
    x2=`echo $x | cut -f2 -d,`
    grep $x2 ~nikolas/dunwich/Projects/IL2/tony-all-2014-01-07.csv | xargs -i echo {}";$x"
    x3=`echo $x | cut -f3 -d,`
    grep $x3 ~nikolas/dunwich/Projects/IL2/tony-all-2014-01-07.csv | xargs -i echo {}";$x"
done | sort -u > ~nikolas/dunwich/Projects/IL2/tony-cb-dgap-2014-01-07.csv



for x in `grep ';CD25,CD4,CD45RA,FOXP3,PSTAT5;' cb-FCS-IL2RA.csv` 
do
    individual=`echo $x | cut -d';' -f1 | cut -d, -f1`
    date=`echo $x | cut -d';' -f4`
    dose=`echo $x | cut -d';' -f5`
    fcs=`echo $x | cut -d';' -f2`
    ln -s $fcs /chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/${individual}_${dose}_${date}.fcs
done


for x in `grep ';CD25,CD4,CD45RA,FOXP3,PSTAT5;' cb-FCS-IL2RA.csv` 
do
    individual=`echo $x | cut -d';' -f1 | cut -d, -f1`
    date=`echo $x | cut -d';' -f4`
    dose=`echo $x | cut -d';' -f5`
    fcs=`echo $x | cut -d';' -f2`
    ln -s $fcs /chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/${individual}_${dose}_${date}.fcs
done

for x in `grep ';CD25,CD3,CD4,CD45RA,CD56,FOXP3,PSTAT5;' ~/dunwich/Projects/IL2/tony-cb-dgap-2014-01-07.csv`;
do
    fcs=`echo $x | cut -d';' -f1`
    date=`echo $x | cut -d';' -f3`
    individual=`echo $x | cut -d';' -f4 | cut -d, -f3`
    #dose=`echo $x | cut -d';' -f5`
    echo $x
    #echo ln -s $fcs /chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/${individual}_${dose}_${date}.fcs
    break
done


for x in  `tail -n+2 ~/dunwich/Projects/IL2/pstat5-il2-cb-dgap.csv`
do 
    fcs=`echo $x | cut -d';' -f1`
    panel=`echo $x | cut -d';' -f2 | sed 's/,/-/g' `
    date=`echo $x | cut -d';' -f3`
    individual=`echo $x | cut -d';' -f4 | cut -d, -f3`
    dose=`echo $x | cut -d';' -f5`
    echo $panel $date $dose
    mkdir -p  /chiswick/data/store/facs/Tony-FCS/$panel/
    echo ln -s $fcs /chiswick/data/store/facs/Tony-FCS/$panel/${individual}_${dose}_${date}.fcs
    rm /chiswick/data/store/facs/Tony-FCS/$panel/${individual}-${dose}-${date}.fcs
done



