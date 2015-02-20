#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3
basedir=~nikolas/dunwich/Projects/IL2/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5
#file containing individual and dates
fcs=${basedir}/individual-date.csv
#fcs=${basedir}/all-FCS-panel-date-ids-dose.csv
#inext=fcs
inext=RData
#indir=${basedir}/Lymphocytes6/RData
#indir=/chiswick/data/store/facs/Tony-FCS/PSTAT5/CD25/CD45RA/CD4/FOXP3/
indir=/chiswick/data/store/facs/Tony-RData/PSTAT5/CD25/CD45RA/CD4/FOXP3/
indir=~nikolas/dunwich/Projects/IL2/Tony-RData/CD25-CD3-CD4-CD45RA-CD56-CD8-FOXP3-PSTAT5/
outext=RData
#outdir=${basedir}/Lymphocytes/pstat5-join/
#outdir=${basedir}/pstat5-join/Lymphocytes6
outdir=${basedir}/pstat5-join/All
#channels on which to join
join=FSCW,SSCW,FSCH,SSCH,FSCA,SSCA,CD4,CD25,CD45RA,FOXP3
join=FSCW,SSCW,FSCH,SSCH,FSCA,SSCA,CD25,CD3,CD4,CD45RA,CD56,CD8,FOXP3,PSTAT5

function usage() {
    echo "syntax: $0"
    echo " -b : basedir [default $basedir]"
    echo " -f : file containing list of FCS files [default $fcs]"
    echo " -h : prints this message"
    exit 1
}

alias qsub="qsub -cwd -b y"
BIN=~nikolas/bin/ann-join2.R
BIN=~nikolas/bin/ann-join.R

function join() {
    echo $indir
    cd $indir
    chmod u+rx $BIN
    #go over individuals, dates and join the four doses files
    #for x in `cat $fcs | tail -n+2 | grep CB01495Z `
    for x in `cat $fcs`
    do
        individual=`echo $x | cut -f1 -d,`
        #individual=`echo $x | cut -f4 -d';' | cut -f3 -d,`
        date=`echo $x | cut -f2 -d,`
        #date=`echo $x | cut -f3 -d';'`
        jobname=${individual}_${date}
        echo qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}.err -o ${outdir}/Out/${individual}_${date}.out $BIN --in.files ${individual}_0U_${date}.${inext},${individual}_01U_${date}.${inext},${individual}_10U_${date}.${inext},${individual}_1000U_${date}.${inext} --out ${outdir}/${individual}_${date}.${outext} --join $join
        qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}.err -o ${outdir}/Out/${individual}_${date}.out $BIN --in.files ${individual}_0U_${date}.${inext},${individual}_01U_${date}.${inext},${individual}_10U_${date}.${inext},${individual}_1000U_${date}.${inext} --out ${outdir}/${individual}_${date}.${outext} --join $join
        sleep 1
    done
}


function main() { 
    #parse the args
    while getopts "i:j:e:b:o:f:h" optionName
    do
        case "$optionName" in
        e) ext=$OPTARG;;
        j) join=$OPTARG;;
        b) basedir=$OPTARG;;
        i) indir=$OPTARG;;
        o) outdir=$OPTARG;;
        f) fcs=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [ "$fcs" = "0" ]
    then
        usage
    fi
    basedir=${basedir}
    echo $basedir
    #
    echo $outdir
    mkdir -p $outdir
    mkdir -p $outdir/Out
    mkdir -p $outdir/Err
    join
}

### MAIN
main $*

