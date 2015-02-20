#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


#ext=pstat5
ext=pstat5-norm
basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes
#outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/mclust-pred
outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/mclust
fcs=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/individual-date.csv

function usage() {
    echo "syntax: $0"
    echo " -b : basedir [default $basedir]"
    echo " -d : outdir [default $outdir]"
    echo " -f : file containing list of FCS files [default $fcs]"
    echo " -e : extension [default $ext]"
    echo " -h : prints this message"
    exit 1
}

alias qsub="qsub -cwd -b y"
BIN=~nikolas/Projects/IL2/bin/mclust.R 

function pstat5response() {
    chmod u+rx $BIN
    for x in `cat $fcs`
    do
        individual=`echo $x | cut -f1 -d,`
        date=`echo $x | cut -f2 -d,`
        jobname=${individual}_${date}_${ext}
        #qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}_${ext}.err -o ${outdir}/Out/${individual}_${date}_${ext}.out $BIN --in.dir ${basedir}/All --out.dir ${outdir} --individual $individual --date $date --ext $ext --parameters CD25,CD45RA,CD4,FOXP3 --predict ~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/mclust/RData/all.RData
        qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}_${ext}.err -o ${outdir}/Out/${individual}_${date}_${ext}.out $BIN --in.dir ${basedir}/All --out.dir ${outdir} --individual $individual --date $date --ext $ext --parameters CD25,CD45RA,CD4,FOXP3 
    done
}


function main() { 
    #parse the args
    while getopts "e:b:o:fh" optionName
    do
        case "$optionName" in
        e) ext=$OPTARG;;
        b) basedir=$OPTARG;;
        d) outdir=$OPTARG;;
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
    mkdir -p $basedir
    mkdir -p $basedir/RData
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err
    #
    echo $outdir
    mkdir -p $outdir
    mkdir -p $outdir/RData
    mkdir -p $outdir/Plot
    mkdir -p $outdir/Out
    mkdir -p $outdir/Err
    pstat5response
}

### MAIN
main $*

