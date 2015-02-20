#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


#ext=pstat5-norm
ext=pstat5
basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes
#outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/Naive
#outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/Memory 
outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/Lymphocytes/pstat5-ann-response
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

BIN=~nikolas/Projects/IL2/bin/pstat5-ann-response.R 
#BIN=~nikolas/Projects/IL2/bin/pstat5-response.R 

function pstat5response() {
    chmod u+rx $BIN
    for x in `cat $fcs`
    do
        individual=`echo $x | cut -f1 -d,`
        date=`echo $x | cut -f2 -d,`
        #grep $individual.*$date PSTAT5-CD25-CD45RA-CD4-FOXP3.csv
        jobname=${individual}_${date}_${ext}
        echo qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}_${ext}.err -o ${outdir}/Out/${individual}_${date}_${ext}.out $BIN --in.dir ${basedir}/RData --plot.dir ${outdir}/Plot --out.dir ${outdir}/RData --individual $individual --date $date --ext $ext --parameter pSTAT5 --doses 0U,01U,10U,1000U
        qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err/${individual}_${date}_${ext}.err -o ${outdir}/Out/${individual}_${date}_${ext}.out $BIN --in.dir ${basedir}/RData --plot.dir ${outdir}/Plot --out.dir ${outdir}/RData --individual $individual --date $date --ext $ext --parameter pSTAT5 --doses 0U,01U,10U,1000U
        sleep 1
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

