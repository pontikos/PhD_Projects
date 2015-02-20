#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob

function usage() {
    echo "syntax: $0"
    echo " -b : indir [default $indir]"
    echo " -K : [default $K]"
    echo " -d : outdir [default $outdir]"
    echo " -h : prints this message"
    exit 1
}

BIN=~nikolas/Projects/IL2/bin/CellPhenotypes/flowclust_pstat5_pheno.R

function pstat5_pheno() {
    chmod u+rx $BIN
    for f in $indir/C*.RData
    do
        x=`basename $f`
        individual=`echo $x | cut -f1 -d_`
        date=`echo $x | cut -f2 -d_`
        jobname=${individual}_${date}
        outfile=$outdir/`basename ${f%.RData}.csv`
        echo qsub -N $jobname -cwd -b y -e /dev/null -o /dev/null $BIN --in.file $f --out.file $outfile 
        qsub -N $jobname -cwd -b y -e /dev/null -o /dev/null $BIN --in.file $f --out.file $outfile 
        sleep 1
    done
}

function main() { 
    #parse the args
    while getopts "b:K:h" optionName
    do
        case "$optionName" in
        b) indir=$OPTARG;;
        K) K=$OPTARG;;
        d) outdir=$OPTARG;;
        ?) usage 0;;
        esac
    done
    indir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/flowClust-$K/RData
    outdir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/flowClust-$K/CSV
    echo indir: $indir
    mkdir -p $outdir
    echo outdir: $outdir
    echo K: $K
    pstat5_pheno

}

### MAIN
main $*


