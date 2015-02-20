#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob

BIN=~nikolas/Projects/IL2/bin/CellPhenotypes/pstat5_pheno_test.R 

function usage() {
    echo "syntax: $0"
    echo " -K : number of flowClust clusters"
    echo " -e : extension"
    echo " -b : basedir"
    echo " -h : prints this message"
    exit 1
}


# args are set in README of each project
function main() {
    K=
    #parse the args
    while getopts "K:e:b:o:fh" optionName
    do
        case "$optionName" in
        K) K=$OPTARG;;
        e) ext=$OPTARG;;
        b) basedir=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [[ "$K" == "" ]]
    then
        usage
    fi
    basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/flowClust-$K
    indir=$basedir/CSV
    pstat5pheno=$basedir/flowclust-$K-$ext-pheno.csv
    if [ ! -e $pstat5pheno ]
    then
        cat $indir/*$ext.csv | sort -ru > $pstat5pheno
    fi
    infile=$pstat5pheno
    plotfile=$basedir/flowclust-$K-$ext-pheno-summary.pdf
    outfile=$basedir/flowclust-$K-$ext-pheno-summary.csv
    Rscript $BIN --in.file $infile --plot.file $plotfile --out.file $outfile
}

### MAIN
main $*

