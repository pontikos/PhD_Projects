#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


function usage() {
    echo "syntax: $0"
    echo " -b : basedir"
    echo " -f : infile"
    echo " -d : infile"
    exit 1
}



alias qsub="qsub -cwd -b y"

function density_threshold() {
    INFILE=$1
    D=$2
    rm --force $basedir/Out/*
    rm --force $basedir/Err/*
    for x in `cat $INFILE`
        do
        chmod ug+rx ~nikolas/bin/FCS/bin/density.filter.R
        y=`basename $x`
        echo qsub -N J$y -cwd -b y -u nikolas -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out ~nikolas/bin/FCS/bin/density.filter.R --in.file $x --plot.dir ${basedir}/Plot --density.threshold $D
        qsub -N J$y -cwd -b y -u nikolas -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out ~nikolas/bin/FCS/bin/density.filter.R --in.file $x --plot.dir ${basedir}/Plot --density.threshold $D
    done
}

function main() {
    basedir=~/dunwich/IL2/DensityThreshold
    infile=
    density="50%"
    #parse the args
    while getopts "b:f:d:" optionName
    do
        case "$optionName" in
        b) basedir=$OPTARG;;
        f) infile=$OPTARG;;
        d) density=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [ "$infile" = "" ]
    then
        usage
    fi
    basedir=$basedir
    echo $basedir
    mkdir -p $basedir
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err
    density_threshold $infile $density
}

### MAIN
main $*
