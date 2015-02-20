### Identify the peaks in the PSTAT5 distribution

#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


function usage() {
    echo "syntax: $0"
    echo " -b : basedir [default pstat5-peaks]"
    echo " -f : file containing list of FCS files."
    echo " -h : prints this message"
    exit 1
}

alias qsub="qsub -cwd -b y"

function peaks() {
    INFILE=$1
    rm --force $basedir/Out/*
    rm --force $basedir/Err/*
    chmod u+rx ~nikolas/Projects/IL2/bin/pstat5.peaks.R 
    for x in `cat $INFILE`
    do
        y=`basename $x`
        echo qsub -N $y -cwd -b y -u nikolas -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out ~nikolas/Projects/IL2/bin/pstat5.peaks.R -f $x --plot $basedir/Plot/${y%.fcs}.png
        qsub -N $y -cwd -b y -u nikolas -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out ~nikolas/Projects/IL2/bin/pstat5.peaks.R -f $x --plot $basedir/Plot/${y%.fcs}.png
        sleep 1
    done
    #cat ${basedir}/Out/*.out | sort -u | sed 's/>,//' | sed 's/>//' 
}



function main() {
    fcs=0
    basedir=`pwd`
    #parse the args
    while getopts "f:b:h" optionName
    do
        case "$optionName" in
        b) basedir=$OPTARG;;
        f) fcs=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [ "$fcs" = "0" ]
    then
        usage
    fi
    basedir=${basedir}/pstat5-peaks
    echo $basedir
    mkdir -p $basedir
    mkdir -p $basedir/RData
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err
    peaks $fcs

}

### MAIN
main $*

