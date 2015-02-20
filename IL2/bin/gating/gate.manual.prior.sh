#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob

GATE=gate.manual.prior.R 
GATE=`dirname $0`/$GATE

PRIOR_GATES=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/PriorGates

function usage() {
    echo "syntax: $0"
    echo " -b : <basedir> Basedir under which to write Lymphocytes/Out, Lymphocytes/Err, Lymphocytes/Plot and Lymphocytes/RData"
    echo " -c : <channels> Channels comma separated."
    echo " -h : Prints this message."
    exit 1
}

alias qsub="qsub -cwd -b y"

function gate() {
    G=$1
    echo $G
    #rm --force $basedir/Out/*
    #rm --force $basedir/Err/*
    for x in $basedir/*.RData
        do
        chmod ug+rx $GATE
        y=`basename $x`
        echo qsub  -N J$y -cwd -b y -u $USER -e ${outdir}/Err/${y%.fcs}.err -o ${outdir}/Out/${y%.fcs}.out $GATE --in.file $x --prior $PRIOR_GATES/$G.prior.RData --plot.file ${outdir}/Plot/$y.png --out.file ${outdir}/$y --channels $channels
        qsub  -N J$y -cwd -b y -u $USER -e ${outdir}/Err/${y%.fcs}.err -o ${outdir}/Out/${y%.fcs}.out $GATE --in.file $x --prior $PRIOR_GATES/$G.prior.RData --plot.file ${outdir}/Plot/$y.png --out.file ${outdir}/$y --channels $channels
        sleep 1
    done
}


# args are set in README of each project
function main() {
    basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All
    gate=Lymphocytes.Single.Cells.CD4.Memory
    channels=FSCA,SSCA,CD4,CD25,CD45RA
    while getopts "b:g:c:h" optionName
    do
        case "$optionName" in
        b) basedir=$OPTARG;;
        g) gate=$OPTARG;;
        c) channels=$OPTARG;;
        ?) usage 0;;
        esac
    done
    outdir=$basedir/$gate
    echo basedir: $basedir
    echo $outdir
    echo gate: $gate
    mkdir -p $basedir
    mkdir -p $outdir/Plot
    mkdir -p $outdir/Out
    mkdir -p $outdir/Err
    gate $gate
}

### MAIN
main $*
