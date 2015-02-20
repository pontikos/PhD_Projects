
#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob


function usage() {
    echo "syntax: $0"
    echo " -K : number of clusters"
    echo " -d : downsample"
    echo " -p : make prior"
    echo " -g : gate"
    echo " -b : build lymph count file"
    echo " -h : prints this message"
    exit 1
}



alias qsub="qsub -cwd -b y"

function gate_lymph() {
    K=$1
    downsample=$2
    rm --force $basedir/Out/*
    rm --force $basedir/Err/*
    for x in `cat ~nikolas/IL2RA/fcsFile.csv`
        do
        chmod ug+rx ~nikolas/bin/FCS/bin/gate.lymphocytes.R
        y=`basename $x`
        qsub -N $y -cwd -b y -u nikolas -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out ~nikolas/bin/FCS/bin/gate.lymphocytes.R --in.file $x --out.dir ${basedir}/RData --plot.dir ${basedir}/Plot --prior ${basedir}/RData/prior.RData --post.threshold .5 --down.sample $downsample -K $K
    done
}


function build() {
    K=$1
    outfile=$basedir/Out/lymph.count$K.csv
    echo "fcsFile,lymph.count$K" > $outfile
    cat $basedir/Out/*.out | grep lymph.count | cut -d'/' -f9 | cut -d, -f1,3 | sort >> $outfile
    echo $outfile
    Rscript ~nikolas/IL2RA/bin/compare.lymph.count.R --file1 ~/IL2RA/CellPhenotypes/manual.csv --file2 $outfile --out.file $basedir/manual-agreement.pdf
}


function build2() {
    K=$1
    outfile=$basedir/Out/lymph.count.mix$K.csv
    echo "fcsFile,lymph.count$K" > $outfile
    for f in $basedir/RData/*.RData
    do
        echo $f
        x=`echo "load('${f}'); Sys.sleep(1); count <- round(dim(d)[1]*res@w[which.max(res@mu[,3])]); print(count);" | R --no-save --no-restore --silent | grep '[1]'`
        x=`echo $x | awk '{print($(NF))}'`
        echo `basename "${f%.RData}"`,$x >> $outfile
    done
    echo $outfile
    #need to delete last line since it contains some artefact
    head -n -1 $outfile > $outfile
    Rscript ~nikolas/IL2RA/bin/compare.lymph.count.R --file1 ~/IL2RA/CellPhenotypes/manual.csv --file2 $outfile --out.file $basedir/manual-agreement-mix.pdf
}




function main() {
    K=
    downsample=0.1
    prior=0
    gate=0
    build=0

    #parse the args
    while getopts "K:d:pgb:h" optionName
    do
        case "$optionName" in
        K) K=$OPTARG;;
        d) downsample=$OPTARG;;
        p) prior=1;;
        g) gate=1;;
        b) build=$OPTARG;;
        ?) usage 0;;
        esac
    done

    basedir=~/dunwich/IL2RA/Lymphocytes$K
    echo $basedir
    mkdir -p $basedir
    mkdir -p $basedir/RData
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err

    if [ "$prior" = 1 ]
    then
        make_prior $K
    elif [ "$gate" = 1 ]
    then
        gate_lymph $K $downsample
    elif [ "$build" = 1 ]
    then
        build $K
    elif [ "$build" = 2 ]
    then
        build2 $K
    fi
}

### MAIN
main $*
