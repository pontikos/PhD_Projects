#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob

#channels=CD25,CD45RA,CD4,FOXP3,FSCA,SSCA
channels=CD25,CD45RA,CD4,FOXP3
channels=FSCW,SSCW,FSCH,SSCH,FSCA,SSCA,CD4,CD25,CD45RA,FOXP3,PSTAT5.1,PSTAT5.2,PSTAT5.3,PSTAT5.4
basedir=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3
lymphdir=${basedir}/Lymphocytes/
outdir=${basedir}/CellPhenotypes/
downsample=.95
fcs=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/lymphocytes-all.csv
fcs=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/lymphocytes.csv
pooldata=$lymphdir/pooled-data-`basename ${fcs%.csv}`.RData

function usage() {
    echo "syntax: $0"
    echo " -b : basedir [default $basedir]"
    echo " -d : outdir [default $outdir]"
    echo " -f : file containing list of FCS files [default $fcs]"
    echo " -h : prints this message"
    exit 1
}

BIN=~nikolas/Projects/IL2/bin/gating/flowclust.R 


function make_pool() {
    plotfile=$lymphdir/pooled-data-`basename ${fcs%.csv}`.png
    Rscript ~nikolas/bin/FCS/pool-data.R --in.file $fcs --channels $channels -T 100000 --plot.file $plotfile --out.file $pooldata
}


function make_prior() {
    downsample=$1
    rm --force $outdir/Err.qsub/prior.err
    rm --force $outdir/Out.qsub/prior.out
    rm --force $outdir/Plot/prior.fcs.png
    chmod ug+rx $BIN
    Rscript $BIN --in.file $pooldata --out.dir $outdir/RData --plot.dir $outdir/Plot --post.threshold .8  -K $K  
}


function gate() {
    fcs=$1
    rm --force $basedir/Out.qsub/C*
    rm --force $basedir/Err.qsub/C*
    chmod u+rx $BIN
    for x in `cat $fcs`
    do
        y=`basename $x`
        individual=`echo $y | cut -f1 -d_`
        date=`echo $y | cut -f2 -d_`
        jobname=${individual}_${date}
        echo qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err.qsub/${individual}_${date}.err -o ${outdir}/Out.qsub/${individual}_${date}.out $BIN --prior ${prior} --in.file $x  --out.dir ${outdir}/RData/ --plot.dir ${outdir}/Plot/ --post.threshold .95 --channels $channels
        qsub -N $jobname -cwd -b y -u nikolas -e ${outdir}/Err.qsub/${individual}_${date}.err -o ${outdir}/Out.qsub/${individual}_${date}.out $BIN --prior ${prior} --in.file $x  --out.dir ${outdir}/RData/ --plot.dir ${outdir}/Plot/ --post.threshold .95 --channels $channels
        sleep 1
    done
}



# args are set in README of each project
function main() {
    #parse the args
    while getopts "K:e:b:o:f:h" optionName
    do
        case "$optionName" in
        K) K=$OPTARG;;
        b) basedir=$OPTARG;;
        o) outdir=$OPTARG;;
        f) fcs=$OPTARG;;
        ?) usage 0;;
        esac
    done
    echo $fcs
    if [ "$fcs" = "0" ]
    then
        usage
    fi
    outdir=$outdir/flowClust-$K
    prior=$outdir/RData/prior.RData
    mkdir -p $outdir/
    mkdir -p $outdir/RData
    mkdir -p $outdir/Out.qsub
    mkdir -p $outdir/Err.qsub
    mkdir -p $outdir/Plot
    if [[ "$K" -le 0 || "$basedir" == "" ]]
    then
        usage
    fi
    echo basedir: $basedir
    echo outdir: $outdir
    echo infile: $fcs
    #check that all files in $fcs are valid
    if [[ "$fcs" != "" ]]
    then
        if [[ ! -e $fcs ]]
        then
            echo $fcs file does not exist!
            exit 1
        fi
        while read -r f
        do
            if [[ ! -e $(eval echo $f) ]]
            then
                echo $f in $fcs does not exist!
                exit 1
            fi
        done < $fcs
    fi
    #make pool file
    if [[ ! -e $pooldata ]]
    then
        echo make_pool $fcs
        make_pool $fcs
    else
        echo $pooldata exists
    fi
    if [ ! -e $prior ]
    then
        echo $prior does not exist!
        echo will create prior
        n=`wc -l $fcs | cut -f1 -d ' '`
        echo $n
        downsample=`echo "scale=0; 10000/$n" | bc -l`
        make_prior $downsample
    elif [ "$fcs" != "" ]
    then
        echo $prior exists
        gate $fcs
    fi
}

### MAIN
main $*
