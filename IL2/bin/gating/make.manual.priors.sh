
BASEDIR=~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3
outdir=$BASEDIR/Prior.RData/
mkdir -p $outdir
Rscript ~nikolas/Projects/IL2/bin/make.manual.priors.R  --clr $BASEDIR/CLR/I022267C_CB00010K_0U.clr --in.file $BASEDIR/pstat5-join/All/CB00010K_2012-11-13.RData --out.dir $outdir

