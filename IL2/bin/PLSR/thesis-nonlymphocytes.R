#### Non-lymphocytes

print(f <- file.path(base.dir,'RData','pstat5-join', sprintf('%s_%s.RData',individual,date)))
print(load(f))
fcs.data <- baseline.relative.pstat5(fcs.data)
#already done
#fcs.data <- transform.scatter(applyTransforms(fcs.data,transforms))
print(load(file.path(base.dir,'CLR','CB00086S_0U_2012-09-18.RData')))
MARKERS <- CORE.MARKERS

# cells identified at lowest dose
par(mfrow=c(1,1))
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',2), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y,outliers=T)
# purple
g <- structure(list(x = c(-3.42781386640277, -2.83331075199373, -2.23880763758468, -2.29825794902559, -3.04138684203689, -4.02231698081181, -4.73572071810266, -4.94379680814583, -4.46819431661859, -3.54671448928458, -3.2197377763596), y = c(0.876407526174889, 0.876407526174889, 0.614897660608943, 0.291856061968656, 0.0303461964027103, 0.153409662551391, 0.322621928505827, 0.584131794071773, 0.799492859831964, 0.922556325980644, 0.861024592906304)), .Names = c("x", "y"))
lines(g$x,g$y,col='purple')
j <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))
# pink
#dput(g <- locator(type='l',lwd=2, col='pink'))
#g <- structure(list(x = c(9.59180433915526, 8.46224842177809, 7.68939437304633, 7.89747046308949, 8.61087420038035, 9.35400309339165, 9.79988042919843, 10.3349332321666, 10.097131986403, 9.59180433915526), y = c(-0.000419670134459622, 0.445685394654507, 0.799492859831964, 1.07638565866649, 0.953322192517815, 0.707195260220453, 0.430302461385922, 0.153409662551391, -0.0773343364773852, -0.0158026034030448)), .Names = c("x", "y"))
g <- structure(list(x = c(7.64296805941769, 7.87307347606061, 8.75514423985847, 9.31762414720784, 9.99515676287866, 10.2636130822954, 9.91845495733102, 9.35597504998166, 8.67844243431083, 7.66853532793357), y = c(-0.195731365610558, 0.735728329573008, -0.00301556660706184, -0.725699812870173, -1.27172791004675, -1.62504020821982, -1.84987530705724, -1.81775600722332, -1.23960861021283, -0.195731365610558)), .Names = c("x", "y"))
lines(g$x,g$y,col='pink')
k <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y)) 
# lightblue
#dput(g <- locator(type='l',lwd=2, col='lightblue'))
#before fixing jon
#g <- structure(list(x = c(2.99281976921489, 4.62770333383976, 5.57890831689423, 6.35176236562598, 7.18406672579864, 7.30296734868045, 6.64901392283051, 5.81670956265785, 5.25193160396925, 4.5385278666784, 3.854849285108, 3.08199523637625, 2.33886634336495, 1.684912917515, 1.0012343359446, 1.0012343359446, 1.80381354039681, 2.30914118764449, 2.96309461349444), y = c(0.799492859831964, 0.99947099232357, 1.01485392559215, 0.968705125786399, 0.830258726369134, 0.368770728311582, 0.0611120629398807, -0.185014869357481, -0.384993001849086, -0.538822334534937, -0.508056467997767, -0.569588201072107, -0.554205267803522, -0.461907668192012, -0.384993001849086, -0.0465684699402148, 0.384153661580167, 0.676429393683283, 0.799492859831964)), .Names = c("x", "y"))
#after fixing join
g <- structure(list(x = c(6.41573917065544, 5.52088477259964, 5.02232303653998, 4.11468500422624, 2.96415792101163, 2.10765442572965, 1.83919810631291, 1.83919810631291, 3.00250882378545, 3.43715238855541, 5.02232303653998, 6.18563375401252, 6.73533002710395, 6.86316636968335, 6.81203183265159, 6.44130643917132), y = c(0.607251130237343, 1.66718802475657, 1.76354592425832, 1.76354592425832, 1.57083012525482, 1.26569677683262, 0.96056342841042, 0.543012530569511, -0.452685764281887, -0.565103313700594, -0.6775208631193, -0.452685764281887, -0.1796717156936, -0.0190752165240198, 0.446654631067762, 0.623310780154301)), .Names = c("x", "y"))
lines(g$x,g$y,col='lightblue')
l <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))
# orange
g <- structure(list(x = c(-1.37677812169157, -0.782275007282528, -0.128321581432581, -0.068871269991677, -0.722824695841624, -1.43622843313247, -1.88210576893926, -1.94155608038016, -1.91183092465971, -1.46595358885293, -1.28760265453021), y = c(-0.323461268774746, -0.308078335506161, -0.523439401266352, -1.20028846508409, -1.38488366430712, -1.36950073103853, -1.01569326586107, -0.815715133369468, -0.492673534729182, -0.323461268774746, -0.338844202043331)), .Names = c("x", "y"))
lines(g$x,g$y,col='orange')
m <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

# cell identfied at highest dose
par(mfrow=c(1,1))
fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',4), paste(MARKERS, collapse='+'), sep='~'))
plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
Y <- p$scores[,1:2]
smoothPlot(Y,outliers=T)
#ywllow
g <- structure(list(x = c(3.90086813263438, 3.51311752200631, 3.34985410700502, 3.39066996075534, 3.63556508325728, 4.10494740138599, 4.65596142701535, 4.737593134516, 4.61514557326503, 4.43147423138858, 4.2069870357618, 3.92127605950954, 3.8396443520089, 3.6763809370076, 3.6763809370076), y = c(6.32942929401107, 5.61260526180101, 4.6316881650925, 3.80168139095453, 3.57531590709872, 4.06577445545298, 4.93350881023358, 5.72578800372891, 6.36715687465371, 6.63124993915215, 6.59352235850952, 6.32942929401107, 6.17851897144054, 5.87669832629946, 5.87669832629946)), .Names = c("x", "y"))
lines(g$x,g$y,col='yellow')
n <- as.logical(point.in.polygon(Y[,1],Y[,2],g$x,g$y))

#identification of pink and blue high density clusters
print(load(file.path(base.dir,'CLR',FILES[[1]])))
CLR <- cbind(CLR,PLSR.1=as.numeric(n))
CLR <- cbind(CLR,PLSR.2=as.numeric(m))
CLR <- cbind(CLR,PLSR.3=as.numeric(l))
CLR <- cbind(CLR,PLSR.4=as.numeric(k))
CLR <- cbind(CLR,PLSR.5=as.numeric(j))

#################### dose response
pdf('~/Thesis/figures/plsr-nonlymphocytes.pdf')
figure.labels <- iter(paste(letters,')',sep=''))
par(mfrow=c(2,2),mai=c(.75,.75,.5,.1))
cols <- c('black','purple','pink','lightblue','orange','yellow')
for (i in 2:4) {
    fm <- as.formula(paste(sprintf('diff.PSTAT5.%d',i), paste(MARKERS, collapse='+'), sep='~'))
    plsr(fm, data=as.data.frame(fcs.data),scale=FALSE)->p
    Y <- p$scores[,1:2]
    smoothPlot(Y,posteriors=CLR[,c('Lymphocytes',paste('PLSR',1:5,sep='.'))],clusters.col=cols,chulls=FALSE,outliers=TRUE,ellipse.lwd=3)
    title(paste(nextElem(figure.labels),DOSES[[i]], sep='\t'), adj=0)
}
plot.new()
legend('center',legend=c('Lymphocytes',paste('PLSR',1:5,sep='.')),fill=cols,cex=2,bty='n')
dev.off()


#################### dose response
pdf('~/Thesis/figures/plsr-nonlymphocytes-dose-response.pdf')
ylim <- range(sapply(c('Lymphocytes',paste('PLSR',1:5,sep='.')), function(cell.type) colMedians(fcs.data[as.logical(CLR[,cell.type]),paste('diff','PSTAT5',1:4,sep='.')]) ))
plot(NULL, xlim=c(0,3), ylim=ylim, xaxt='n', xlab='dose', ylab='pSTAT5 MFI')
axis(1, at=0:3, labels=DOSES)
i <- 1
cols <- c('black','purple','pink','lightblue','orange','yellow')
for (cell.type in c('Lymphocytes',paste('PLSR',1:5,sep='.'))) {
    mfi <- fcs.data[which(as.logical(CLR[,cell.type])),paste('diff.PSTAT5',1:4,sep='.')]
    lines(0:3, colMedians( mfi ), col=cols[[i]], lwd=3)
    i <- i+1
} 
legend('topleft', c('Lymphocytes',paste('PLSR',1:5,sep='.')), fill=cols, bty='n')
dev.off()

#################### univariate cluster plots
#univariate view
#yellow and pink are interesting
pdf('~/Thesis/figures/plsr-nonlymphocytes-clusters.pdf')
#univariate view
r <- c(1,1,2,2,3,3)
layout(rbind(r,3+r,6+r,rep(10,length(r))),heights=c(4,4,4,1)) 
par(mai=rep(0.5, 4)) 
for (marker in CORE.MARKERS)
smoothPlot1D(fcs.data[,marker],posteriors=CLR[,c('Lymphocytes','PLSR.1','PLSR.2')],chulls=FALSE,clusters.col=c('black','purple','pink'), outliers=TRUE,ellipse.lwd=3, main=marker,col='white')
par(mai=c(0,0,0,0))
plot.new()
legend('center',ncol=3,legend=c('Lymphocytes','PLSR.1','PLSR.2'),fill=c('black','purple','pink'),bty='n',cex=1.5)
dev.off()

# only keep the yellow population
# 'purple','pink','lightblue', 'orange', 'yellow'
CLR <- CLR[,c(!colnames(CLR) %in% c('PLSR.2','PLSR.3','PLSR.4','PLSR.5'))]
#now remove all those which have already been included as part of the lymphocytes
CLR[which(CLR[,'Lymphocytes']>1),grep('MMPART',colnames(CLR))]<-0 
#################### newcells
save(CLR, file=file.path(base.dir,'newcells',sprintf('PLSR-nonlymphocytes-%s',FILES[[1]])))


