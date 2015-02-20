library(flowCore)
library(flowStats)
library(flowViz)

f <- flowCore::read.flowSet(files=c('~/dunwich/FCS.Tony/FCS.repeats/CB00396E/day1/0U.fcs','~/dunwich/FCS.Tony/FCS.repeats/CB00396E/day2/0U.fcs'))

wf <- workFlow(f)
tl <- transformList(colnames(f)[c(9,12)], logicleTransform(), transformationId="lt")
add(wf, tl)

pars <- colnames(Data(wf[["base view"]]))[c(9,12)]
norm <- normalization(normFun=function(x, parameters, ...) warpSet(x, parameters, ...),
                      parameters=pars,
                      #arguments=list(grouping="GroupID", monwrd=TRUE),
                      normalizationId="Warping")

add(wf, norm)

