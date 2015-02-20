#d = read.csv('/chiswick/data/hui/hui_archive/flow/FACS/Tony_ls_dgap_pairs/LS_DGAP_paired_data_CW_copy_tony_csv', header=T, as.is=T)
d = read.csv('~/dunwich/Projects/IL2/LS_DGAP_paired_data_CW_copy_tony_csv', header=T, as.is=T)

d = d[1:51, -c(5,6)]
d$HC[7] = 'pair 40'
colnames(d) = c('pairs', 'naive', 'memory', 'actv', 'naive.1', 'memory.1', 'actv.1')



## only for the paired data, do paired t-test
d.p = na.omit(d)
d.p$dif.naive = with(d.p, naive.1-naive)
d.p$dif.memory = with(d.p, memory.1 - memory)
d.p$dif.actv = with(d.p, actv.1 - actv) 

plot(d.p$dif.naive)
plot(d.p$dif.memory)
plot(d.p$dif.actv)

t.test(d.p$dif.naive, alternative = "two.sided")
# t = 1.6879, df = 39, p-value = 0.09941

t.test(d.p$dif.memory, alternative = "two.sided")
# t = -0.0077, df = 39, p-value = 0.9939

t.test(d.p$dif.actv, alternative = "two.sided")
# t = -0.7752, df = 39, p-value = 0.4429


## assign Row i with empty pair info with the pair info of Row i-1
## for each trio group, if phenotype measured from 2 IDs, use the average to form 1 pair
  
n = dim(d)[1]
for (i in 1:n) {
  d$pairs[i] = ifelse(d$pairs[i] != '', d$pairs[i], d$pairs[i-1])
}  

d.dup = d[c(which(duplicated(d$pairs)), which(duplicated(d$pairs))-1), ]

d.dup.split = split(d.dup, d.dup$pairs)

d.dup.avrg = vector('list', length(d.dup.split))
for (j in names(d.dup.split)) {
  d.dup.avrg[[j]]$pairs = names(d.dup.split)[j]

  d.dup.avrg[[j]]$naive = mean(na.omit(d.dup.split[[j]]$naive))
  d.dup.avrg[[j]]$memory = mean(na.omit(d.dup.split[[j]]$memory)) 
  d.dup.avrg[[j]]$actv = mean(na.omit(d.dup.split[[j]]$actv)) 

  d.dup.avrg[[j]]$naive.1 = mean(na.omit(d.dup.split[[j]]$naive.1))
  d.dup.avrg[[j]]$memory.1 = mean(na.omit(d.dup.split[[j]]$memory.1)) 
  d.dup.avrg[[j]]$actv.1 = mean(na.omit(d.dup.split[[j]]$actv.1)) 
}


## combine with orignal paired data, do paired t-test  
d.dup.avrg = do.call('rbind.data.frame', d.dup.avrg)
d.dup.avrg$pairs = rownames(d.dup.avrg)
d.dup.avrg$dif.naive = with(d.dup.avrg, naive.1-naive)
d.dup.avrg$dif.memory = with(d.dup.avrg, memory.1 - memory)
d.dup.avrg$dif.actv = with(d.dup.avrg, actv.1 - actv) 

d.p.dup.avrg = rbind(subset(d.p, !pairs %in% d.dup.avrg$pairs), d.dup.avrg)

t.test(d.p.dup.avrg$dif.naive, alternative = "two.sided")
# t = 2.197, df = 39, p-value = 0.03403

t.test(d.p.dup.avrg$dif.memory, alternative = "two.sided")
# t = 0.4235, df = 39, p-value = 0.6743

t.test(d.p.dup.avrg$dif.actv, alternative = "two.sided")
# t = -0.5241, df = 39, p-value = 0.6032





## for each trio group, use data of 1 ID twice to form two pairs

d.dup.pair = vector('list', length(d.dup.split))
for (j in names(d.dup.split)) {
  d.dup.pair[[j]] = d.dup.split[[j]]
  if (is.na(d.dup.split[[j]]$naive[1])) {
    d.dup.pair[[j]][1, 2:4] = d.dup.pair[[j]][2, 2:4]
  }
  if (is.na(d.dup.split[[j]]$naive.1[1])) {
    d.dup.pair[[j]][1, 5:7] = d.dup.pair[[j]][2, 5:7]
  }
}

d.dup.pair = do.call('rbind', d.dup.pair)
d.dup.pair$dif.naive = with(d.dup.pair, naive.1-naive)
d.dup.pair$dif.memory = with(d.dup.pair, memory.1 - memory)
d.dup.pair$dif.actv = with(d.dup.pair, actv.1 - actv)

d.dup.pair.comb = rbind(subset(d.p, !pairs %in% d.dup.avrg$pairs), d.dup.pair) 
d.dup.pair.comb$pairs = as.numeric(gsub('pair ', '', d.dup.pair.comb$pairs))

 



## bootstapping in 2 strata (pairs, trios), use weight 3/4 for trios

library(boot)
set.seed(10000)

d.dup.pair.comb$group = c(rep(1, 29), rep(2, 22))

d.pair = data.frame(d.dup.pair.comb[1:29, -c(2:7,11)], dif.naive.1 = NA, dif.memory.1 = NA, dif.actv.1 = NA, group = d.dup.pair.comb[1:29, 11])

d.trio = d.dup.pair.comb[-(1:29), ]
trio.for.boot = cbind(d.trio[seq(1,dim(d.trio)[1],by=2), c(1,8:10)],d.trio[seq(2,dim(d.trio)[1],by=2), 8:11])
colnames(trio.for.boot)[5:7] = colnames(d.pair)[5:7]

d.boot = rbind(d.pair, trio.for.boot)

w = c(rep(1, 29), rep(3/4, 11))/sum(c(rep(1, 29), rep(3/4, 11)))
grp1 = 1:table(d.boot$group)[1]


avrg = function(d, weight=w) {
     mn.n <- sum(d$dif.naive * weight, d$dif.naive.1[-grp1] * weight[-grp1])  
     mn.m <- sum(d$dif.memory * weight, d$dif.memory.1[-grp1] * weight[-grp1]) 
     mn.a <- sum(d$dif.actv * weight, d$dif.actv.1[-grp1] * weight[-grp1])
     #c(mn.n, sd.n, mn.m)
     c(mn.n, mn.a, mn.m)
}

bstp = boot(data=d.boot, statistic=avrg, R=100000, weights = w, stype = "w", strata = d.boot$group)
m = bstp$t0
s = sqrt(diag(var(bstp$t)))

# naive
t.n = m[1]/s[1]
pval.n = 2*pt(abs(t.n), 50, lower.tail=FALSE)
# t = 1.586395, df = 50, p-value = 0.1189545
 
# memory
t.m = m[2]/s[2]
pval.m = 2*pt(abs(t.m), 50, lower.tail=FALSE)
# t = 1.577787, df = 50, p-value = 0.1209212

# activated
t.a = m[3]/s[3]
pval.a = 2*pt(abs(t.a), 50, lower.tail=FALSE)
# t = -0.9558598, df = 50, p-value = 0.343742
