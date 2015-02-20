CELL.TYPES <- c("Memory Eff", "Memory Treg", "Naive Eff", "Naive Treg")

counts <- data.frame()
for (f in list.files('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR', pattern='.*.RData', full.names=TRUE)) {
    print(f)
    print(load(f))
    counts <- rbind(counts,data.frame(name=gsub('.RData','',basename(f)), total=nrow(CLR), t(apply(CLR,2,function(x)length(which(x==1))))))
}
save(counts, file='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/counts.RData')

print(load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/CellPhenotypes/magnetic-manual-gates2/counts.RData'))

counts$Lymphocytes/counts$total



