Data:
/home/hui/flow/FACS/CD4CTLA-4.from.Laura/data/LS_DGAP_paired_data_CW_copy_tony_csv
-- Pair 5 include 2 pairs, so name them as 'pair 5' and 'pair 40'
-- 29 pairs and 11 trios, data collected on the same day
 


R code:
/home/hui/flow/FACS/scripts/ip/tony.ls.dgap.paired.treg.r



Statistical analysis:

To test if there is a difference, for each of the three phenotypes (naive/memory/activated Tregs), in matched pairs (healthy controls vs T1D)

Three approaches:

1. for each of the trios, remove the 3rd ID, paired t-test (40 pairs)

-- naive:  
   t = 1.6879, df = 39, p-value = 0.09941
-- memory
   t = -0.0077, df = 39, p-value = 0.9939
-- activated
   t = -0.7752, df = 39, p-value = 0.4429



2. for each of the trios, if phenotype measured from 2 IDs, use the average to form 1 pair, paired t-test (40 pairs)

-- naive:  
   t = 2.197, df = 39, p-value = 0.03403
-- memory
   t = 0.4235, df = 39, p-value = 0.6743
-- activated
   t = -0.5241, df = 39, p-value = 0.6032



3. for each of the trios, use 1 ID twice to form 2 pairs, weighted (1 for original pairs, 3/4 for pairs formed from trios) and stratified (pairs, trios) bootstrapping and paired t-test (51 pairs)

data format for bootstrap:
                diff.pair1  diff.pair2
Group.pairs     *           NA
.               *           NA
.               *           NA
Group.trios     *           *
.               *           *
.               *           *


-- naive:  
   t = 1.586, df = 50, p-value = 0.119
-- memory
   t = 1.578, df = 50, p-value = 0.121
-- activated
   t = -0.956, df = 50, p-value = 0.344
