# --- HCV PACKAGE EXAMPLES --- #

# install package if necessary #

#devtools::install_github('bbroyle/evo3D')

library(evo3D)
library(tidyverse)

############# HEPC E1E2 ###############

# Fig 2a, b) E1 and E2 mapped to 8fsj ----

# define dataset
msa1 = 'HCV/e1hits_dropgap.fa'
msa2 = 'HCV/e2hits_dropgap.fa'
pdb1 = 'HCV/8fsj.cif'

# see wrapper defaults
show_evo3d_defaults()
show_evo3d_defaults('stat') # and specific module defaults

# run evo3d for single site entropy and block entropy #
res1 = run_evo3d(msa = list(msa1, msa2), detail = 2,
                 pdb = pdb1,
                 stat_controls = list(calc_site_entropy = TRUE, calc_block_entropy = TRUE,
                                      calc_avg_patch_entropy = TRUE),
                 pdb_controls = list(rsa_cutoff = NA, sasa_cutoff = 10,
                                     dist_cutoff = NA, max_patch = 15,
                                     distance_method = 'all')
)

# this dataframe can be inspected in R #
res1_df = res1$evo3d_df

# quick summary of stats #
summary(res1_df$site_entropy) # some invariant positions, some as high as 3.4 bits
summary(res1_df$mean_site_entropy) # some invariant pathes, some as high as 1.2 bits / residue
summary(res1_df$block_entropy) # some invariant patches, some as high as 7.7 bits

# median surface block entropy #
median(res1_df$block_entropy, na.rm = TRUE) # 3.7 bits

# save these results #
saveRDS(res1, 'HCV/e1e2_entropy_fig2ab.rds')

# write to b factor #
write_stat_to_pdb(res1, stat_name = 'site_entropy', outfile = 'HCV/e1e2_site_entropy_fig2a.pdb')
write_stat_to_pdb(res1, stat_name = 'block_entropy', outfile = 'HCV/e1e2_block_entropy_fig2b.pdb')

# Fig 2c) spatial null model testing ----

# read in results from fig2a,b
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')

# results can go straight to generate null models #
null_hap = generate_null_model(res1, n = 10000, len = 15, seed = 1219)

# calculate block entropy on these null haplotypes #
null_be = calculate_patch_stats(final_msa_subsets = null_hap$msa_subsets,
                                residue_df = null_hap$null_df,
                                stat_controls = list(calc_block_entropy = TRUE))

# save null block entropy distribution #
saveRDS(null_be, 'HCV/spatial_null_block_entropy.rds')

### ------- TEST FOR PATCH SIGNIFICNACE ---------  ###

# only track surface positions #
res1_df = res1$evo3d_df
surf = res1_df[!is.na(res1_df$block_entropy),]

# load in null block entropy results #
null_be = readRDS('HCV/spatial_null_block_entropy.rds')
null_vals = as.numeric(null_be$block_entropy) # 10K null values

# compare true block entropy to null dist #
true_vals = surf$block_entropy     # 276 observed patches

#  two-sided p value
m = length(null_vals)
p_emp = vapply(true_vals, function(x){
  plo = (sum(null_vals <= x, na.rm=TRUE)+1)/(m+1)
  phi = (sum(null_vals >= x, na.rm=TRUE)+1)/(m+1)
  min(1, 2*min(plo, phi))
}, 0.0)

surf$pval = p_emp

# Direction label
mu0 = mean(null_vals, na.rm=TRUE)
surf$direction = ifelse(surf$block_entropy < mu0, "conserved", "diverse")

# Standardized effect
sd0 = sd(null_vals, na.rm=TRUE)
surf$z = (surf$block_entropy - mu0)/sd0

# order table on abs(zscore)
hold = surf[order(-abs(surf$z)),]
nonover = filter_overlaps(hold, overlap = 0.33)

# 4 cons / 1 div -- 0.33 over
# 3 cons / 1 div -- 0 over (collapse dimer interface into one patch)

# multiplicity on pruned set
nonover$qval = p.adjust(nonover$pval, "BH")
nonover$bonf = p.adjust(nonover$pval, "bonferroni")

# quick counts
table(BH_5 = nonover$qval <= 0.05, Bonf_10 = nonover$bonf <= 0.10)
table(BH_5 = nonover$qval <= 0.05)

# final calls
calls = subset(nonover, qval <= 0.05)
calls[c("residue_id","block_entropy","direction","z","pval","qval","bonf")]

#    residue_id block_entropy direction         z       pval       qval       bonf
#494     685_E_     0.0000000 conserved -2.813062 0.00019998 0.00259974 0.00519948
#471     662_E_     0.2455437 conserved -2.639459 0.00019998 0.00259974 0.00519948
#285     476_E_     7.6539763   diverse  2.598397 0.00079992 0.00519948 0.02079792
#316     507_E_     0.3219429 conserved -2.585444 0.00039996 0.00346632 0.01039896
#415     606_E_     0.5867107 conserved -2.398250 0.00339966 0.01767823 0.08839116

# Fig 2d) linear block entropy, and comparison with spatial ----

# comparable linear sliding window (resolved region +/- 7 for window completeness) #
# evo3D does not currently have dedicated linear sliding window functionality #
# could add but is simply outside focus of the package #

# load in results from fig2 ab
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')
res1_df = res1$evo3d_df

# need to trim linear model to resolved range #
res1_df %>% group_by(msa) %>% filter(residue_id != '-') %>%
  summarise(mi = min(as.numeric(codon)),
            ma = max(as.numeric(codon)))

# first and last resolved pos of e1 in 8fsj (1-121)
# first and last resolved pos of e2 in 8fsj (38 - 321)

# need linear windows of 15 aa, with extended termini for maximum residue coverage #
e1 = res1_df %>% filter(msa == 'msa1', codon %in% 1:(121+7))
e2 = res1_df %>% filter(msa == 'msa2', codon %in% (38-7):(321 + 7))

##### --------- start with e1 ---------- ###
start_pos = 1
end_pos = 128
win_size = 15

# generate windows
windows = lapply(start_pos:(end_pos - win_size + 1), function(i) {
  seq(i, i + win_size - 1)
})

# collapse to "codon_patch" strings
codon_patches = vapply(windows, function(x) paste0(paste0(x, '_msa1'), collapse = "+"), character(1))

e1_patches = tibble(codon = 8:121, codon_patch = codon_patches) # 8 to 121 for central positions of window
e1_patches$msa_subset_id = e1_patches$codon
e1_patches$msa = 'msa1'

msa_set = .extract_msa_subsets(res1$msa_info_sets,
                               e1_patches)

e1_blent = sapply(msa_set, function(x) {
  # compute block entropy #
  block_entropy(x, translate = TRUE)
})

e1_patches$block_entropy = e1_blent

####--------- do e2 as well ------------ ###
start_pos = 31
end_pos = 328
win_size = 15

# generate windows
windows = lapply(start_pos:(end_pos - win_size + 1), function(i) {
  seq(i, i + win_size - 1)
})

# collapse to "codon_patch" strings
codon_patches = vapply(windows, function(x) paste0(paste0(x, '_msa2'), collapse = "+"), character(1))

e2_patches = tibble(codon = 38:321, codon_patch = codon_patches) # 38 to 321 for central positions of window
e2_patches$msa_subset_id = e2_patches$codon
e2_patches$msa = 'msa2'

msa_set = .extract_msa_subsets(res1$msa_info_sets,
                               e2_patches)

e2_blent = sapply(msa_set, function(x) {
  # compute block entropy #
  block_entropy(x, translate = TRUE)
})

e2_patches$block_entropy = e2_blent

# save these for other analysis
saveRDS(e1_patches, 'HCV/e1_linear_be.rds')
saveRDS(e2_patches, 'HCV/e2_linear_be.rds')

## --- calculate (spatial minus linear) block entropy --- ###

e1_patches = readRDS('HCV/e1_linear_be.rds')
e2_patches = readRDS('HCV/e2_linear_be.rds')
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')

hold = rbind(e1_patches, e2_patches)
hold$codon = as.character(hold$codon)
hold$linear_block_entropy = hold$block_entropy
hold$linear_patch = hold$codon_patch

res1_df = left_join(res1_df, hold %>% select(msa,codon,linear_block_entropy,linear_patch), by = c('msa', 'codon'))

res1_df$diff = res1_df$block_entropy - res1_df$linear_block_entropy
max(res1_df$diff, na.rm = TRUE) # 3.946
min(res1_df$diff, na.rm = TRUE) # -5.821

# we can add new statistics back to evo3d_results and write to pdb #
res1$evo3d_df = res1_df
write_stat_to_pdb(res1, stat_name = 'diff', outfile = 'HCV/e1e2_entropy_diff_fig2d.pdb')

# set up linear null model ----

# evo3D doesn't nativly support linear null model testing - yet we can
# approach either by rebuilding a evo3D_results list #
# or setup null model outside of generate_null_model() #

# testing separately #

# read in linear results
e1_patches = readRDS('HCV/e1_linear_be.rds')
e2_patches = readRDS('HCV/e2_linear_be.rds')

# read in res1 as scaffold #
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')

# leave MSA in place - replace evo3d_df #
# codon patches and msa column along with msa_info_sets build null model #
res1$evo3d_df = e1_patches

# results can go straight to generate null models #
e1_null = generate_null_model(res1, n = 10000, len = 15, seed = 1219, match_codon_frequency = FALSE)

# need to provide null haplotypes #
msa_set = e1_null$msa_subsets

# calculate block entropy on these null haplotypes #
# we will just use block_entropy() directly
e1_null_be = sapply(msa_set, function(x) {
  # compute block entropy #
  block_entropy(x, translate = TRUE)
})

# save null block entropy distribution #
saveRDS(e1_null_be, 'HCV/e1_null_block_entropy.rds')

# ---- lets generate e2 null model as well --- ###
res1$evo3d_df = e2_patches
e2_null = generate_null_model(res1, n = 10000, len = 15, seed = 1219)
msa_set = e2_null$msa_subsets
e2_null_be = sapply(msa_set, function(x) {block_entropy(x, translate = TRUE)})
saveRDS(e2_null_be, 'HCV/e2_null_block_entropy.rds')

# *** test linear empirical block entropy against null models ----

# load res1 to get residue id #
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')
res1_df = res1$evo3d_df

# do one at a time #
e1_patches = readRDS('HCV/e1_linear_be.rds')
e1_patches$codon = as.character(e1_patches$codon)
e1_patches = left_join(e1_patches, res1_df %>% select(codon, msa, residue_id))

# load in null block entropy results #
null_vals = readRDS('HCV/e1_null_block_entropy.rds')

# compare true block entropy to null dist #
true_vals = e1_patches$block_entropy

#  two-sided p value
m = length(null_vals)
p_emp = vapply(true_vals, function(x){
  plo = (sum(null_vals <= x, na.rm=TRUE)+1)/(m+1)
  phi = (sum(null_vals >= x, na.rm=TRUE)+1)/(m+1)
  min(1, 2*min(plo, phi))
}, 0.0)

e1_patches$pval = p_emp

# Direction label
mu0 = mean(null_vals, na.rm=TRUE)
e1_patches$direction = ifelse(e1_patches$block_entropy < mu0, "conserved", "diverse")

# Standardized effect
sd0 = sd(null_vals, na.rm=TRUE)
e1_patches$z = (e1_patches$block_entropy - mu0)/sd0

# order table on abs(zscore)
hold = e1_patches[order(-abs(e1_patches$z)),]
nonover = filter_overlaps(hold, overlap = 0.33)

# multiplicity on pruned set
nonover$qval = p.adjust(nonover$pval, "BH")
nonover$bonf = p.adjust(nonover$pval, "bonferroni")

table(BH_5 = nonover$qval <= 0.05, Bonf_10 = nonover$bonf <= 0.10)
table(BH_5 = nonover$qval <= 0.05)

# no significant E1 regions under strict threshold
# none under lower significance either

# ------ moving to E2 --- ###
# do one at a time #
e2_patches = readRDS('HCV/e2_linear_be.rds')
e2_patches$codon = as.character(e2_patches$codon)
e2_patches = left_join(e1_patches, res1_df %>% select(codon, msa, residue_id))

# load in null block entropy results #
null_vals = readRDS('HCV/e2_null_block_entropy.rds')

# compare true block entropy to null dist #
true_vals = e2_patches$block_entropy

#  two-sided p value
m = length(null_vals)
p_emp = vapply(true_vals, function(x){
  plo = (sum(null_vals <= x, na.rm=TRUE)+1)/(m+1)
  phi = (sum(null_vals >= x, na.rm=TRUE)+1)/(m+1)
  min(1, 2*min(plo, phi))
}, 0.0)

e2_patches$pval = p_emp

# Direction label
mu0 = mean(null_vals, na.rm=TRUE)
e2_patches$direction = ifelse(e2_patches$block_entropy < mu0, "conserved", "diverse")

# Standardized effect
sd0 = sd(null_vals, na.rm=TRUE)
e2_patches$z = (e2_patches$block_entropy - mu0)/sd0

# order table on abs(zscore)
hold = e2_patches[order(-abs(e2_patches$z)),]
nonover <- filter_overlaps(hold, overlap = 0.33)

# multiplicity on pruned set
nonover$qval = p.adjust(nonover$pval, "BH")
nonover$bonf = p.adjust(nonover$pval, "bonferroni")

#
table(BH_5 = nonover$qval <= 0.05, Bonf_10 = nonover$bonf <= 0.10)
table(BH_5 = nonover$qval <= 0.05)

# final calls
calls = subset(nonover, qval <= 0.05)
calls[c("block_entropy","direction","z","pval","qval","bonf")]

# examine spatial versus linear codon patches ----
e1_patches = readRDS('HCV/e1_linear_be.rds')
e2_patches = readRDS('HCV/e2_linear_be.rds')
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')

# combine linear data and add residue_id and block_entropy
full_set = rbind(e1_patches, e2_patches)
full_set$linear_patch = full_set$codon_patch
full_set$linear_entropy = full_set$block_entropy
full_set$codon_patch = NULL
full_set$block_entropy = NULL
full_set$codon = as.character(full_set$codon)

full_set = left_join(full_set, res1$evo3d_df %>% select(msa, codon, residue_id, codon_patch, block_entropy))

keep = full_set %>% filter(!is.na(linear_entropy), !is.na(block_entropy))
ggplot(keep, aes(block_entropy, linear_entropy))+geom_point()
cor(keep$linear_entropy, keep$block_entropy) # 0.6587

keep$diff = keep$block_entropy - keep$linear_entropy
min(keep$diff) # range -5.8
max(keep$diff) # range 3.9

# find patch overlap
keep$overlap = NA
for(i in 1:nrow(keep)){
  p1 = keep$codon_patch[i]
  p2 = keep$linear_patch[i]

  p1 = strsplit(p1, '\\+')[[1]]
  p2 = strsplit(p2, '\\+')[[1]]

  overlap = length(intersect(p1, p2))

  keep$overlap[i] = overlap
}

mean(keep$overlap) # 7.5 residues
sd(keep$overlap)   # 2.1 residues
median(keep$overlap) # 8 residues
min(keep$overlap)   # 2
max(keep$overlap)   # 13


# Extra analyses - pulling out true haplotype counts ----
res1 = readRDS('HCV/e1e2_entropy_fig2ab.rds')

# haplotypes per window
set = res1$final_msa_subsets$`476_E_` # 476_E_
set = apply(set, 1, paste0, collapse = '')

# loop through and convert to aa #
hi = lapply(set, function(x){
  aa = seqinr::translate(seqinr::s2c(x))
  aa = paste0(aa, collapse = '')
}) %>% unlist()

length(hi) # 271 seqs
table(hi)
length(unique(hi)) # 224 sequences
max(table(hi)) # with most frequent occupying 9/271 only 3.3%

# haplotypes per window
set = res1$final_msa_subsets$`606_E__pdb1` # 606_E_
set = apply(set, 1, paste0, collapse = '')

# loop through and convert to aa #
hi = lapply(set, function(x){
  aa = seqinr::translate(seqinr::s2c(x))
  aa = paste0(aa, collapse = '')
}) %>% unlist()

length(hi) # 271 seqs
table(hi) # 251 out of 271
251 / 271 # 92.6%

set = res2$final_msa_subsets$`561_E__pdb1` # 561_E_
set = apply(set, 1, paste0, collapse = '')

# loop through and convert to aa #
hi = lapply(set, function(x){
  aa = seqinr::translate(seqinr::s2c(x))
  aa = paste0(aa, collapse = '')
}) %>% unlist()


length(hi) # 271 seqs
table(hi) # 247 out of 271
247 / 271  # 91.1%


