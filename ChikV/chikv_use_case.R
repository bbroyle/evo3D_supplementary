# --- HCV PACKAGE EXAMPLES --- #

# install package if necessary #

#devtools::install_github('bbroyle/evo3D')

library(evo3D)
library(tidyverse)

############# # # # # # ChikV E1E2 # # # # # ###############

# fig3 (ab) residue wise block entropy ----

# set dataset
msa1 = 'ChikV/chikv_e1hits.fa'
msa2 = 'ChikV/chikv_e2hits.fa'
pdb1 = 'ChikV/8fcg.pdb'

# same as HCV, but this time specify residue patches
# HCV residue and codon patch_mode is equivalent b/c
# there are no codon duplicates
pc = list(rsa_cutoff = NA, sasa_cutoff = 10,
          dist_cutoff = NA, max_patch = 15,
          distance_method = 'all', patch_mode = 'residue')

res2 = run_evo3d(list(msa1, msa2), pdb1, verbose = 2,
                 pdb_controls = pc,
                 stat_controls = list(calc_block_entropy = TRUE),
                 analysis_mode = 'residue') # also specifiy residue analysis mode

# count codon per window #
res2_df = res2$evo3d_df
table(res2_df$codon_len)
table(res2_df$unique_codon)  # as few as 10 unique codons in 15 positions #

# some codon duplicates captured #
table(res2_df$unique_codon == 15)
30/(30 + 2543) # 1.2%

# unify to codon - average #
cods = unique(res2_df$codon_id)
cods = cods[!is.na(cods) & cods != '-']

res2_df$res_n = NA
res2_df$avg_be = NA
res2_df$min_be = NA
res2_df$max_be = NA

for(c in cods){
  ro = which(res2_df$codon_id == c)

  scores = res2_df$block_entropy[ro]

  res2_df$res_n[ro] = length(scores[!is.na(scores)])
  res2_df$avg_be[ro] = mean(scores, na.rm = TRUE)
  res2_df$min_be[ro] = min(scores, na.rm = TRUE)
  res2_df$max_be[ro] = max(scores, na.rm = TRUE)
}

# add new stats back to results list #
res2$evo3d_df = res2_df

# THESE ARE FIGURE 3 C and D #
write_stat_to_pdb(res2, stat_name = 'avg_be', outfile = 'ChikV/chikv_avg_block_ent.pdb')
saveRDS(res2, 'ChikV/chikungunya_results.rds')

# -- fig 3b ----
res2 = readRDS('ChikV/chikungunya_results.rds')
res2_df = res2$evo3d_df

# how are results on e2 146 #
look = res2_df %>% filter(grepl('146_._', residue_id), msa == 'msa2')
look
# E - 0.40
# F - 0.38
# G - 0.40
# H - 1.00

# summarize at codon level
coddf = res2_df %>% select(codon_id, codon, msa, res_n, avg_be, min_be, max_be) %>% distinct()

# drop codons with no spatial windows (buried or outside PDB coverage)
coddf = coddf %>% filter(res_n != 0)

# difference between max and min values
coddf$diff = coddf$max_be - coddf$min_be
summary(coddf$diff) # median 0.03 IQR 0 to 0.12, max 1.7

table(coddf$diff >= 1) # 14 at least 1

# summary of average block entropy
summary(coddf$avg_be)  # median 0.4, min 0, max 2.2

# view by gene
e1 = coddf %>% filter(msa == 'msa1')
summary(e1$avg_be) # median 0.2

e2 = coddf %>% filter(msa == 'msa2')
summary(e2$avg_be) # median 0.6

# top diversity positions
coddf %>% arrange(-avg_be) %>% head()

# codon 182 msa2 (residue 119 E2)
res2_df %>% filter(codon == 182, msa == 'msa2')

# need site entropy of each codon as well #
coddf = calculate_site_entropy(msa_info_sets = res2$msa_info_sets,
                               residue_df = coddf)

# view those top diversity patches again #
coddf %>% arrange(-avg_be) %>% head()

write_tsv(coddf, 'ChikV/chikungunya_avg_be_per_codon.tsv')

# figS2: check interface of e1e2-mxra8 ----
msa1 = 'ChikV/chikv_e1hits.fa'
msa2 = 'ChikV/chikv_e2hits.fa'
pdb = 'ChikV/6jo8.cif'
pp2 = list(interface_dist_cutoff = 4, distance_method = 'ca', rsa_cutoff = NA, sasa_cutoff = 10, max_patch = 35, dist_cutoff = NA)


res_int = run_evo3d(list(msa1, msa2),
                    pdb = pdb,
                    analysis_mode = 'residue',
                    interface_chain = c('M', 'N', 'O'),
                    stat_controls = list(calc_site_entropy = TRUE),
                    pdb_controls = pp2, verbose = 2, detail = 2)

saveRDS(res_int, 'ChikV/chikungunya_interface_results.rds')


# ----------- save previous block entropy with this structure ------------ #

# three copies of interface -- what were there site entropies #
df_int = res_int$evo3d_df
int = tail(df_int, 3)

p1 = int$codon_patch[1]
p1 = strsplit(p1, '\\+')[[1]]
ents1 = tibble(codon_id = p1, int = 'E1/E2 - MXRA8: chain M')

p2 = int$codon_patch[2]
p2 = strsplit(p2, '\\+')[[1]]
ents2 = tibble(codon_id = p2, int = 'E1/E2 - MXRA8: chain N')


p3 = int$codon_patch[3]
p3 = strsplit(p3, '\\+')[[1]]
ents3 = tibble(codon_id = p3, int = 'E1/E2 - MXRA8: chain O')

ents = rbind(ents1, ents2, ents3)

ents = left_join(ents, df_int %>% select(codon_id, site_entropy) %>% distinct())

hold = tibble(site_entropy = coddf$site_entropy)
hold$int = 'Full E1/E2 Surface (8fcg)'

ents$codon_id = NULL

ents = rbind(ents, hold)

# ADD SITE ENTROPY FROM 8fcg full surface #

p = ggplot(ents, aes(site_entropy))+
  geom_histogram()+
  facet_wrap(~int, ncol = 1, scales = 'free_y')+
  theme_bw()+
  xlab('Shannon Entropy per Site')+
  ylab('Count')

ggsave("ChikV/mxra8_site_ent.svg", p, device = "svg", width = 3, height = 5)

# add is in interface to full codon table #
coddf$in_int = ifelse(coddf$codon_id %in% c(p1, p2, p3), T, F)

summary(coddf %>% filter(in_int) %>% select(site_entropy)) # median 0, IQR 0 to 0.03
#
summary(coddf %>% filter(!in_int) %>% select(site_entropy)) # median 0, IQR 0 to 0.0

table(coddf$in_int, coddf$site_entropy > 1)
coddf %>% filter(site_entropy > 1)
res2_df %>% filter(codon == 211, msa == 'msa1') # highest entropy single site - residue 211 E1
res2_df %>% filter(codon == 137, msa == 'msa2') # highest entorpy site within interface - residue 74 E2
