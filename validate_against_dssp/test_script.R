# This script is to validate evo3D DSSP SASA logic implementation  #
# against MKDSSP ---- MKDSSP 4.4.10 was installed and ran locally  #
# but would need to be installed locally on your end to rerun that #
# part of the analysis .

library(evo3D)
library(tidyverse)

# functions ----

# function to load mkdssp outputs #
parse_dssp = function(x){
  # look for header line # *RESIDUE *AA
  ro = grep('# *RESIDUE *AA', x)
  headers = strsplit(x[ro], ' +')[[1]][-1]

  headers[1] = 'id'

  tab = x[(ro+1):(length(x))]

  id = as.numeric(trimws(str_sub(tab, 1, 5)))
  res = as.numeric(trimws(str_sub(tab, 6, 10)))
  chain = trimws(str_sub(tab, 11, 12))
  aa = trimws(str_sub(tab, 13, 14))
  str = str_sub(tab, 15, 26)
  bp1 = trimws(str_sub(tab, 27, 30))
  bp2 = trimws(str_sub(tab, 31, 34))
  acc = as.numeric(trimws(str_sub(tab, 35, 38)))
  nh_o1 = trimws(str_sub(tab, 39, 52))
  o_nh1 = trimws(str_sub(tab, 53, 62))
  nh_o2 = trimws(str_sub(tab, 63, 74))
  o_nh2 = trimws(str_sub(tab, 75, 85))
  tco = as.numeric(trimws(str_sub(tab, 86, 91)))
  kap = as.numeric(trimws(str_sub(tab, 92, 97)))
  alp = as.numeric(trimws(str_sub(tab, 98, 103)))
  phi = as.numeric(trimws(str_sub(tab, 104, 109)))
  psi = as.numeric(trimws(str_sub(tab, 110, 115)))
  xcor = as.numeric(trimws(str_sub(tab, 116, 122)))
  ycor = as.numeric(trimws(str_sub(tab, 123, 129)))
  zcor = as.numeric(trimws(str_sub(tab, 130, 136)))

  tab = data.frame(id = id, residue = res, chain = chain, aa = aa, structure = str, bp1 = bp1, bp2 = bp2,
                   acc = acc, nh_o1 = nh_o1, o_nh1 = o_nh1, nh_o2 = nh_o2, o_nh2 = o_nh2,
                   tco = tco, kappa = kap, alpha = alp, phi = phi, psi = psi, x = xcor, y = ycor, z = zcor)

  return(tab)

}

# gathered sample set of 95 pdb structures ----
pdb_set = read_lines('validate_against_dssp/sample_95_structures.txt')

# download data
#bio3d::get.pdb(pdb_set, path = 'validate_against_dssp/test_cif', format = 'cif')
list.files('validate_against_dssp/test_cif/') %>% length()  # all downloaded

# grab information about these structures #

#pdb_info = bio3d::pdb.annotate(pdb_set)
#write_tsv(pdb_info, 'validate_against_dssp/pdb_95_annotation.tsv')
pdb_info = read_tsv('validate_against_dssp/pdb_95_annotation.tsv')

# Process PDB INFO
chain_num = pdb_info %>% group_by(structureId) %>% summarise(n = n())
table(chain_num$n)

# 1 to 23 chains #

# what kinds of structural features are present #
chain_num$structureId = tolower(chain_num$structureId)

chain_num$pdb_h1 = NA
chain_num$pdb_h2 = NA
chain_num$pdb_het = NA
chain_num$pdb_alt = NA
chain_num$pdb_ins = NA
chain_num$pdb_pca = NA   # i know 6 structures had differing number of residue positions from dssp and .calculate_accessibility() -
                         # inspecting them shows PCA in seqres but labelled as hetatm, so dropped in evo3D but included in DSSP -

for(i in 1:nrow(chain_num)){
  id = paste0(chain_num$structureId[i], '.cif')
  print(i)
  pdb = bio3d::read.cif(paste0('validate_against_dssp/test_cif/', id), rm.alt = F)

  x = table(is.na(pdb$atom$alt))
  if(!is.na(x['FALSE'])){
    chain_num$pdb_alt[i] = x['FALSE']
  }

  x = table(is.na(pdb$atom$insert))
  if(!is.na(x['FALSE'])){
    chain_num$pdb_ins[i] = x['FALSE']
  }

  x = table(pdb$atom$type)
  if(!is.na(x['HETATM'])){
    chain_num$pdb_het[i] = x['HETATM']
  }

  x = table(pdb$atom$elesy[pdb$atom$type == 'ATOM'])
  if(!is.na(x['H'])){
    chain_num$pdb_h1[i] = x['H']
  }

  x = table(pdb$atom$elety[pdb$atom$type == 'ATOM'])
  if(any(grepl('^H', names(x)))){
    chain_num$pdb_h2[i] = length(grep('^H', pdb$atom$elety[pdb$atom$type == 'ATOM']))
  }

  x = table(pdb$atom$resid)
  if(!is.na(x['PCA'])){
    chain_num$pdb_pca[i] = x['PCA']
  }
}


# look at spread
table(is.na(chain_num$pdb_ins)) # 55 have insert codes
table(is.na(chain_num$pdb_alt)) # 35 have alt atoms
table(is.na(chain_num$pdb_het)) # 80 have hetatm
table(is.na(chain_num$pdb_h1))  # 13 have hydrogen in elesy
table(is.na(chain_num$pdb_h2))  # 13 have hydrogen in elety
table(is.na(chain_num$pdb_pca))  # 6 have hydrogen in PCA

# Compare the two results ----
f3 = list.files('validate_against_dssp/dssp_res_cif/')
f3 = gsub('.dssp', '', f3)
chain_num$dssp2 = NA
chain_num$dssp2 = f3[match(chain_num$structureId, tolower(f3))]

# takes a while to load all 95 cif files #
# can parallelize
library(parallel)

cl = makeCluster(8)

# load packages
clusterEvalQ(cl, {
  library(tidyverse)
  library(bio3d)
  library(evo3D)
})

# export dataframe and helper function
clusterExport(cl, c("chain_num", "parse_dssp"), envir = environment())

res = parLapply(cl, seq_len(nrow(chain_num)), function(i) {
  tryCatch({
    id1 = paste0("validate_against_dssp/test_cif/", chain_num$structureId[i], ".cif")
    id2 = paste0("validate_against_dssp/dssp_res_cif/", chain_num$dssp2[i], ".dssp")

    # MKDSSP results
    x = read_lines(id2)
    x = parse_dssp(x) %>% select(residue, chain, acc)

    # ---- PATH 1 (altLoc kept and non-canonical AA fixed -- not evo3D standard) ----
    pdb = bio3d::read.cif(id1, rm.alt = FALSE)

    # fix PCA
    ro = pdb$atom$resid == "PCA"
    if (any(ro)) pdb$atom$type[ro] = "ATOM"

    # compare with evo3D
    sasa1 = evo3D::.calculate_accessibility(pdb)
    sasa1$chain2 = paste0(sasa1$orig_insert, sasa1$orig_chain)

    # join to DSSP res so we dont miss any original DSSP results
    y = left_join(x, sasa1, by = c("residue" = "orig_resno", "chain" = "chain2"))
    y$diff = abs(y$sasa - y$acc)
    err2 = max(y$diff, na.rm = TRUE)

    # ---- PATH 2 (altLoc removed; "no_edit_sasa") ----
    pdb2 = bio3d::read.cif(id1, rm.alt = TRUE)

    sasa2 = evo3D::.calculate_accessibility(pdb2)
    sasa2$chain2 = paste0(sasa2$orig_insert, sasa2$orig_chain)
    colnames(sasa2)[6] <- "no_edit_sasa"

    # join so we can compare across 3 methods
    y = left_join(y, sasa2 %>% select(residue_id, no_edit_sasa), by = "residue_id")
    y$diff2 = abs(y$no_edit_sasa - y$acc)
    err3 = max(y$diff2, na.rm = TRUE)

    y$pdb = chain_num$structureId[i]
    keep = y[, c("pdb", "sasa", "acc", "no_edit_sasa", 'residue', 'chain', 'aa')]

    list(ok = TRUE, i = i, err2 = err2, err3 = err3, keep = keep)
  }, error = function(e) {
    list(ok = FALSE, i = i, err2 = NA_real_, err3 = NA_real_,
         keep = NULL, error = conditionMessage(e))
  })
})

stopCluster(cl)

# any missed? #
ok  = vapply(res, function(x) x$ok, logical(1))
table(ok) # all processed

# write results back into chain_num (per pdb max error)
chain_num$err2 = vapply(res, `[[`, numeric(1), "err2")
chain_num$err3 = vapply(res, `[[`, numeric(1), "err3")

# this is the MAXIMUM difference across a whole PDB -- residue specific histogram and IQR is better #
hist(chain_num$err2)
summary(chain_num$err2) # bounded between -0.5 and 0.5 (this is with hetatm and altloc "corrections")

hist(chain_num$err3)
summary(chain_num$err3) # Max error all the way to 88.356 (could be AltLOC or hetatm) - median error 0.5

# collect the residue table
keep_res2 = dplyr::bind_rows(lapply(res, `[[`, "keep"))

# if residue is NA - that is from MKDSSP and signifies chain breaks #
keep_res2 = keep_res2 %>% filter(!is.na(residue))

# write these two chain_num and keep_res2 results #
write_tsv(keep_res2, 'validate_against_dssp/residue_sasa_vs_dssp.tsv')
write_tsv(chain_num, 'validate_against_dssp/pdb_issue_annotations.tsv')

# inspect results ----
keep_res2 = read_tsv('validate_against_dssp/residue_sasa_vs_dssp.tsv')
chain_num = read_tsv('validate_against_dssp/pdb_issue_annotations.tsv')

# As described in methods of paper
keep_res2$diff1 = keep_res2$sasa - keep_res2$acc
hist(keep_res2$diff1) # differences bounded -0.5 and 0.5 from rounding values #
summary(keep_res2$diff1) #(median 0)#
cor(keep_res2$sasa, keep_res2$acc) # 0.9999784 across 114140 residues

# compare to native evo3D implement
# important to note - most disagreement comes from altLoc atoms, and
# its reasonable to argue using one atom location versus multiple is
# a better approach for sasa calc #
keep_res2$diff2 = keep_res2$no_edit_sasa - keep_res2$acc
hist(keep_res2$diff2) # differences extend to -88 and +34 # largely from AltLOC - some from non-canonical hetatm #
summary(keep_res2$diff2) # 7 NA (times PCA was present) - median 0
hold = keep_res2 %>% filter(!is.na(no_edit_sasa))
cor(hold$no_edit_sasa, hold$acc) # 0.9997519 across 114113 residues

quantile(abs(keep_res2$diff2), probs = seq(0.99,1,0.001), na.rm = TRUE)  # 99.4% 0.5 or less diff #


# groub by pdb issue - alt loc and noncanonical AA #
chain_num$pdb_alt[is.na(chain_num$pdb_alt)] = 0
noalt = chain_num %>% filter(pdb_alt == 0) %>% pull(structureId)
keep_res2$has_altLoc = ifelse(keep_res2$pdb %in% noalt, FALSE, TRUE)

# add if pca is not 0
chain_num$pdb_pca[is.na(chain_num$pdb_pca)] = 0
nopca = chain_num %>% filter(pdb_pca == 0) %>% pull(structureId)
keep_res2$haspca = ifelse(keep_res2$pdb %in% nopca, FALSE, TRUE)

keep_res2$label = NA

ro1 = which(keep_res2$has_altLoc)
ro2 = which(keep_res2$haspca)

ro3 = intersect(ro1, ro2) # both problems
ro4 = setdiff(ro1, ro2)   # only altloc
ro5 = setdiff(ro2, ro1)   # only noncanonical

keep_res2$label = 'no_issue'
keep_res2$label[ro3] = 'altLoc_noncanonical_AA'
keep_res2$label[ro4] =  'altLoc'
keep_res2$label[ro5] = 'noncanonical_AA'

# order by mean(abs(diff2)) -- or just abs(max(diff2)) #
hold = keep_res2 %>% group_by(pdb) %>%
  summarise(nres = n(),
            o1 = mean(abs(diff2), na.rm = TRUE),
            o2 = abs(max(diff2, na.rm = TRUE)),
            o3 = sum(abs(diff2) > 5, na.rm = TRUE), # adding a third number of residues missed by 5 Angstrom +
            o4 = o3/nres) # or o3 normalized by number of residues

# both essentially highlight the same issue -- (altLoc is more common - non canocial maybe has more effect (we need to check amount of impacted residues))
keep_res2$pdb = factor(keep_res2$pdb, levels = hold %>% arrange(o1) %>% pull(pdb))
ggplot(keep_res2, aes(pdb, diff2, color = label))+geom_boxplot()+coord_flip()

keep_res2$pdb = factor(keep_res2$pdb, levels = hold %>% arrange(o2) %>% pull(pdb))
ggplot(keep_res2, aes(pdb, diff2, color = label))+geom_boxplot()+coord_flip()

# sorting by percent of residues with greater than 5 ang difference
keep_res2$pdb = factor(keep_res2$pdb, levels = hold %>% arrange(o4) %>% pull(pdb))
ggplot(keep_res2, aes(pdb, diff2, fill = label, color = label))+
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = -5, ymax = 5, fill = 'grey82')+
  geom_boxplot()+coord_flip()+
  theme_bw()+
  ylab('Default evo3D SASA minus MKDSSP4.4.10')+
  scale_y_continuous(breaks = seq(-90, 40, 10))+
  xlab('PDB id sorted by percent of residues with larger than 5 Angstrom^2 difference') +
  theme(
    legend.position = c(0.09, 0.02),
    legend.justification = c("left", "bottom")
  )

# add label type to hold
hold$label = keep_res2$label[match(hold$pdb, keep_res2$pdb)]
hold %>% filter(label == 'noncanonical_AA') # only 3 out of 817 residues affected when noncanonical was the only issue

hold %>% arrange(-o4) # 7d9z as high as 10% of residues are off - yet only has 1 PCA atom -- same with 4toy (5% off)
# take away is dssp uses altloc but I disagree with its use for spatial windows and sasa calc

