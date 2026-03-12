# Runtime and peak ram usage #

#install.packages("peakRAM")
library(peakRAM)
library(evo3D)
library(tidyverse)

# read in analysis data #
msa1 = 'ChikV/chikv_e1hits.fa'
msa2 = 'ChikV/chikv_e2hits.fa'
pdb = 'ChikV/8fcg.pdb'

##### first record times and memory for centroid ####

pc = list(
  max_patch = 15,
  dist_cutoff = NA,
  distance_method = 'centroid',
  rsa_cutoff = NA,
  sasa_cutoff = 10
)

peakRAM({
  res = run_evo3d(list(msa1, msa2), pdb, chain = c('A', 'E'),
                  pdb_controls = pc,
                  detail = 1, analysis_mode = 'residue')
})

system.time({
  res = run_evo3d(list(msa1, msa2), pdb, chain = c('A', 'E'),
                  pdb_controls = pc,
                  detail = 1, analysis_mode = 'residue')
})

# recording for each run time (t1, t2) and peak memory usage (ram_mb1, ram_mb2) #
# each chain combination was run twice to record the average runtime and usage #
tib = tibble(distance_method = 'centroid',
             chains = c('A,E', 'AB,EF', 'ABC,EFG', 'ABCD,EFGH'),
             residues = c(858,1716,2574,3432),
             atoms = c(6622,13244,19866,26488),
             windows = c(653,1298,1925,2573),
             ram_mb1 = c(259.8,389.4,440.4,539.7),
             ram_mb2 = c(240.8,389.4,443.9,577.7),
             t1 = c(5.69,7.32,10.44,13.84),
             t2 = c(5.49,7.65,10.25,13.79)
)

tib$mean = (tib$t1 + tib$t2) / 2
tib$mem = (tib$ram_mb1 + tib$ram_mb2) / 2

# these records go to supplementary table S4
tib %>% select(distance_method, residues, atoms, windows, mean, mem)

#### record times and usage for all atom ####

pc = list(
  max_patch = 15,
  dist_cutoff = NA,
  distance_method = 'all',
  rsa_cutoff = NA,
  sasa_cutoff = 10
)

peakRAM({
  res = run_evo3d(list(msa1, msa2), pdb, chain = c('A', 'E'),
                  pdb_controls = pc,
                  detail = 1, analysis_mode = 'residue')
})

system.time({
  res = run_evo3d(list(msa1, msa2), pdb, chain = c('A', 'E'),
                  pdb_controls = pc,
                  detail = 1, analysis_mode = 'residue')
})


# recording for each run time (t1, t2) and peak memory usage (ram_mb1, ram_mb2) #
# each chain combination was run twice to record the average runtime and usage #
tib2 = tibble(distance_method = 'all',
              chains = c('A,E', 'AB,EF', 'ABC,EFG', 'ABCD,EFGH'),
              residues = c(858,1716,2574,3432),
              atoms = c(6622,13244,19866,26488),
              windows = c(653,1298,1925,2573),
              ram_mb1 = c(747.8,3281.6,7305.6,11076.5),
              ram_mb2 = c(862.2,3275.9,6273.3,11076.5),
              t1 = c(7.72,16.47,30.92,86.50),
              t2 = c(7.61,15.39,30.22,95.83)
)

tib2$mean = (tib2$t1 + tib2$t2) / 2
tib2$mem = (tib2$ram_mb1 + tib2$ram_mb2) / 2

# these records go to supplementary table S4
tib2 %>% select(distance_method, residues, atoms, windows, mean, mem)

#### check size of residues and atoms per chain combination ####

# each chain combination was checked as follows #
check = res$pdb_info_sets$pdb1$pdb
check2 = bio3d::trim.pdb(check, chain = c('A', 'E',
                                    'B', 'F',
                                    'C', 'G',
                                    'D', 'H'))
