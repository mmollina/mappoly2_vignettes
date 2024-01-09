## ----setup, include = FALSE, eval = TRUE-------------------------------------------------------------------
setwd("~/repos/training/mappoly2_vignettes/alfalfa/")
require(mappoly2)
## ----loading dataset---------------------------------------------------------------------------------------
alfalfa.f1 <- read_geno_csv(file.in = "I195_x_J432.csv",
                            ploidy.p1 = 4,
                            ploidy.p2 = 4,
                            name.p1 = "I195",
                            name.p2 = "J432")


## ----print plot dataset------------------------------------------------------------------------------------
alfalfa.f1
plot(alfalfa.f1)

## ----bundled dataset, eval=FALSE---------------------------------------------------------------------------
## alfalfa.f1 <- alfa_f1

## ----QA_QC_miss--------------------------------------------------------------------------------------------
alfalfa.f1 <- filter_data(alfalfa.f1, mrk.thresh = 0.2, ind.thresh = 0.1)
plot(alfalfa.f1)
alfalfa.f1

## ----filter ind, eval=FALSE--------------------------------------------------------------------------------
alfalfa.f1 <- filter_individuals(alfalfa.f1)

## ----filter ind4-------------------------------------------------------------------------------------------
alfalfa.f1

## ----------------------------------------------------------------------------------------------------------
plot(alfalfa.f1, type = "density")

## ----two_points,fig.width = 6, fig.height= 6---------------------------------------------------------------
## 55 sec
t1 <- system.time(alfalfa.f1.all <- pairwise_rf(alfalfa.f1, mrk.scope = "all", ncpus = 8))
plot(alfalfa.f1.all)
## 8 sec
t2 <- system.time(alfalfa.f1.chr <- pairwise_rf(alfalfa.f1, mrk.scope = "per.chrom", ncpus = 8))
plot(alfalfa.f1.chr)
rbind(t1, t2)

## ----rf_based_filter1--------------------------------------------------------------------------------------
alfalfa.f1.all <- rf_filter(alfalfa.f1.all,
                            thresh.LOD.ph = 5,
                            thresh.LOD.rf = 5,
                            thresh.rf = 0.15,
                            probs = c(0.025, 0.975))
alfalfa.f1.all
plot(alfalfa.f1.all)

## ----groupping---------------------------------------------------------------------------------------------
g <- group(x = alfalfa.f1.all, expected.groups = 8, comp.mat = TRUE)
g
plot(g)

## ----make_seq----------------------------------------------------------------------------------------------
s <- make_sequence(g,
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8),
                   lg = list(1, 7, 6, 4, 2, 8, 3, 5))
print(s, type = "mds")
print(s, type = "genome")

## ----ordering----------------------------------------------------------------------------------------------
s <- order_sequence(s, type = "mds")
print(s, type = "mds")
s <- order_sequence(s, type = "genome")
print(s, type = "genome")

## ----comparing orders--------------------------------------------------------------------------------------
plot_rf_matrix(s, type = "mds", fact = 2)
plot_rf_matrix(s, type = "genome", fact = 2)
plot_mds_vs_genome(s)

## ----rf_based_filter2--------------------------------------------------------------------------------------
#### RF-based filter per groups ####
s <- rf_filter(s, type = "mds", probs = c(0.025, 0.975), diag.markers = 50)
s <- rf_filter(s, type = "genome", probs = c(0.025, 0.975), diag.markers = 50)
mappoly2:::plot_rf_matrix(s, type = "genome", fact = 2)


## ----rf_phasing--------------------------------------------------------------------------------------------
s <- pairwise_phasing(s, type = "mds",
                      thresh.LOD.ph = 3,
                      thresh.LOD.rf = 3,
                      thresh.rf = 0.5,
                      max.search.expansion.p1 = 10,
                      max.search.expansion.p2 = 10)
print(s, type = "mds")
s <- pairwise_phasing(s,
                      type = "genome",
                      thresh.LOD.ph = 3,
                      thresh.LOD.rf = 3,
                      thresh.rf = 0.5,
                      max.search.expansion.p1 = 10,
                      max.search.expansion.p2 = 10)
print(s, type = "genome")


## ----hmm_single_p1-----------------------------------------------------------------------------------------
s <- mapping(s, type = "mds", parent = "p1", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p1", ncpus = 8)
print(s, type = "mds")
print(s, type = "genome")
plot_map(s, lg = 1, type = "mds", parent = "p1")
plot_map(s, lg = 1, type = "genome", parent = "p1")


## ----hmm_single_p2-----------------------------------------------------------------------------------------
s <- mapping(s, type = "mds", parent = "p2", ncpus = 8)
s <- mapping(s, type = "genome", parent = "p2", ncpus = 8)
print(s, type = "mds")
print(s, type = "genome")
plot_map(s, lg = 1, type = "mds", parent = "p2")
plot_map(s, lg = 1, type = "genome", parent = "p2")


## ----merge_maps--------------------------------------------------------------------------------------------
## 28  sec
system.time(s <- merge_single_parent_maps(s, type = "mds", ncpus = 8, error = 0.05))
## 35  sec
system.time(s <- merge_single_parent_maps(s, type = "genome", ncpus = 8, error = 0.05))
plot_map_list(s, type = "mds", parent = "p1p2")
plot_map_list(s, type = "genome", parent = "p1p2")
map_summary(s)
plot_map(s, lg = 1, type = "mds", parent = "p1p2")
plot_map(s, lg = 1, type = "genome", parent = "p1p2")
s <- calc_haplotypes(s, type = "mds", ncpus = 8)
s <- calc_haplotypes(s, type = "genome", ncpus = 8)


## ----augment_merged_maps-----------------------------------------------------------------------------------
## 51 sec
system.time(s <- augment_phased_map(s, type = "mds", ncpus = 8))
## 56 sec
system.time(s <- augment_phased_map(s, type = "genome", ncpus = 8))
plot_map_list(s, type = "mds", parent = "p1p2")
plot_map_list(s, type = "genome", parent = "p1p2")
map_summary(s)
plot_map(s, lg = 1, type = "mds", parent = "p1p2")
plot_map(s, lg = 1, type = "genome", parent = "p1p2")


## ----compare_mds_geno--------------------------------------------------------------------------------------
maps.comp <- compare_order(s)
maps.comp
plot(maps.comp)


## ----recompute_haplotype_probs-----------------------------------------------------------------------------
## 37 sec
s <- calc_haplotypes(s, type = "mds", ncpus = 8)
s <- calc_haplotypes(s, type = "genome", ncpus = 8)
plot_haplotypes(s, lg = 1, ind = "F1.85.21")
I195xJ432_map <- s


## ----final_f1----------------------------------------------------------------------------------------------
map_summary(I195xJ432_map, type = "genome")
plot_map_list(I195xJ432_map, type = "genome", col = mp_pal(8))


## ----load_bc-----------------------------------------------------------------------------------------------
download.file("https://github.com/mmollina/mappoly2_vignettes/raw/main/I195_x_F1-85-209_map.rda",
              destfile = "temp_file.rda")
load("temp_file.rda")
## 23 sec
I195xF1_85_209_map <- mapping(I195xF1_85_209_map, type = "genome", error = 0.05, ncpus = 8)
print(I195xF1_85_209_map, type = "genome")
map_summary(I195xF1_85_209_map, type = "genome")
plot_map_list(I195xF1_85_209_map, type = "genome", col = mp_pal(8))


## ----prepare1----------------------------------------------------------------------------------------------
MAPs <- list(I195xJ432 = I195xJ432_map,
             I195xF1_85_209 = I195xF1_85_209_map)
plot_multi_map(MAPs)
prep.maps <- prepare_to_integrate(MAPs)
plot(prep.maps)
plot_shared_markers(prep.maps)

## ----est_consensus------------------------------------------------------------------------------------------
##93 sec
consensus.map <- estimate_consensus_map(prep.maps, ncpus = 8, err = 0.05)
consensus.map
plot(consensus.map)
plot(consensus.map, only.consensus = TRUE, col = mp_pal(8))


## ----comp_cons_hap_prob------------------------------------------------------------------------------------
consensus.map <- calc_consensus_haplo(consensus.map, ncpus = 8)
plot_consensus_haplo(consensus.map, lg = 1, ind = "F1.85.70")
plot_consensus_haplo(consensus.map, lg = 1, ind = "AphBC.21")


## ----save_image, echo=FALSE, results='hide'----------------------------------------------------------------
save.image(file = "all.rda")

