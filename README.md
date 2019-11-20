# FamANC
FamANC estimates local ancestry in large pedigrees by: 1) using an existing software (e.g. SABER+ and HAPMIX) to infer local ancestries for all individuals in a family, temporarily assuming they are unrelated, and then 2) using FamANC to improve the local ancestry inference accuracy by using the known pedigree structure to correct inference errors.

# Examples
## Load map file in PLINK format
load("sim.map.Rdata")
head(sim.map)

