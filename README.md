# FamANC
FamANC estimates local ancestry in large pedigrees by: 
1. Using an existing software (e.g. SABER+ and HAPMIX) to infer local ancestries for all individuals in a family, temporarily assuming they are unrelated

2. Using FamANC to improve the local ancestry inference accuracy by using the known pedigree structure to correct inference errors that arise from
   - Mendelian inconsistencies identified from the pedigree structure
   - double crossovers occurring within 2cM


# Examples
## Load required package
```
library("kinship2")
```
## Import data
Load map and fam file in PLINK format
```
load("sim.map.Rdata") ### FID IID PAT MAT SEX PHE
load("data/sim.map.Rdata") ### Chr rsid cM bp
```
Load local ancestry file inferred by existing software (e.g. SABER+ and HAPMIX)
```
load("data/sim.anc.Rdata")
```
Load reference genetic map
```
load("data/genetic_map_GRCh37_chr22.RData") 
```
## Run FamANC
Generate genetic distance between markers
```
Morgan=Genetic.Distance(sim.map,ref.gmap)
```
Generate Mendelian paths for a given pedigree
```
mPath=Mendelian.Path(sim.fam)
```
Correct errors in local ancestry
```
ancs1=Mendelian.Anc(sim.anc,mPath,Morgan,t=6,thres=2,epsilon=0.1)
```
