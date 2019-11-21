# FamANC
FamANC estimates local ancestry in large pedigrees by: 
1. Using an existing software (e.g. SABER+ and HAPMIX) to infer local ancestries for all individuals in a family, temporarily assuming they are unrelated

2. Using FamANC to improve the local ancestry inference accuracy by using the known pedigree structure to correct inference errors that arise from
   - Mendelian inconsistencies identified from the pedigree structure
   - double crossovers occurring within 2cM


# Examples
### Load required package
```
library("kinship2")
```
### Import data
Load map file in PLINK format
```
load("data/sim.map.Rdata") ### Chr rsid cM bp

head (sim.map) 
##   V1         V2           V3       V4
## 1 22  rs4819973 0.0002247961 16054117
## 2 22  rs5992637 0.0003491192 16055900
## 3 22 rs11912507 0.0003883777 16056590
## 4 22  rs5747018 0.0004514758 16057699
## 5 22  rs1076103 0.0005432497 16059312
## 6 22 rs17807317 0.0006019482 16060519
```
Load fam file in PLINK format
```
load("data/sim.fam.Rdata") ### FID IID PAT MAT SEX PHE

sim.fam 
##    V1  V2  V3  V4 V5 V6
## 21  2 201   0   0  1 -9
## 22  2 202   0   0  2 -9
## 23  2 203   0   0  1 -9
## 24  2 207   0   0  1 -9
## 25  2 204 201 202  2 -9
## 26  2 205 201 202  1 -9
## 27  2 206 201 202  2 -9
## 28  2 208 203 204  1 -9
## 29  2 209 203 204  1 -9
## 30  2 210 207 206  2 -9
```
Load local ancestry file inferred by existing software (e.g. SABER+ and HAPMIX)
```
load("data/sim.anc.Rdata") ### Rows: individuals; Columns: variants

dim(sim.anc)
## [1]    10 18210

sim.anc[1:5,1:5]
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    1    1    1    1
## [2,]    2    2    2    2    2
## [3,]    1    1    1    1    1
## [4,]    2    2    2    2    2
## [5,]    2    2    2    2    2
```
Load reference genetic map
```
load("data/genetic_map_GRCh37_chr22.RData")

head(ref.gmap)
##   Chromosome Position.bp. Rate.cM.Mb.  Map.cM.
## 1      chr22     16051347    8.096992 0.000000
## 2      chr22     16052618    8.131520 0.010291
## 3      chr22     16053624    8.131967 0.018472
## 4      chr22     16053659    8.132625 0.018756
## 5      chr22     16053758    8.129606 0.019561
## 6      chr22     16054713    8.024772 0.027325
```
### Run FamANC
Generate genetic distance between markers
```
Morgan=Genetic.Distance(sim.map,ref.gmap)

Morgan[1:10]
##  [1] 0.0002247961 0.0003491192 0.0003883777 0.0004514758 0.0005432497
##  [6] 0.0006019482 0.0006746992 0.0006804014 0.0007710180 0.0009503945

```
Generate Mendelian paths for a given pedigree
```
mPath=Mendelian.Path(sim.fam)
```
Correct errors in local ancestry
```
ancs1=Mendelian.Anc(sim.anc,mPath,Morgan,t=6,thres=2,epsilon=0.1)

dim(ancs1)
## [1]    10 18210
```
