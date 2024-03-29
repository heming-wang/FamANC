# FamANC
FamANC estimates local ancestry in large pedigrees by: 
1. Using an existing software (e.g. SABER+ and HAPMIX) to infer local ancestries for all individuals in a family, temporarily assuming they are unrelated

2. Using FamANC to improve the local ancestry inference accuracy by using the known pedigree structure to correct inference errors that arise from
   - Mendelian inconsistencies identified from the pedigree structure
   - Double crossovers occurring within 2cM

# R Functions
- Genetic.Distance
- Mendelian.Path
- Prob.Error
- Select.Path
- Mendelian.Anc
- Pedigree.Divide
- plotped

### Require R package kinship2
```
library("kinship2")
```

### Genetic.Distance(map,gmap)
Estimate genetic distance between neighboring markers for a set of variants on the same chromosome (map). 
A reference genetic map (gmap) is needed.

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
Generate genetic distance between markers
```
Morgan=Genetic.Distance(sim.map,ref.gmap)

Morgan[1:10]
##  [1] 0.0002247961 0.0003491192 0.0003883777 0.0004514758 0.0005432497
##  [6] 0.0006019482 0.0006746992 0.0006804014 0.0007710180 0.0009503945
```

### Mendelian.Path(ped)
Mendelian paths are all possible local ancestry markers at a position in a pedigree (ped) satisfying Mendelian inheritance.
For example, in a nuclear family with two parents and one offspring (n=3), the Mendelian paths corresponding to (Father, Mother, Child) = {(0,0,0), (0,1,0), (0,1,1), (0,2,1), (1,0,0), (1,0,1), (1,1,0), (1,1,1), (1,1,2), (1,2,1), (1,2,2), (2,0,1), (2,1,1), (2,1,2,), (2,2,2)}. Function Mendelian.Path returns to a reordered ped file and the Mendelian paths with the corresponding order.

Load fam file in PLINK format
```
load("data/sim.fam.Rdata") ### FID IID PAT MAT SEX PHE

sim.fam ## 10 individuals
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
Generate Mendelian path for sim.fam
```
mPath=Mendelian.Path(sim.fam) 

ls(mPath)
## [1] "ped"      "ped.path"

dim(mPath$ped.path)
## [1] 3615   10

mPath$ped.path[1:5,]
##    201 202 203 207 204 205 206 208 209 210
## 1    0   0   0   0   0   0   0   0   0   0
## 2    1   0   0   0   0   0   0   0   0   0
## 4    0   1   0   0   0   0   0   0   0   0
## 5    1   1   0   0   0   0   0   0   0   0
## 10   0   0   1   0   0   0   0   0   0   0
```
### Prob.yx(xt,yt,epsilon)
The probability of an inferred local ancestry path (yt) given a Mendelian path (xt) given the average allelic local ancestry inference error rate epsilon.
```
xt=mPath$ped.path[1,]
yt=c(2,2,1,2,2,1,2,1,2,2)

pr=Prob.yx(xt,yt,epsilon=0.01); pr
## [1] 7.684768e-32
```
### Select.Path(yt,ped.path,thres,epsilon)
The probabilities of observing many >=3 individuals with ancestry error at the same locus were small in simulated pedigrees. To save computation time, we select a smaller set of Mendelian paths (ped.path) with a threshold number (thres=2) of values different from the observed path (yt). Function Select.Path returns to the selected Mendelian paths and the corrsponding probabilities of observing yt.
```
sPath=Select.Path(yt,mPath$ped.path,thres=2,epsilon=0.01)

ls(sPath)
## [1] "pr"          "select.path"
```
### Mendelian.Anc(anc,mPath,Morgan,t,thres,epsilon)
Correct local ancestries estimated from existing software (anc). t is the number of generations since admixture.

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
Correct errors in local ancestry
```
ancs1=Mendelian.Anc(sim.anc,mPath,Morgan,t=8,thres=2,epsilon=0.01)
```
### Pedigree.Divide(big.ped)
For some pedigrees with missing first-generation genotype, we removed the first generation and divided them into smaller pedigrees.
```
ped1=Pedigree.Divide(big.ped)
```
### plotped(ped)
Plot pedigree
```
plotped(ped)
```

# Examples
### Source R functions
```
source("FamANC.functions.R")
```
### Import data
```
load("data/sim.fam.Rdata")
load("data/sim.map.Rdata") 
load("data/sim.anc.Rdata")
load("data/genetic_map_GRCh37_chr22.RData")
```
### Run FamANC
```
Morgan=Genetic.Distance(sim.map,ref.gmap)
mPath=Mendelian.Path(sim.fam)
ancs1=Mendelian.Anc(sim.anc,mPath,Morgan,t=8,thres=2,epsilon=0.01)

```
