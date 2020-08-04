### PSCAN v1.1

**PSCAN** Protein-structure-guided scan methods for gene-level association test and signal region detection

**PSCAN** package has the main PSCAN function that implements protein-structure-guided scan (PSCAN)
methods for detecting gene-level associations and signal variants. PSCAN methods leverage the tendency of
functional variants to cluster in 3D protein space. PSCAN methods are built upon flexibly shaped spatial
scan statistics, with scan windows adaptively defined to accommodate diverse topologies of variant positions
in protein space. PSCAN performs fast gene-level association tests by combining SNP-set-based testing
p-values across windows using the Cauchy method. In addition, PSCAN implements an efficient search
algorithm for the detection of multiple signal regions in protein space. The method details are described in Tang et
al. (2020).

### References

Tang ZZ, Sliwoski GR, Chen G, Jin B, Bush WS, Li B, and Capra JA (2020). PSCAN: Spatial scan tests guided by protein structures improve complex disease gene discovery and signal variant detection. Genome Biology.


### Installation

**PSCAN** can be installed from github directly as follows:

```r
install.packages("devtools")
library(devtools)
install_github("tangzheng1/PSCAN")
```

### Tutorial

A vignette with examples illustrating the usage of **PSCAN** is available at: https://github.com/tangzheng1/PSCAN/blob/master/vignette

### Contact

**Zheng-Zheng Tang** (UW-Madison) tang@biostat.wisc.edu.

