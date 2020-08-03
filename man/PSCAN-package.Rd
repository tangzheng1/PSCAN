\name{PSCAN-package}
\alias{PSCAN-package}
\alias{PSCAN}
\docType{package}
\title{
Protein-structure-guided scan tests for detecting gene-level associations and signal variants
}
\description{
PSCAN package has the main PSCAN function that implements protein-structure-guided scan (PSCAN)
methods for detecting gene-level associations and signal variants.
}
\details{
\tabular{ll}{
Package: \tab PSCAN\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2020-07-30\cr
License: \tab GPL (>=2)\cr
}
PSCAN
}
\author{
Zheng-Zheng Tang

Maintainer: Zheng-Zheng Tang <tang@biostat.wisc.edu>

}
\references{
Tang ZZ, Sliwoski GR, Chen G, Jin B, Bush WS, Li B, and Capra JA. (2020). Spatial scan tests guided by protein structures improve complex disease gene discovery and signal variant detection. Genome Biology.


}

\keyword{ package }
\examples{
data(one.gene)
U = one.gene$U
V = one.gene$V
MAC = one.gene$MAC
weight = one.gene$weight
# perform PSCAN mean-based test and signal region detection
PSCAN(type="mean", U=U, V=V, MAC=MAC, weight=weight, N.MC=1000, FWER=0.05, plot3D=TRUE)
}
