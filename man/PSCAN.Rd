\name{PSCAN}
\alias{PSCAN}

\title{
Protein-structure-guided scan tests for detecting gene-level associations and signal variants
}
\description{
This function allows users to (a) perform PSCAN gene-level association tests; (b) run PSCAN procedures to identify trait-associated variants.
}
\usage{
PSCAN(type="mean", p.comb="Cauchy", U, V, MAC, weight=NULL, N.MC = 1000, FWER = 0.05, 
f=0, details=FALSE, plot3D = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
type="mean" is the scan test for testing the mean of the effects and type="variance" is the scan test for testing the variance
}
  \item{p.comb}{
p.comb="Cauchy" or "minP" is to use Cauchy method or minimum p-value method to obtain the gene-level p-value
}
  \item{U}{
a vector of score statistics with names as SNP ID in chr:position format
}
  \item{V}{
a covariance matrix for the score statistics
}
  \item{MAC}{
a vector of minor allele count
}
  \item{weight}{
a vector contains the weights for SNPs (default = NULL, flat weight)
}
  \item{N.MC}{
the number of Monte Carlo simulations for minP gene-level test and signal region detection. 
}
  \item{FWER}{
family wise error rate for the PSCAN signal region detection
}
  \item{f}{
overlap fraction cutoff in the signal detection algorithm (0<=f<=1). When f = 0, algorithm identifies nonoverlapping signal regions. When f = 1, algorithm essentially
keeps every region passing the significance threshold as the detected signal regions.
}
  \item{details}{
whether output additional information, including p-values for individual signal regions, the type of protein sturcture used in the algorithm (experimentally determined protein structure ("PDB") or computationally predicted structure ("Modbase")), and the 3D coordinates of variants in the protein space.
}
  \item{plot3D}{
whether generate interactive 3D plot for signal regions in protein space. In the plot, the purple spheres represent signal variants and grey spheres represent non-signal variants. The plot will not be generated if no signal region is detected
}
}
\value{
\item{pscan.pval}{a p-value of PSCAN mean-based or variance-based test}
\item{signal}{a list of vector with each vector contains SNP IDs in a signal region}
\item{signal.pval}{a vector of p-values with each p-value corresponds to the mean-based or variance-based test p-value for a signal region}
}
\references{
Tang ZZ, Sliwoski GR, Chen G, Jin B, Bush WS, Li B, and Capra JA. (2020). Spatial scan tests guided by protein structures improve complex disease gene discovery and signal variant detection. Genome Biology.

}
\author{
Zheng-Zheng Tang
}

\examples{
data(one.gene)
U = one.gene$U
V = one.gene$V
MAC = one.gene$MAC
weight = one.gene$weight
# perform PSCAN mean-based test and signal region detection
PSCAN(type="mean", U=U, V=V, MAC=MAC, weight=weight, N.MC=1000, FWER=0.05, plot3D=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PSCAN }
