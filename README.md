# KeBABS - An R Package for Kernel-Based Analysis of Biological Sequences
The 'kebabs' package provides functionality for kernel based analysis of biological sequences via Support Vector Machine (SVM) based methods. Biological sequences include DNA, RNA, and amino acid (AA) sequences. Sequence kernels define similarity measures between sequences. The package implements some of the most important kernels for sequence analysis in a very flexible and efficient way and extends the standard position-independent functionality of these kernels in a novel way to take the position of patterns in the sequences into account for the similarity measure.

Although the package is maintained by Ulrich Bodenhofer, the package itself
has been implemented by Johannes Palme.

## Installation

The package can be installed from
[Bioconductor](https://bioconductor.org/). Therefore, the the simplest way to install the package is to enter
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("kebabs")
```
into your R session. If, for what reason ever, you prefer to install the package manually, follow the instructions in the [user manual](https://bioconductor.org/packages/release/bioc/vignettes/kebabs/inst/doc/kebabs.pdf).

## User support

If you encounter any issues or if you have any question that might be of interest also for other users, before writing a private message to the package developers/maintainers, please create an issue in this repository and also consider posting on [Bioconductor Support](https://support.bioconductor.org/) or on [StackOverflow](https://stackoverflow.com/). For other matters regarding the package, please contact the package author(s).

## Citing this package

If you use this package for research that is published later, you are kindly asked to cite it as follows:

- J. Palme, S. Hochreiter, and U. Bodenhofer (2015). KeBABS: an R package for kernel-based analysis of biological sequences. *Bioinformatics* **31**:2574-2576. DOI: [10.1093/bioinformatics/bty176](http://doi.org/10.1093/bioinformatics/bty176](http://doi.org/10.1093/bioinformatics/btr406)