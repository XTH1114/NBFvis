<!-- README.md is generated from README.Rmd. Please edit that file -->

# NBFvis

<!-- badges: start -->
<!-- badges: end -->

Dimensionality reduction of spatial omic data can reveal shared,
spatially structured patterns of expression across a collection of
genomic features. We study strategies for discovering and interactively
visualizing low-dimensional structure in spatial omic data based on the
construction of neighborhood features. We design quantile and
network-based spatial features that result in spatially consistent
embeddings. A simulation compares embeddings made with and without
neighborhood-based featurization, and a re-analysis of \[Keren et al.,
2019\] illustrates the overall workflow. We provide an R package,
NBFvis, to support computation and interactive visualization for the
proposed dimensionality reduction approach.

![UMAP Embedding Plot and Spatial Plot](../example_plot_for_pkgdown.png)

## Installation

You can install the released version of NBFvis from
[GitHub](https://github.com/XTH1114/NBFvis).

    devtools::install_github("XTH1114/NBFvis")
