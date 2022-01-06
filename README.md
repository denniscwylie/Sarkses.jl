# Sarkses.jl
Julia implementation of suffix array kernel smoothing (SArKS) algorithms for identifying sequence motifs correlated with numeric scores, as well as for motif-based feature extraction for predictive modeling of scores based on sequences.

For details on the original SArKS approach for *de novo* correlative motif identification see:

https://academic.oup.com/bioinformatics/article-abstract/35/20/3944/5418797

A preprint of the article is also available on biorxiv at:

https://www.biorxiv.org/content/early/2018/10/25/133934

Note that there is also an Bioconductor package in R implementing the original SArKS motif identification algorithm available at http://www.bioconductor.org/packages/release/bioc/html/sarks.html, though there are additional features in the Julia implementation here not available in that package.

## Installation
After downloading or cloning this git repository, change directories into the Sarkses.jl folder, then start up the Julia REPL and type `]add .` and then `using Sarkses` to access SArKS functionality.
