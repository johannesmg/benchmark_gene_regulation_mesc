# benchmark_gene_regulation_mesc

This repository contains code for benchmarking gene regulatory network
reconstruction algorithms in R. It contains code from raucpr
(https://github.com/kboyd/raucpr) in the raucpr directory. This
directory also contains the license under which reuse is
permitted. Additionally code from the ROCR package
(https://github.com/ipa-tys/ROCR) is used. Here the GPL (>=2) applies.

Be warned that because of the memory inefficient implementation, the
computation may require 10s of gigabytes of memory, depending on the
number of cores used.

In order to compute the benchmark, an expression set matrix is
needed. A collection of prenormalized mouse embryonic stem cell
expression data from the GEO database can be downloaded at:
