# Snakemake-SNP-Calling-Pipeline
A Snakemake version of the i18 SNP calling pipeline.

# Description:
The pipeline runs snippy (and associated snippy-core/snippy-clean), removes reference/outgroups as required, iqtree, cfml, maskrc-svg, and snp-dists.

Note that snp-dists produces a SNP matrix but not p-distances. Use MegaX for that (or perhaps a future update if I find a tool for it).

# Tutorial:
To come...


## Resources for snakemake and examples of other pipelines:
\nBacDist: A relatively new snp pipeline (https://github.com/MigleSur/BacDist)
\nExample of a VERY clean snakemake workflow (https://github.com/snakemake-workflows/rna-seq-star-deseq2)
\nA custom snp pipeline running similar tools (an example of how to incorporate them; for learning purposes) (https://github.com/CJREID/snplord)
