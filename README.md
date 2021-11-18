# Snakemake-SNP-Calling-Pipeline
A Snakemake version of the i18 SNP calling pipeline.

# Description:
The pipeline runs snippy (and associated snippy-core/snippy-clean), removes reference/outgroups as required, iqtree, cfml, maskrc-svg, and snp-dists.

Note that snp-dists produces a SNP matrix but not p-distances. Use MegaX for that (or perhaps a future update if I find a tool for it).

# Tutorial:
To come...


## Resources for snakemake and examples of other pipelines:
BacDist: A relatively new snp pipeline (https://github.com/MigleSur/BacDist)<br>
Example of a VERY clean snakemake workflow (https://github.com/snakemake-workflows/rna-seq-star-deseq2)<br>
A custom snp pipeline running similar tools (an example of how to incorporate them; for learning purposes) (https://github.com/CJREID/snplord)
<br>
<br>
A great tutorial for the basics of snakemake: https://lachlandeer.github.io/snakemake-econ-r-tutorial/initial-steps-with-snakemake.html
<br>
Another tutorial that gets into some nitty-gritty: https://carpentries-incubator.github.io/workflows-snakemake/
<br>
Some mentions of how to handle wildcard constraints and when you need to: https://edwards.sdsu.edu/research/wildcards-in-snakemake/
<br><br>
Another carpentries tutorial: https://carpentries-incubator.github.io/snakemake-novice-bioinformatics/
<br><br>
A list of tutorials for snakemake: https://edwards.sdsu.edu/research/snakemake-tutorial/
