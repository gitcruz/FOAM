# FOAM
FOAM (Flye ONT Assembly of Mitogenomes). Is a pipeline we are using to assemble the mitogenome of different organisms using ONT and Illumina reads. As each genome and sample poses different challenges is in continuous development

# Current release
Version 0.2 was successfully used to assemble the mitogenomes of the pearl razorfish (_Xyrichtys novacula_) and the common octopus (_Octopus vulgaris_).

To make it run:
1. Clone this repository.
2. Go to your local utils/refseq_dir.
3. Uncompress the refseq databases locally, example:  tar jvxf refseq81m.tar.bz2
4. Create config files for your genome (check .yaml and .spec for slurm in example/configs).
5. Run the pipeline with snakemake pointing to the snakefile in bin.
   
# References
A chromosome-level reference genome for the common octopus, _Octopus vulgaris_ (Cuvier, 1797)
Dalila Destanović, ProfileDarrin T. Schultz, Ruth Styfhals, Fernando Cruz, Jèssica Gómez-Garrido, Marta Gut, Ivo Gut, Graziano Fiorito, Oleg Simakov, Tyler S. Alioto, Giovanna Ponte, Eve Seuntjens.
doi: https://doi.org/10.1101/2023.05.16.540928 
