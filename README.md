<br><b># PanGeniePipeline</b></br> 
- a Snakemake pipeline that uses PanGenie - a graph-based genotyper - to infer genotypes of variants based on a pan-genome graph and short reads. 
- This pipeline first takes a unmerged, phased VCF to make a pan-genome VCF, named pangenome.vcf.gz, which is required for PanGenie as its input.
- For HGSVC + HPRC (n=77 diploid, phased genomes) and a 30X Illumina short-read genome, it takes 12 cores, 25GB per core, and ~20hrs to genotype the genome.
<br></br>
To run the pipeline:
  snakemake --drmaa " -l centos=7 -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -j 20 -w 60 -p
