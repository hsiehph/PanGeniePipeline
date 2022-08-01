<br><b># PanGeniePipeline</b></br> 
- a Snakemake pipeline that uses PanGenie - a graph-based genotyper - to infer genotypes of variants based on a pan-genome graph and short reads. 
- This pipeline first takes a unmerged, phased VCF to make a pan-genome VCF, named pangenome.vcf.gz, which is required for PanGenie as its input.
- For a graph from the HGSVC + HPRC (n=77) diploid, phased assemblies, it takes 12 cores, 25GB per core, and ~20hrs to genotype a 30X Illumina short-read genome
<br></br>
<b>To run the pipeline:</b>
<br> snakemake --drmaa " -l centos=7 -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -j 20 -w 60 -p </br>
