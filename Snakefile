import os

#shell.prefix("source config.sh; set -eo pipefail ; ")
shell.prefix("source env.cfg ; set -eo pipefail ;")

configfile:'config.yaml'

if not os.path.exists("log"):
	os.makedirs("log")

## parameters
vcf = config["vcf"]
reference = config["reference"]
reference_shortRead = config["reference_shortRead"]
genotypeRegion = config["genotypeRegion"]

scripts = "scripts"
# allow max 20 % missing alleles at a variant position
frac_missing = config["frac_missing"]
kmer_size = config["kmer_size"]

dict_sample_type = {}
for k in config["shortReadGenome"].keys():
	for sample in config["shortReadGenome"][k]:
		dict_sample_type[sample] = k


rule all:
	input:
		expand('output/{SAMPLE}_genotyping.vcf.gz', SAMPLE=dict_sample_type.keys()),


###########################################################
##  Convert VCF into pangenome graph by merging variants ##
##  that are overlapping into multi-allelic positions.   ##
###########################################################

# check that VCF is correct and remove positions with more than frac_missing missing alleles
rule prepare_vcf:
	input:
		vcf
	output:
		'mergeInputVcf/input-missing-removed.vcf'
	benchmark:
		'mergeInputVcf/prepare-vcf.benchmarks'
	log:
		'log/prepare-vcf.log'
	threads: 1
	conda:
		'env/merging.yml'
	resources:
		mem=200,
		hrs=48,
	shell:
		"which python 2> {log}; zcat {input} | python3 {scripts}/prepare-vcf.py --missing {frac_missing} 2> {log} 1> {output}"


# assign IDs to all alleles
rule add_ids:
	input:
		'mergeInputVcf/input-missing-removed.vcf'
	output:
		'mergeInputVcf/callset.vcf'
	benchmark:
		'mergeInputVcf/add-ids.benchmarks'
	log:
		'log/callset.log'
	conda:
		'env/merging.yml'
	resources:
		mem=20,
		hrs=48,
	shell:
		"cat {input} | python {scripts}/add-ids.py 2> {log} 1> {output}"


# create biallelic VCF with one record per ALT allele
rule normalize:
	input:
		'mergeInputVcf/callset.vcf'
	output:
		'mergeInputVcf/callset-biallelic.vcf'
	benchmark:
		'mergeInputVcf/callset-biallelic.benchmarks'
	log:
		'log/callset-biallelic.log'
	conda:
		'env/merging.yml'
	resources:
		mem=20,
		hrs=48,
	shell:
		"bcftools norm -m- {input} 2>> {log} 1> {output}"


# merge variants into a pangenome graph
rule merge_haplotypes:
	input:
		vcf = 'mergeInputVcf/callset-biallelic.vcf',
		reference = reference
	output:
		'mergeInputVcf/pangenome.vcf.gz'
	benchmark:
		'mergeInputVcf/merge-haplotypes.benchmarks'
	log:
		 'log/pangenome.log'
	conda:
		"env/merging.yml"
	resources:
		mem=20,
		hrs=48,
	shell:
		"""
		python3 {scripts}/merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2  2> {log} | bgzip -c > {output}
		"""

rule PanGenie:
	input:
		pangenome='mergeInputVcf/pangenome.vcf.gz',
		reference = reference,
		reference_shortRead = reference_shortRead
	output:
		temp('output/{SAMPLE}_genotyping.vcf'),
		temp('output/{SAMPLE}_histogram.histo'),
		temp('output/{SAMPLE}_path_segments.fasta'),
	log:
		 'log/genotype.{SAMPLE}.log'
	benchmark:
		'benchmark/genotype.{SAMPLE}.benchmarkds'
	threads: 12
	resources:
		mem=25,
		hrs=72,
	run:
		sample = list({wildcards.SAMPLE})[0]
		inFileType = dict_sample_type[sample]
		shortread = config["shortReadGenome"][inFileType][sample]

		if genotypeRegion:
			print("genotype by region")
			if inFileType not in ["BAM", "CRAM"]:
				raise TypeError("Region-based genotyping must take BAM or CRAM alignment files as input. Exit now...")
			else:
				region = config["region"]
				if inFileType == "BAM":
					shell("~hsiehph/bin/pangenie/PanGenie  -k %s  -i <(samtools view -h %s --threads {threads} %s | samtools fastq --threads {threads} -)  -r {input.reference}  -v <(zcat {input.pangenome})  -s {wildcards.SAMPLE}  -o  output/{wildcards.SAMPLE}  -t {threads}  -j {threads}" % (kmer_size, shortread, region))
				elif inFileType == "CRAM":
					shell("~hsiehph/bin/pangenie/PanGenie  -k %s  -i <(samtools view -h --reference {input.reference_shortRead} --threads {threads} %s %s | samtools fastq --threads {threads} -)  -r {input.reference}  -v <(zcat {input.pangenome})  -s {wildcards.SAMPLE}  -o  output/{wildcards.SAMPLE}  -t {threads}  -j {threads}" % (kmer_size, shortread, region))

		else:
			print("genotype the whole genome")
			if inFileType == "FASTQ":
				shell("~hsiehph/bin/pangenie/PanGenie  -k %s  -i %s  -r {input.reference}  -v <(zcat {input.pangenome})  -s {wildcards.SAMPLE}  -o  output/{wildcards.SAMPLE}  -t {threads}  -j {threads} " % (kmer_size, shortread))
			elif inFileType == "BAM":
				shell("~hsiehph/bin/pangenie/PanGenie  -k %s  -i <(samtools fastq --threads {threads} %s)  -r {input.reference}  -v <(zcat {input.pangenome})  -s {wildcards.SAMPLE}  -o  output/{wildcards.SAMPLE}  -t {threads}  -j {threads}" % (kmer_size, shortread))
			elif inFileType == "CRAM":
				shell("~hsiehph/bin/pangenie/PanGenie  -k %s  -i <(samtools fastq --reference {input.reference_shortRead} --threads {threads} %s)  -r {input.reference}  -v <(zcat {input.pangenome})  -s {wildcards.SAMPLE}  -o  output/{wildcards.SAMPLE}  -t {threads}  -j {threads}" % (kmer_size, shortread))



rule bgzip:
	input:
		'output/{SAMPLE}_genotyping.vcf',
		'output/{SAMPLE}_histogram.histo',
		'output/{SAMPLE}_path_segments.fasta',
	output:
		'output/{SAMPLE}_genotyping.vcf.gz'
	resources:
		mem=10,
		hrs=48,
	shell:
		"""
		cat {input} | bgzip -c > {output} && tabix -p vcf {output}

		"""


