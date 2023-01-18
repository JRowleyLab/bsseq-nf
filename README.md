# bsseq-nf
A Nextflow pipeline for the end-to-end data processing of Bisulfite Sequencing (BS-seq) paired-end data. 

## How to use

1. Clone the repository: `git clone <repo.git>`
2. Create a singularity image from the nf-core Docker container: `singularity pull bsseq.sif docker://nfcore/methylseq`
3. Create a samplesheet csv (`samples.csv`) that contains `sample`, `read 1` and `read 2` information like below:

**NOTE:** Use this exact format.

```
sample_id,fastq_1,fastq_2
Control_A,<path/to/read_1.fq.gz>,<path/to/read_2.fq.gz>
Control_B,<path/to/read_1.fq.gz>,<path/to/read_2.fq.gz>
```

4. Run the pipeline with the following code:

`nextflow run main.nf --samplesheet <samples.csv>`

## Replicates

With the `--merge` and `--plusmerge` arguments, replicates are merged. This is done by using the sample name and replicate separated by "_". (eg. **Sample_Replicate** OR **Control_A**)

## Parameters

* `--aligner`: `bwa-meth` OR `bismark` [default: `bwa-meth`]
* `--samplesheet`: Samplesheet csv containing `sample_id,fastq_1,fastq_2`.
* `--plusmerge`: Process all individual samples, AND combine replicates (eg. `Control_A`, `Control_B`).
* `--merge`: Combine replicates (eg. `Control_A`, `Control_B`).
* `--outdir`: Output directory
* `--index`: path to index [make sure it matches `--aligner`]
* `--genome`: path to genome, an index will be built based on the `--aligner` [WORK IN PROGRESS|MAY NOT WORK]

## Test

The pipeline can be tested with the following command `nextflow run main.nf --index <path/to/index>`. Test reads are found in the `data/` folder.

## Singularity

The singularity image generated can be used with `-with-singularity <path/to/image>` or by placing the path in the `nextflow.config` file.
