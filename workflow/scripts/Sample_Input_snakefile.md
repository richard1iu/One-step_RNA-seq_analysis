## snakefile

## method 1
import pandas

def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]

_samples = parse_samples("sample.txt")

## method 2
import glob
SAMPLES = []
for filename in glob.glob(r'sra_file/*.sra'):
    SAMPLES.append(filename.replace('.sra','').replace('sra_file/',''))
rule all:
    input:
        expand('fastq_file/{sample}.fastq', sample=SAMPLES)

rule all:
    input:
        expand("2.rmhost/{sample}_{R}.rmhost.fq.gz", sample=_samples.index, R=["1","2"])

rule trim:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        trim_r1 = "1.trimmed/{sample}_1.trimmed.fastq.gz",
        trim_r2 = "1.trimmed/{sample}_2.trimmed.fastq.gz"
    threads: 2
    params:
        min_len = 70
    log:
        json = "1.trimmed/{sample}.fastp.json",
        html = "1.trimmed/{sample}.fastp.html",
        fastp_log = "1.trimmed/{sample}.fastp.log"
    shell:
    """
    fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} -w {threads} --length_required {params.min_len} -j {log.json} -h {log.html} 2> {log.fastp_log}
    """
        