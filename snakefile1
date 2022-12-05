configfile: "config.yaml"
DATA=config["DATA"]
print(DATA)
#
input_dir="data/"
#
log_dir="log/"
#
fastqc_dir="results/fastqc_results/"
#
fastp_dir="results/trim/fastp/"
trammomatic_dir="results/trim/trimmomatic/"
sickle_dir="results/trim/sickle/"
#
reference_genome="data/Genome/Gmax_275_v2.0.fa"
index_genome_dir="data/Genome/"
#
hisat_dir="results/align/hisat2/"
#
reference_gtf = "data/Genome/Gmax_275_Wm82.a2.v1.gene_exons.gff3"
stringtie_dir = "results/quant/stringtie/"

rule all:
    input:
        expand(stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf",data = DATA),
        stringtie_dir+"diff_gene_list.csv"

rule fastqc_fq:
    input:
        input_dir+"{data}.fq"
    output:
        multiext(fastqc_dir+ "{data}",".zip",".html")
    log:
        log_dir+"{data}.fastqc.log"
    params:
        fastqc_dir=fastqc_dir,
        prefix = "{data}"
    threads: 4
    shell:
        "fastqc {input} -o {params.fastqc_dir} 2> {log}"

rule multiqc_fq:
    input:
        fastqc_dir
    output:
        fastqc_dir+"multiqc_report.html"
    threads: 4
    shell:
        "multiqc {input}/. -o {input} "

rule fastp_trim:
    input:
        fq=input_dir+"{data}.fq"
    output:
        fq=fastp_dir+"{data}.trim_fastp.fq",
        html=fastp_dir+"{data}.trim_fastp.html",
        json=fastp_dir+"{data}.trim_fastp.json"
    log:
        log_dir+"{data}.trim_fastp.log"
    threads: 4
    shell:
        "fastp -i {input.fq} -o {output.fq} --thread {threads} --html {output.html} --json {output.json} 2>{log}"

rule trimmomatic_trim:
	input:
		fq=input_dir+"{data}.fq"
	output:
		trim_fq=trammomatic_dir+"{data}.trim_momatic.fq"
	shell:
		"""
        trimmomatic SE {input.fq} {output.trim_fq} 
        SLIDINGWINDOW:4:15 
        MINLEN:25 
        ILLUMINACLIP:TruSeq3-SE:2:30:10
        """

rule sickle_trim:
	input:
		fq=input_dir+"{data}.fq"
	output:
		trim_fq=sickle_dir+"{data}.trim_sickle.fq"
	shell:
		"""
        sickle se -t sanger -f {input} \
        -o {output} \
        -q 35 -l 45
        """

rule hisat2_index:
    input:
        reference_genome = reference_genome
    output:
        index_genome = index_genome_dir+"Hisat_index"
    threads:
        4
    shell:
        """
        hisat2-build -p 4 {input} {output}
        """

rule hisat2_align:
    input:
        fq = fastp_dir+"{data}.trim_fastp.fq"
    output:
        sam = hisat_dir + "{data}.align_hisat2.sam",
        bam = hisat_dir + "{data}.align_hisat2.bam"
    threads:
        4
    params:
        reference_genome = reference_genome,
        index_genome = rules.hisat2_index.output.index_genome
    shell:
        """
        hisat2 -p {threads} --dta -x {params.index_genome} -U {input.fq} -S {output.sam} --dta-cufflinks --no-unal
        samtools sort -@ {threads} -o {output.bam} {output.sam}
        """

rule stringtie_align:
    input:
        bam = rules.hisat2_align.output.bam
    output:
        gtf = stringtie_dir + "align/{data}.align_stringtie.gtf"
    params:
        l = "{data}",
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        stringtie -p {threads} -G {params.reference_gtf} -o {output.gtf} -l {params.l} {input.bam}
        """

rule stringtie_merge:
    input:
        gtf = expand(stringtie_dir + "align/{data}.align_stringtie.gtf",data=DATA)
    output:
        gtf = "results/quant/stringtie/stringtie_merged.gtf"
    params:
        l = "merge",
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        ls -1 {input.gtf} >> {params.stringtie_dir}/sample_gtf_list.txt
        stringtie --merge -p 4 -G {params.reference_gtf} -o {output.gtf}  {params.stringtie_dir}/sample_gtf_list.txt
        """

rule gffcompare:
    input:
        gtf = rules.stringtie_merge.output.gtf
    output:
        gffcompare = stringtie_dir+"gffcompare.annotated.gtf"
    params:
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        gffcompare -r {params.reference_gtf} -o {params.stringtie_dir}gffcompare {input.gtf}
        """

rule stringtie_quant:
    input:
        bam = rules.hisat2_align.output.bam,
        merged_gtf=rules.stringtie_merge.output.gtf
    output:
        gtf = stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf"
    threads:
        4
    shell:
        """
        stringtie -e -B -p {threads} -G {input.merged_gtf} -o {output.gtf} {input.bam}
        """

rule stringtie_count:
    input:
        gtf = expand(stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf",data=DATA)
    output:
        stringtie_dir+"gene_count_matrix.csv",
        stringtie_dir+"transcript_count_matrix.csv"
    params:
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        ls {input.gtf} > {stringtie_dir}gtf_path.txt
        # find {stringtie_dir} -name *.quant*.gtf > {stringtie_dir}gtf_path.txt
        cat {stringtie_dir}gtf_path.txt | cut -d / -f5 > {stringtie_dir}sampleid.txt
        paste {stringtie_dir}sampleid.txt {stringtie_dir}gtf_path.txt > {stringtie_dir}gtf_list.txt
        rm {stringtie_dir}sampleid.txt {stringtie_dir}gtf_path.txt

        python workflow/scripts/prepDE_py3.py \
        -i {stringtie_dir}gtf_list.txt \
        -t {stringtie_dir}transcript_count_matrix.csv \
        -g {stringtie_dir}gene_count_matrix.csv
        """

rule DEseq2:
    input:
        gene_count = stringtie_dir+"gene_count_matrix.csv",
        transcript_count = stringtie_dir+"transcript_count_matrix.csv"
    output:
        diff_gene_list = stringtie_dir+"diff_gene_list.csv"
    params:
        sample_group = input_dir+"sample_group.txt"
    #conda:
        #"envs/ballgown.yaml"
    script:
        "workflow/scripts/Deseq2.R"