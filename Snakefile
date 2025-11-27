import pandas as pd
import os
import numpy as np
import sys
import re
from snakemake.io import glob_wildcards

configfile: "config/config.yaml"

DS = config.get("DS", ["10"])
DS_METHOD = config.get("DS_METHOD",["random"])
GENOME_SIZE = config.get("GENOME_SIZE", 3.1e9)
MANIFEST = config.get("MANIFEST", "config/manifest.tab")
QUALITY_THRESHOLD = config.get("QUALITY_THRESHOLD", ["10"])
read_min_length_dict = config["READ_MIN_LENGTH"]
CHUNK_COV = config.get("CHUNK_COV", 5)

wildcard_constraints:
    method = "|".join(["random","long"])

def get_read_paths(fofn):
    with open(fofn) as finp:
        read_paths = finp.read().strip().split("\n")
    return read_paths

def get_run_id(read_path):
    ont_runid_pattern = r"\d{8}_[0-9]{4}_[A-Z0-9]{2}_[A-Z]+[0-9]+_[a-z0-9]{8}"
 
    try:
        return re.findall(ont_runid_pattern, read_path)[0]
    except:
        return read_path.split("/")[-1].replace(".fastq.gz","").replace(".","_")


def expand_to_rows(df, column_name):
    rows = []
    for _, row in df.iterrows():
        values = row[column_name]
        for value in values:
            new_row = row.copy()
            new_row[column_name] = value
            rows.append(new_row)
    return pd.DataFrame(rows)

def find_raw_fastq(wildcards):
    return trim_df.at[wildcards.run_id, "READ_PATH"]


def find_trimmed_fastq(wildcards):
    run_ids = ds_df.loc[(wildcards.sample, wildcards.seq_type), "RUN_ID"]
    if isinstance(run_ids, pd.Series):
        run_ids = run_ids.tolist()
    else:
        run_ids = [run_ids]
    return [ f"quality_trimmed_reads/fastq/{wildcards.sample}.{wildcards.seq_type}.{run_id}.Q{wildcards.quality_threshold}.fastq.gz" for run_id in run_ids ]

def get_window_size(wildcards):
    return read_min_length_dict[wildcards.seq_type]

def chunk_fastqs_for_downsampled_reads(wildcards):
    check_point = checkpoints.sample_reads.get(
        sample=wildcards.sample,
        seq_type=wildcards.seq_type,
        method=wildcards.method,
        quality_threshold=wildcards.quality_threshold,
        ds=wildcards.ds,
    )
    chunk_dir = check_point.output.regions  # downsampled_reads/tmp/{sample}/{...}_{ds}

    pattern = os.path.join(chunk_dir, "{chunk}.tab")
    chunks = glob_wildcards(pattern).chunk

    # return for each chunk fastq
    return expand(
        "downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X_{chunk}.fastq.gz",
        sample=wildcards.sample,
        seq_type=wildcards.seq_type,
        method=wildcards.method,
        quality_threshold=wildcards.quality_threshold,
        ds=wildcards.ds,
        chunk=chunks,
    )

def get_chunk_tabs(wildcards):
    check_point = checkpoints.sample_reads.get(
        sample=wildcards.sample,
        seq_type=wildcards.seq_type,
        method=wildcards.method,
        quality_threshold=wildcards.quality_threshold,
        ds=wildcards.ds,
    )
    chunk_dir = check_point.output.regions
    pattern = os.path.join(chunk_dir, "{chunk}.tab")
    chunks = glob_wildcards(pattern).chunk

    if not chunks:
        return []

    return expand(
        "downsampled_reads/tmp/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X/{chunk}.tab",
        sample=wildcards.sample,
        seq_type=wildcards.seq_type,
        method=wildcards.method,
        quality_threshold=wildcards.quality_threshold,
        ds=wildcards.ds,
        chunk=chunks,
    )

manifest_df = pd.read_csv(MANIFEST, sep="\t", comment="#")
manifest_df["READ_PATH"] = manifest_df["FOFN"].apply(get_read_paths)
manifest_df = expand_to_rows(manifest_df, "READ_PATH").reset_index(drop=True)
manifest_df["RUN_ID"] = manifest_df["READ_PATH"].apply(get_run_id).reset_index(drop=True)
trim_df = manifest_df.copy().set_index(["RUN_ID"], drop=True).sort_index()
ds_df = manifest_df.copy().set_index(["SAMPLE", "SEQ_TYPE"], drop=True).sort_index()


rule all:
    input:
        expand(
            "downsampled_reads/stats/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X.stats",
            sample = ds_df.index.get_level_values("SAMPLE").unique(),
            seq_type = ds_df.index.get_level_values("SEQ_TYPE").unique(),
            ds = DS,
            method = DS_METHOD,
            quality_threshold = QUALITY_THRESHOLD,
            )

rule trim_reads:
    input:
        fastq = find_raw_fastq
    output:
        fastq = "quality_trimmed_reads/fastq/{sample}.{seq_type}.{run_id}.Q{quality_threshold}.fastq.gz",
        fai = "quality_trimmed_reads/fastq/{sample}.{seq_type}.{run_id}.Q{quality_threshold}.fastq.gz.fai",
        report = "quality_trimmed_reads/report/{sample}.{seq_type}.{run_id}.Q{quality_threshold}.report",
        html = "quality_trimmed_reads/report/html/{sample}.{seq_type}.{run_id}.Q{quality_threshold}.html",
        json = "quality_trimmed_reads/report/json/{sample}.{seq_type}.{run_id}.Q{quality_threshold}.json",
    threads: 8
    resources:
        mem = 8,
        hrs = 96,
    params:
        window_size = get_window_size,
    singularity:
        "docker://eichlerlab/fastqtrim:0.1"
    shell: """
        fastplong -i {input.fastq} -5 -3 -M {wildcards.quality_threshold} -q {wildcards.quality_threshold} -w {threads} -l {params.window_size} \
        -R {wildcards.sample}.{wildcards.seq_type}.{wildcards.run_id}.Q{wildcards.quality_threshold} -h {output.html} -j {output.json} --stdout 2> {output.report} | bgzip -c > {output.fastq}
        samtools fqidx {output.fastq}
        """

checkpoint sample_reads:
    input:
        trimmed_fastq_files = find_trimmed_fastq
    output:
        regions = directory("downsampled_reads/tmp/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X"),
        flag = "downsampled_reads/tmp/flags/{sample}.{seq_type}_{method}_Q{quality_threshold}_{ds}X.done"
    threads: 1
    resources:
        mem = 24,
        hrs = 72
    params:
        chunk_cov = CHUNK_COV
    run:
        outdir = output.regions
        outflag = output.flag

        os.makedirs(outdir, exist_ok=True)

        df = pd.DataFrame()        
        for fastq_file in input.trimmed_fastq_files:
            fai_df = pd.read_csv(fastq_file+".fai", sep="\t", header=None, usecols=[0,1], names=["read_name", "len"])
            fai_df["source"] = fastq_file.rstrip()
            df = pd.concat([df, fai_df], ignore_index=True)

        exp_cov = float(eval(wildcards.ds)) * GENOME_SIZE
        if wildcards.method == "long":
            cov_df = df.sort_values("len", ascending=False).reset_index(drop=True)
        else:
            cov_df = df.sample(frac=1, random_state=42).reset_index(drop=True) ## 42 = Arbitrary random seed.    
        cov_df["cum_cov"] = np.cumsum(cov_df["len"])
        cov_df = cov_df.loc[cov_df["cum_cov"] <= exp_cov ].copy() # until the total coverage for DS.
        
        chunk_cov_bases = float(params.chunk_cov) * GENOME_SIZE

        cov_df["chunk"] = ((cov_df["cum_cov"] - 1) // chunk_cov_bases).astype(int)

        for chunk_id, sub_df in cov_df.groupby("chunk"):
            out_path = os.path.join(outdir, f"{chunk_id}.tab")
            sub_df.to_csv(out_path, sep="\t", index=False)

        with open(outflag, "w") as f_flag:
            print (f"{wildcards.sample}.{wildcards.seq_type}_{wildcards.method}_Q{wildcards.quality_threshold}_{wildcards.ds}", file=f_flag)


rule extract_chunk_reads:
    input:
        regions = "downsampled_reads/tmp/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X/{chunk}.tab"
    output:
        reads = temp("downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X_{chunk}.fastq")
    params:
        reg = lambda wildcards, resources: f"{resources.tmpdir}/{wildcards.sample}/{wildcards.seq_type}_{wildcards.method}_Q{wildcards.quality_threshold}_{wildcards.ds}_{wildcards.chunk}_reg.tab",
    threads: 1
    resources:
        mem = 96,
        hrs = 120,
        heavy_io = 1
    run:
        ## in case the temp fastq was not deleted due to an interruption.
        if os.path.isfile(f"{resources.tmpdir}/{wildcards.sample}/{os.path.basename(output.reads)}"): 
            out_read_base = os.path.basename(output.reads)
            print(f"delete tmp : {out_read_base}")
            shell(f"rm -f {resources.tmpdir}/{wildcards.sample}/$(basename {output.reads})")
        ## <===============================

        df = pd.read_csv(input.regions, sep="\t")
        
        os.makedirs(f"{resources.tmpdir}/{wildcards.sample}", exist_ok=True)
        

        for file in df["source"].unique():
            file_base = os.path.basename(file) # raw fastq source name
            reg_df = df.loc[df["source"] == file]
            reg_df[["read_name"]].to_csv(params.reg, sep="\t", header=False, index=False)
            shell(
                """. /usr/share/Modules/init/bash;"""
                """module load modules modules-init modules-gs/prod modules-eichler/prod;"""
                """module load seqtk/1.4 samtools/1.19;"""
                """samtools fqidx -r {params.reg} {file} | seqtk seq -l0 >> {resources.tmpdir}/{wildcards.sample}/$(basename {output.reads})"""
                )

#            shell(f"samtools fqidx -r {params.reg} {file} | seqtk seq -l0 >> {resources.tmpdir}/{wildcards.sample}/$(basename {output.reads})")
        shell(f"rsync -av {resources.tmpdir}/{wildcards.sample}/$(basename {output.reads}) {output.reads}")
        shell(f"rm -f {params.reg}")
        shell(f"rm -f {resources.tmpdir}/{wildcards.sample}/$(basename {output.reads})")

rule compress_and_index:
    input:
        reads = rules.extract_chunk_reads.output.reads
    output:
        reads = temp("downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X_{chunk}.fastq.gz"),
        fai = "downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X_{chunk}.fastq.gz.fai"
    threads: 16,
    resources:
        mem = 16,
        hrs = 72
    shell: """
        bgzip -@ {threads} -c {input.reads} > {output.reads}
        samtools fqidx {output.reads}
        """

rule merged_and_index:
    input:
        chunk_fastqs_for_downsampled_reads
    output:
        reads = "downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X.fastq.gz",
        fai = "downsampled_reads/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X.fastq.gz.fai",
    threads:1,
    resources:
        mem = 128,
        hrs = 120,
    run:
        chunk_files = list(input)
        if len(chunk_files) == 1:
            shell(f"cp {chunk_files[0]} {output.reads}")
        else:
            shell(f"cat {' '.join(chunk_files)} > {output.reads}")
        shell(f"samtools fqidx {output.reads}")

rule calculate_stats:
    input:
        reads = rules.merged_and_index.output.reads,
        fai = rules.merged_and_index.output.fai
    output:
        stats = "downsampled_reads/stats/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X.stats",
        plot = "downsampled_reads/stats/{sample}/{seq_type}_{method}_Q{quality_threshold}_{ds}X.logged.len_dist.png"
    threads: 1,
    resources:
        mem = 16,
        hrs = 4
    shell: """
        /net/eichler/vol28/software/pipelines/compteam_tools/ont_stats -f {input.reads} -s {wildcards.sample}_Q{wildcards.quality_threshold}_{wildcards.ds}X -o {output.stats} -p {output.plot} -l
        """
        
