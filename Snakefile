import pandas as pd
import os
import numpy as np
import sys
import re
from collections import defaultdict

configfile: "config/config.yaml"

DS = config.get("DS", ["10"])
DS_METHOD = config.get("DS_METHOD",["random"])
GENOME_SIZE = config.get("GENOME_SIZE", 3.1e9)
MANIFEST = config.get("MANIFEST", "config/manifest.tab")
QUALITY_THRESHOLD = config.get("QUALITY_THRESHOLD", ["20"])
READ_MIN_LENGTH = config.get("READ_MIN_LENGTH",21) # meryk k=21


wildcard_constraints:
    __default__ = r"[^.]+",
    method = "|".join(["random","long"]),
    ds = r"\d+"


def get_read_paths(fofn):
    read_paths = []
    with open(fofn) as finp:
        read_paths = finp.read().strip().split("\n")
    return read_paths


def expand_to_rows_paired(df, column_names):
    rows = []
    for _, row in df.iterrows():
        col1 = row[column_names[0]]
        col2 = row[column_names[1]]
        if not isinstance(col1, (list, tuple)):
            col1 = [col1]
        if not isinstance(col2, (list, tuple)):
            col2 = [col2]
        for value1, value2 in zip(col1, col2):
            new_row = row.copy()
            new_row[column_names[0]] = value1
            new_row[column_names[1]] = value2
            rows.append(new_row)
    return pd.DataFrame(rows)


def find_raw_fastq(which_one):
    def inner(wildcards):
        return trim_df.at[wildcards.basename, f"{which_one}_READ_PATH"]
    return inner


def find_trimmed_fastq_files(which_one):
    def inner(wildcards):
        basenames = ds_df.loc[wildcards.sample, "BASENAME"]
        if isinstance(basenames, pd.Series):
            basenames = basenames.tolist()
        else:
            basenames = [basenames]
        return [ f"quality_trimmed_reads/fastq/{wildcards.sample}.PE.{basename}.Q{wildcards.quality_threshold}.{which_one}.fastq.gz" for basename in basenames ]
    return inner


def find_read_ids(which_one):
    def inner(wildcards):
        basenames = ds_df.loc[wildcards.sample, "BASENAME"]
        if isinstance(basenames, pd.Series):
            basenames = basenames.tolist()
        else:
            basenames = [basenames]
        read_dir = f"downsampled_reads/tmp/{wildcards.sample}/PE_{wildcards.method}_Q{wildcards.quality_threshold}_{wildcards.ds}_{which_one}_lists"
        return [ f"{read_dir}/{wildcards.sample}.PE.{basename}.Q{wildcards.quality_threshold}.{which_one}.reads.txt" for basename in basenames]
    return inner
    

def get_basename(read_path):
    return "_".join(os.path.basename(read_path).replace(".fastq.gz","").replace(".","_").split("_")[:-1])
        
    
df = pd.read_csv(MANIFEST, sep="\t", comment="#")
df["R1_READ_PATH"] = df["R1_FOFN"].apply(get_read_paths)
df["R2_READ_PATH"] = df["R2_FOFN"].apply(get_read_paths)
df = expand_to_rows_paired(df,["R1_READ_PATH","R2_READ_PATH"])
df["BASENAME"] = df["R1_READ_PATH"].apply(get_basename)
df = df.drop_duplicates().reset_index(drop=True)

trim_df = df.copy().set_index(["BASENAME"], drop=True).sort_index()
ds_df = df.copy().set_index(["SAMPLE"], drop=True).sort_index()


rule all:
    input:
        expand("downsampled_reads/stats/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_{read}.stats",
            sample = ds_df.index.get_level_values('SAMPLE').unique(),
            ds = DS,
            method = DS_METHOD,
            quality_threshold = QUALITY_THRESHOLD,
            read = ["R1","R2"],
        ),
rule trim_reads:
    input:
        r1_raw_fastq = find_raw_fastq(which_one="R1"),
        r2_raw_fastq = find_raw_fastq(which_one="R2")
    output:
        r1_trim_fastq = temp("quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R1.fastq"),
        r2_trim_fastq = temp("quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R2.fastq"),
        html = "quality_trimmed_reads/report/html/{sample}.PE.{basename}.Q{quality_threshold}.html",
        json = "quality_trimmed_reads/report/json/{sample}.PE.{basename}.Q{quality_threshold}.json",
    params:
        read_min_length = READ_MIN_LENGTH,
        dup_calc_accuracy = 6,
    threads: 16
    resources:
        mem = 8,
        hrs = 96,
    singularity:
        "docker://eichlerlab/fastqtrim:0.2"
    shell: """
        fastp -w {threads} \
        -i {input.r1_raw_fastq} -I {input.r2_raw_fastq} \
        -o {output.r1_trim_fastq} -O {output.r2_trim_fastq} \
        -l {params.read_min_length} \
        --dedup --dup_calc_accuracy {params.dup_calc_accuracy} \
        --html {output.html} --json {output.json} \
        --report_title "{wildcards.sample}_{wildcards.basename}"
        """

rule compress_r1:
    input:
        r1_trim_fastq = rules.trim_reads.output.r1_trim_fastq
    output:
        r1_trim_fastq_gz = "quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R1.fastq.gz",
        r1_trim_fai = "quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R1.fastq.gz.fai"
    threads: 16
    resources:
        mem = 8,
        hrs = 96,
    shell: """
        bgzip -@ {threads} -c {input.r1_trim_fastq} > {output.r1_trim_fastq_gz}
        samtools fqidx {output.r1_trim_fastq_gz}
        """

rule compress_r2:
    input:
        r2_trim_fastq = rules.trim_reads.output.r2_trim_fastq
    output:
        r2_trim_fastq_gz = "quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R2.fastq.gz",
        r2_trim_fai = "quality_trimmed_reads/fastq/{sample}.PE.{basename}.Q{quality_threshold}.R2.fastq.gz.fai"
    threads: 16
    resources:
        mem = 8,
        hrs = 96,
    shell: """
        bgzip -@ {threads} -c {input.r2_trim_fastq} > {output.r2_trim_fastq_gz}
        samtools fqidx {output.r2_trim_fastq_gz}
        """

rule build_pairs_index:
    input:
        trim_r1_fastq_files = find_trimmed_fastq_files(which_one="R1"),
        trim_r2_fastq_files = find_trimmed_fastq_files(which_one="R2")
    output:
        pairs_index = "downsampled_reads/pairs_index/{sample}/PE_Q{quality_threshold}.tsv.gz",
        flag = "downsampled_reads/pairs_index/{sample}/.PE_Q{quality_threshold}.done"
    threads: 1
    resources:
        mem = 320,
        hrs = 96,
        heavy_io = 3,        
    benchmark: "downsampled_reads/pairs_index/{sample}/PE_Q{quality_threshold}.benchmark.txt"
    run:
        def load_fai_table(fq):
            fai = fq + ".fai"
            df = pd.read_csv(fai, sep="\t", header=None, usecols=[0,1], names=["read_name", "len"], dtype={"read_name":"string","len":"int32"})
            df["source"] = pd.Categorical([os.path.basename(fq).replace(".R1.fastq.gz","").replace(".R2.fastq.gz","")]*len(df))
            return df

        print ("Loading fai files...")
        r1_df_all = pd.concat([load_fai_table(fq) for fq in input.trim_r1_fastq_files], ignore_index=True)
        r2_df_all = pd.concat([load_fai_table(fq) for fq in input.trim_r2_fastq_files], ignore_index=True)

        r1_df_all["end"] = "R1"
        r2_df_all["end"] = "R2"
        key_cols = ["read_name","source"]
        merged = r1_df_all.merge(r2_df_all, on=key_cols, suffixes=("_r1","_r2"))
        del r1_df_all, r2_df_all

        merged["pair_len"] = (merged["len_r1"].astype("int64") + merged["len_r2"].astype("int64"))
        merged = merged[["read_name","source","pair_len"]]

        print(f"Saving merged pairs_index to {output.pairs_index}")
        merged.to_csv(output.pairs_index, sep="\t", index=False, compression="gzip")
        with open(output.flag,"w") as fout:
            fout.close()

rule sample_reads_r1:
    input:
        pairs_index = rules.build_pairs_index.output.pairs_index,
        flag = rules.build_pairs_index.output.flag,
    output:
        selected_list = directory("downsampled_reads/tmp/{sample}/PE_{method}_Q{quality_threshold}_{ds}_R1_lists"),
        flag = "downsampled_reads/tmp/{sample}/.PE_{method}_Q{quality_threshold}_{ds}.R1.sample_done"
    threads: 1
    resources:
        mem = 120,
        hrs = 96,
        heavy_io = 3,
    benchmark: "downsampled_reads/tmp/{sample}/.PE_{method}_Q{quality_threshold}_{ds}.R1.benchmark.txt"
    run:
        os.makedirs(output.selected_list, exist_ok=True)
        merged = pd.read_csv(input.pairs_index, sep="\t")

        if wildcards.method == "long":
            cov_df = merged.sort_values("pair_len", ascending=False).reset_index(drop=True)
        else:  # random
            cov_df = merged.sample(frac=1, random_state=42).reset_index(drop=True)

        exp_cov = float(eval(wildcards.ds)) * float(GENOME_SIZE)

        cov_csum = np.cumsum(cov_df["pair_len"])
        keep_idx = cov_csum <= exp_cov
        out_df = cov_df.loc[keep_idx, ["read_name","source"]]

        grouped_by_source = {src: sub for src, sub in out_df.groupby("source", sort=False)}

        for src_name, sub_df in grouped_by_source.items():
            list_path = os.path.join(output.selected_list, f"{src_name}.R1.reads.txt")
            sub_df[["read_name"]].to_csv(list_path, sep="\t", header=False, index=False)
        with open(output.flag,"w") as fout:
            fout.close()

rule sample_reads_r2:
    input:
        pairs_index = rules.build_pairs_index.output.pairs_index,
        flag = rules.build_pairs_index.output.flag,
    output:
        selected_list = directory("downsampled_reads/tmp/{sample}/PE_{method}_Q{quality_threshold}_{ds}_R2_lists"),
        flag = "downsampled_reads/tmp/{sample}/.PE_{method}_Q{quality_threshold}_{ds}.R2.sample_done"
    threads: 1
    resources:
        mem = 120,
        hrs = 96,
        heavy_io = 3,
    benchmark: "downsampled_reads/tmp/{sample}/.PE_{method}_Q{quality_threshold}_{ds}.R2.benchmark.txt"
    run:
        os.makedirs(output.selected_list, exist_ok=True)
        merged = pd.read_csv(input.pairs_index, sep="\t")

        if wildcards.method == "long":
            cov_df = merged.sort_values("pair_len", ascending=False).reset_index(drop=True)
        else:  # random
            cov_df = merged.sample(frac=1, random_state=42).reset_index(drop=True)

        exp_cov = float(eval(wildcards.ds)) * float(GENOME_SIZE)

        cov_csum = np.cumsum(cov_df["pair_len"])
        keep_idx = cov_csum <= exp_cov
        out_df = cov_df.loc[keep_idx, ["read_name","source"]]

        grouped_by_source = {src: sub for src, sub in out_df.groupby("source", sort=False)}

        for src_name, sub_df in grouped_by_source.items():
            list_path = os.path.join(output.selected_list, f"{src_name}.R2.reads.txt")
            sub_df[["read_name"]].to_csv(list_path, sep="\t", header=False, index=False)
        with open(output.flag,"w") as fout:
            fout.close()
        

rule extract_reads_r1:
    input:
        flag = rules.sample_reads_r1.output.flag
    output:
        reads = temp('downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.fastq')
    params:
        regions = find_read_ids(which_one="R1"),
        reg = '/data/scratch/tmp/{sample}/PE_{method}_Q{quality_threshold}_{ds}_R1_reg.tab',
    threads: 1
    resources:
        mem = 96,
        hrs = 72
    run:
        ## in case the temp fastq was not deleted due to an interruption.
        if os.path.isfile(f"{resources.tmpdir}/{os.path.basename(output.reads)}"): 
            out_read_base = os.path.basename(output.reads)
            print(f"delete tmp : {out_read_base}")
            shell(f'rm -f {resources.tmpdir}/$(basename {output.reads})')
        ## <===============================

        os.makedirs(f"/data/scratch/tmp/{wildcards.sample}", exist_ok=True)
        for region in params.regions:
            file_base = os.path.basename(region).replace(".reads.txt","")
            token = file_base.split(".")
            print (token)
            sample, seqtype, basename, quality, read = token
            source = f"quality_trimmed_reads/fastq/{sample}.PE.{basename}.{quality}.{read}.fastq.gz"
            reg_df = pd.read_csv(region, sep='\t')
            reg_df.to_csv(params.reg, sep='\t', header=False, index=False)
            shell(
                """. /usr/share/Modules/init/bash;"""
                """module load modules modules-init modules-gs/prod modules-eichler/prod;"""
                """module load seqtk/1.4 samtools/1.19;"""
                """samtools fqidx -r {params.reg} {source} | seqtk seq -l0 >> {resources.tmpdir}/$(basename {output.reads})"""
                )
#            shell(f'samtools fqidx -r {params.reg} {file} | seqtk seq -l0 >> {resources.tmpdir}/$(basename {output.reads})')
            shell(f'rsync -av {resources.tmpdir}/$(basename {output.reads}) {output.reads}')
            shell(f'rm -f {params.reg}')
            shell(f'rm -f {resources.tmpdir}/$(basename {output.reads})')


rule extract_reads_r2:
    input:
        flag = rules.sample_reads_r2.output.flag
    output:
        reads = temp('downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.fastq')
    params:
        regions = find_read_ids(which_one="R2"),
        reg = '/data/scratch/tmp/{sample}/PE_{method}_Q{quality_threshold}_{ds}_R2_reg.tab',
    threads: 1
    resources:
        mem = 96,
        hrs = 72
    run:
        ## in case the temp fastq was not deleted due to an interruption.
        if os.path.isfile(f"{resources.tmpdir}/{os.path.basename(output.reads)}"): 
            out_read_base = os.path.basename(output.reads)
            print(f"delete tmp : {out_read_base}")
            shell(f'rm -f {resources.tmpdir}/$(basename {output.reads})')
        ## <===============================

        os.makedirs(f"/data/scratch/tmp/{wildcards.sample}", exist_ok=True)
        for region in params.regions:
            file_base = os.path.basename(region).replace(".reads.txt","")
            token = file_base.split(".")
            print (token)
            sample, seqtype, basename, quality, read = token
            source = f"quality_trimmed_reads/fastq/{sample}.PE.{basename}.{quality}.{read}.fastq.gz"
            reg_df = pd.read_csv(region, sep='\t')
            reg_df.to_csv(params.reg, sep='\t', header=False, index=False)
            shell(
                """. /usr/share/Modules/init/bash;"""
                """module load modules modules-init modules-gs/prod modules-eichler/prod;"""
                """module load seqtk/1.4 samtools/1.19;"""
                """samtools fqidx -r {params.reg} {source} | seqtk seq -l0 >> {resources.tmpdir}/$(basename {output.reads})"""
                )
#            shell(f'samtools fqidx -r {params.reg} {file} | seqtk seq -l0 >> {resources.tmpdir}/$(basename {output.reads})')
            shell(f'rsync -av {resources.tmpdir}/$(basename {output.reads}) {output.reads}')
            shell(f'rm -f {params.reg}')
            shell(f'rm -f {resources.tmpdir}/$(basename {output.reads})')


rule compress_and_index_r1:
    input:
        reads = rules.extract_reads_r1.output.reads
    output:
        reads = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.fastq.gz',
        fai = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.fastq.gz.fai',
    threads: 16
    resources:
        mem = 8,
        hrs = 72
    shell: """
        bgzip -@ {threads} -c {input.reads} > {output.reads}
        samtools fqidx {output.reads}
        """

rule compress_and_index_r2:
    input:
        reads = rules.extract_reads_r2.output.reads
    output:
        reads = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.fastq.gz',
        fai = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.fastq.gz.fai',
    threads: 16
    resources:
        mem = 8,
        hrs = 72
    shell: """
        bgzip -@ {threads} -c {input.reads} > {output.reads}
        samtools fqidx {output.reads}
        """


rule calculate_stats_r1:
    input:
        reads = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.fastq.gz',
    output:
        stats = 'downsampled_reads/stats/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.stats',
        plot = 'downsampled_reads/stats/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R1.logged.len_dist.png'
    threads: 1
    resources:
        mem = 64,
        hrs = 72
    shell: """
        /net/eichler/vol28/software/pipelines/compteam_tools/ont_stats -f {input.reads} -s {wildcards.sample}_Q{wildcards.quality_threshold}_{wildcards.ds}X_R1 -o {output.stats} -p {output.plot} -l -w 1 -x 151
        """
        
rule calculate_stats_r2:
    input:
        reads = 'downsampled_reads/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.fastq.gz',
    output:
        stats = 'downsampled_reads/stats/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.stats',
        plot = 'downsampled_reads/stats/{sample}/PE_{method}_Q{quality_threshold}_{ds}X_R2.logged.len_dist.png'
    threads: 1
    resources:
        mem = 64,
        hrs = 72
    shell: """
        /net/eichler/vol28/software/pipelines/compteam_tools/ont_stats -f {input.reads} -s {wildcards.sample}_Q{wildcards.quality_threshold}_{wildcards.ds}X_R2 -o {output.stats} -p {output.plot} -l -w 1 -x 151
        """


