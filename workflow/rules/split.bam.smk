# workflow/rules/split.bam.smk

rule split_single_bam:
    '''
    Split a single BAM file (pool) by celltype and donor.
    Runs once per unique BAM file.
    '''
    input:
        rules.filter_metadata.output.filtered_metadata
    output:
        config['output_folder'] + '/split-bams/{pool_id}.done'
    params:
        celltype_col=config['celltype_column'],
        output_bams=config['output_folder'] + '/split-bams',
        pool_id='{pool_id}',
        output_metas=config['output_folder'] + '/split-metadatas',
    threads: 1
    resources:
        mem_mb=2000,
        time="02:00:00"
    params:
        scar=config['scar_dir'] + '/build/scar'
    shell:
        '''
        # Create per-pool metadata
        mkdir -p {params.output_metas}
        head -1 {input} > {params.output_metas}/{params.pool_id}.metadata.filtered.tsv
        grep -w "{params.pool_id}" {input} >> {params.output_metas}/{params.pool_id}.metadata.filtered.tsv || true
        
        # Run scar split
        {params.scar} \
        split \
        -i {params.output_metas}/{params.pool_id}.metadata.filtered.tsv \
        -o {params.output_bams} \
        -c {params.celltype_col}
        
        touch {output}
        '''


def get_pool_ids_from_metadata(wildcards):
    '''
    Extract unique pool_ids from the original metadata (not filtered).
    This runs at DAG construction time using the INPUT metadata.
    '''
    # Use the ORIGINAL metadata file which exists
    original_metadata = pd.read_table(config['metadata_file'])
    return original_metadata["pool_id"].unique().tolist()


checkpoint aggregate_split_bams:
    '''
    Checkpoint that waits for all BAM splits to complete.
    Enables dynamic discovery of celltype/donor combinations.
    '''
    input:
        filtered_metadata=rules.filter_metadata.output.filtered_metadata,
        flags=expand(
            config['output_folder'] + '/split-bams/{pool_id}.done',
            pool_id=get_pool_ids_from_metadata(None)  # Use original metadata
        )
    output:
        done=touch(config['output_folder'] + '/split-bams/all.done')
    threads: 1
    shell:
        '''
        echo "All BAM splits completed"
        echo "Split BAMs location: {config[output_folder]}/split-bams/"
        '''


rule flagstat_split_bam:
    '''
    Run samtools flagstat on each split BAM file to collect initial QC metrics.
    '''
    input:
        bam=config['output_folder'] + '/split-bams/{celltype}/{donor}.bam'
    output:
        flagstat=config['output_folder'] + '/qc/flagstat-split-bams/{celltype}/{donor}.flagstat.txt',
        qc_tsv=config['output_folder'] + '/qc/flagstat-split-bams/{celltype}/{donor}.qc.tsv'
    threads: 1
    resources:
        mem_mb=1000,
        time="00:05:00"
    shell:
        '''
        mkdir -p $(dirname {output.flagstat})
        
        # Run flagstat
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        
        # Parse flagstat to TSV format
        echo -e "sample\\tcelltype\\tstep\\treads\\tpairs\\tmapped\\tproperly_paired\\tsingletons\\tread1\\tread2" > {output.qc_tsv}
        
        reads=$(grep "in total" {output.flagstat} | awk '{{print $1}}')
        pairs=$(grep "paired in sequencing" {output.flagstat} | awk '{{print $1}}')
        mapped=$(grep "mapped (" {output.flagstat} | head -1 | awk '{{print $1}}')
        proper=$(grep "properly paired" {output.flagstat} | awk '{{print $1}}')
        singletons=$(grep "singletons" {output.flagstat} | awk '{{print $1}}')
        read1=$(grep "read1" {output.flagstat} | awk '{{print $1}}')
        read2=$(grep "read2" {output.flagstat} | awk '{{print $1}}')
        
        echo -e "{wildcards.donor}\\t{wildcards.celltype}\\t0_split\\t$reads\\t$pairs\\t$mapped\\t$proper\\t$singletons\\t$read1\\t$read2" >> {output.qc_tsv}
        '''


def get_all_flagstat_outputs(wildcards):
    '''
    Get all flagstat outputs after BAM splitting checkpoint.
    This reads the FILTERED metadata after the checkpoint completes.
    '''
    # Ensure checkpoint has run
    checkpoints.aggregate_split_bams.get(**wildcards).output.done
    
    # NOW read the filtered metadata (it exists because checkpoint completed)
    filtered_df = pd.read_table(config['output_folder'] + '/metadata-filtered/metadata.filtered.tsv')
    
    # Get unique celltype/donor combinations
    combinations = filtered_df[[config["celltype_column"], "donor_id"]].drop_duplicates()
    
    # Build QC paths
    qc_paths = []
    for _, row in combinations.iterrows():
        path = f"{config['output_folder']}/qc/flagstat-split-bams/{row[config['celltype_column']]}/{row['donor_id']}.qc.tsv"
        qc_paths.append(path)
    
    return qc_paths


rule aggregate_flagstats:
    '''
    Aggregate all initial flagstat QC metrics into a single report.
    '''
    input:
        get_all_flagstat_outputs
    output:
        config['output_folder'] + '/qc/flagstat_split_bams_qc.tsv'
    threads: 1
    shell:
        '''
        # Combine all QC files
        head -1 {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        echo "Aggregated QC for $(wc -l < {output}) samples"
        '''
