# workflow/rules/featurecounts.count.smk

def get_bams_for_celltype(wildcards):
    """
    Get all processed BAM files for a specific celltype.
    This aggregates across all donors within a celltype.
    """
    # Ensure checkpoint has run
    checkpoints.aggregate_split_bams.get(**wildcards).output.done
    
    # Read filtered metadata
    filtered_df = pd.read_table(config['output_folder'] + '/metadata-filtered/metadata.filtered.tsv')
    
    # Get donors for this celltype
    celltype_data = filtered_df[filtered_df[config["celltype_column"]] == wildcards.celltype]
    donors = celltype_data["donor_id"].unique()
    
    # Build BAM paths
    bam_paths = []
    for donor in donors:
        path = f"{config['output_folder']}/processed-bams/{wildcards.celltype}/{donor}.bam"
        bam_paths.append(path)
    
    return bam_paths


rule featurecounts_per_celltype:
    """
    Run featureCounts on all BAMs for a specific celltype.
    Aggregates read counts across all donors within the celltype.
    """
    input:
        saf=config['output_folder'] + '/cell-peaks-saf/{celltype}_peaks.saf',
        bams=get_bams_for_celltype
    output:
        counts=config['output_folder'] + '/counts/{celltype}.count_matrix',
        summary=config['output_folder'] + '/counts/{celltype}.count_matrix.summary'
    params:
        featurecounts=config['featurecounts_dir'] + '/featureCounts'
    threads: 10
    resources:
        mem_mb=10000,
        time="00:20:00"
    shell:
        '''
        mkdir -p $(dirname {output.counts})
                
        {params.featurecounts} \
            -p \
            -T {threads} \
            -F SAF \
            --donotsort \
            -a {input.saf} \
            -o {output.counts} \
            {input.bams}
        '''


def get_all_celltype_counts(wildcards):
    """
    Get count matrices for all celltypes.
    Called after split_bams checkpoint to dynamically determine celltypes.
    """
    # Ensure checkpoint has run
    checkpoints.aggregate_split_bams.get(**wildcards).output.done
    
    # Read filtered metadata to get actual celltypes that passed filtering
    filtered_df = pd.read_table(config['output_folder'] + '/metadata-filtered/metadata.filtered.tsv')
    celltypes = filtered_df[config["celltype_column"]].unique()
    
    # Build count matrix paths
    count_paths = []
    for celltype in celltypes:
        path = f"{config['output_folder']}/counts/{celltype}.count_matrix"
        count_paths.append(path)
    
    return count_paths


rule aggregate_counts:
    """
    Aggregate all featureCounts outputs.
    Creates a completion flag to track when all counting is done.
    """
    input:
        get_all_celltype_counts
    output:
        done=touch(config['output_folder'] + '/counts/all_counts.done'),
        summary=config['output_folder'] + '/counts/summary_report.txt'
    threads: 1
    shell:
        '''
        echo "FeatureCounts completed for all celltypes" > {output.summary}
        echo "Total count matrices: $(echo {input} | wc -w)" >> {output.summary}
        echo "" >> {output.summary}
        echo "Count matrices:" >> {output.summary}
        for matrix in {input}; do
            celltype=$(basename $matrix .count_matrix)
            n_features=$(tail -n +3 $matrix | wc -l)
            n_samples=$(head -1 $matrix | awk '{{print NF - 6}}')
            echo "  $celltype: $n_features features x $n_samples samples" >> {output.summary}
        done
        
        cat {output.summary}
        '''
