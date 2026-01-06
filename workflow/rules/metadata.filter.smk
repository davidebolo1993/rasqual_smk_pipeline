# workflow/rules/metadata.filter.smk

rule filter_metadata:
    '''
    Filter metadata to retain cell types with sufficient representation.
    Generate per-cell-type peak BED files.
    '''
    input:
        metadata=config['metadata_file'],
        peaks=config['peaks_file']
    output:
        filtered_metadata=config['output_folder'] + '/metadata-filtered/metadata.filtered.tsv',
        bed_files=expand(config['output_folder'] + '/cell-peaks-bed/{celltype}_peaks.bed',celltype=metadata_df[config["celltype_column"]].unique())
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        time_min=lambda wildcards, attempt: attempt * 1
    params:
        celltype_col=config['celltype_column'],
        min_cells=config['min_cells'],
        bed_dir=config['output_folder'] + '/cell-peaks-bed',
        scar=config['scar_dir'] + '/build/scar'
    shell:
        '''
        {params.scar} \
        filter \
        -i {input.metadata} \
        -o {output.filtered_metadata} \
        -c {params.celltype_col} \
        -m {params.min_cells} \
        -p {input.peaks} \
        -b {params.bed_dir}
        '''
