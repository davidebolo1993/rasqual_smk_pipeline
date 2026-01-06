rule filter_metadata:
    '''
    Filter metadata to retain cell types with sufficient representation
    Generate per-cell-type peak BED files.
    #Need to add singularity container in the end
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
        mem_mb=1000,
        time="0:01:00"
    params:
        celltype_col=config['celltype_column'],
        min_cells=config.get('min_cells', 10),
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
