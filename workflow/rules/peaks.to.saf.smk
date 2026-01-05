rule bed_to_saf:
    '''
    Convert BED peak files to SAF format for featureCounts
    '''
    input:
        rules.filter_metadata.output.bed_files
    output:
        expand(config['output_folder'] + '/cell-peaks-saf/{celltype}_peaks.saf',celltype=metadata_df[config["celltype_column"]].unique())
    threads:
        1
    params:
        bed_dir=config['output_folder'] + '/cell-peaks-saf'
    resources:
        mem_mb=1000,
        time="0:02:00"
    shell:
        '''
        mkdir -p {params.bed_dir}
        for bed in {input}; do
            saf=$(basename $bed | sed 's/.bed/.saf/')
            awk 'OFS="\\t" {{print $1"_"$2"_"$3, $1, $2, $3, "."}}' \
            $bed > {params.bed_dir}/$saf
        done
        '''
