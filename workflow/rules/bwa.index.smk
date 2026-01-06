# workflow/rules/bwa.index.smk


rule bwa_mem2_index:
    '''
    Index the reference genome with bwa-mem2
    '''
    input:
        config['reference_genome']
    output:
        multiext(config['output_folder'] + '/wasp-processing/bwa_index/genome_index', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')
    threads: 1
    resources:
        mem_mb=100000,
        time="0:50:00"
    params:
        prefix=config['output_folder'] + '/wasp-processing/bwa_index/genome_index'
    shell:
        '''
        bwa-mem2 \
        index \
        -p {params.prefix} \
        {input}
        '''
