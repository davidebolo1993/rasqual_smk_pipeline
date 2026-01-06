# workflow/rules/bwa.index.smk

rule bwa_mem2_index:
    '''
    Index the reference genome with bwa-mem2.
    '''
    input:
        config['reference_genome']
    output:
        multiext(config['output_folder'] + '/wasp-processing/bwa_index/genome_index', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100000,
        time_min=lambda wildcards, attempt: attempt * 50
    params:
        prefix=config['output_folder'] + '/wasp-processing/bwa_index/genome_index'
    conda:
        '../envs/wasp.yaml'
    shell:
        '''
        bwa-mem2 \
        index \
        -p {params.prefix} \
        {input}
        '''
