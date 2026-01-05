# workflow/rules/wasp.run.smk


rule make_chrom_info:
    '''
    Create chromosome info file - CHROMOSOME_NAME\tCHROMOSOME_LENGTH
    '''
    input:
        config['reference_genome']
    output:
        config['output_folder'] + '/wasp-processing/index/chrom_info.txt.gz'
    threads: 1
    params:
        genome_fasta_index=config['reference_genome'] + '.fai'
    resources:
        mem_mb=500,
        time="0:05:00"
    shell:
        '''
        if [ ! -f {params.genome_fasta_index} ]; then
            samtools faidx {input}
        fi
        cut -f1,2 {params.genome_fasta_index} | gzip > {output}
        '''

rule wasp_index:
    '''
    Generate h5 files required by WASP
    '''
    input:
        info=rules.make_chrom_info.output,
        vcf_files=VCF_FILES
    output:
        haplotypes=config['output_folder'] + '/wasp-processing/index/haplotypes.h5',
        snp_index=config['output_folder'] + '/wasp-processing/index/snp_index.h5',
        snp_tab=config['output_folder'] + '/wasp-processing/index/snp_tab.h5'
    threads: 1
    resources:
        mem_mb=50000,
        time="01:00:00"
    params:
        wasp_dir=config['wasp_dir']
    shell:
        '''
        {params.wasp_dir}/snp2h5/snp2h5 \
        --chrom {input.info} \
        --haplotype {output.haplotypes} \
        --snp_index {output.snp_index} \
        --snp_tab {output.snp_tab} \
        --format vcf \
        {input.vcf_files}
        '''

