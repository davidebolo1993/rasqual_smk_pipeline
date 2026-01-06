# workflow/rules/wasp.run.smk

if config.get('run_wasp', True):
    # WASP is enabled - use full pipeline
    ruleorder: wasp_complete_pipeline > copy_split_bam_no_wasp
else:
    # WASP is disabled - copy to output
    ruleorder: copy_split_bam_no_wasp > wasp_complete_pipeline

rule wasp_complete_pipeline:
    '''
    Complete WASP pipeline for a single celltype/donor BAM.
    '''
    input:
        bam=config['output_folder'] + '/split-bams/{celltype}/{donor}.bam',
        qc_done=config['output_folder'] + '/qc/flagstat-split-bams/{celltype}/{donor}.qc.tsv',
        haplotype=rules.wasp_index.output.haplotypes,
        snp_index=rules.wasp_index.output.snp_index,
        snp_tab=rules.wasp_index.output.snp_tab, 
        bwa_index=rules.bwa_mem2_index.output
    output:
        final_bam=config['output_folder'] + '/processed-bams/{celltype}/{donor}.bam',
        final_bai=config['output_folder'] + '/processed-bams/{celltype}/{donor}.bam.bai',
        qc_report=config['output_folder'] + '/qc/wasp-complete/{celltype}/{donor}.qc.tsv',
        rmdup_stats=config['output_folder'] + '/qc/wasp-complete/{celltype}/{donor}.rmdup_stats.txt'
    params:
        workdir=config['output_folder'] + '/wasp-processing/{celltype}/{donor}',
        wasp_dir=config['wasp_dir'],
        genome_index=config['output_folder'] + '/wasp-processing/bwa_index/genome_index',
        donor='{donor}',
        celltype='{celltype}'
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100000,
        time_min=lambda wildcards, attempt: attempt * 150
    conda:
        '../envs/wasp.yaml'
    shell:
        '''        
        mkdir -p {params.workdir}
        
        # Initialize QC file
        echo -e "sample\\tcell\\tstep\\treads\\tpairs\\tmapped\\tproperly_paired\\tsingletons\\tread1\\tread2" > {output.qc_report}
        
        # Function to collect QC stats
        collect_qc() {{
            local bam_file=$1
            local step_name=$2
            
            flagstat=$(samtools flagstat -@ 2 "$bam_file")
            
            reads=$(echo "$flagstat" | grep "in total" | awk '{{print $1}}')
            pairs=$(echo "$flagstat" | grep "paired in sequencing" | awk '{{print $1}}')
            mapped=$(echo "$flagstat" | grep "mapped (" | head -1 | awk '{{print $1}}')
            proper=$(echo "$flagstat" | grep "properly paired" | awk '{{print $1}}')
            singletons=$(echo "$flagstat" | grep "singletons" | awk '{{print $1}}')
            read1=$(echo "$flagstat" | grep "read1" | awk '{{print $1}}')
            read2=$(echo "$flagstat" | grep "read2" | awk '{{print $1}}')
            
            echo -e "{params.donor}\\t{params.celltype}\\t$step_name\\t$reads\\t$pairs\\t$mapped\\t$proper\\t$singletons\\t$read1\\t$read2" >> {output.qc_report}
        }}
        
        # QC: Original input
        collect_qc "{input.bam}" "0_original"
        
        # ========================================================================
        # STEP 1: Name-sort, fixmate, filter pairs, remap
        # ========================================================================

        STEP1DIR="{params.workdir}/01.remap"
        mkdir -p $STEP1DIR
        
        samtools collate \
        -@ {threads} \
        -u \
        -O \
        {input.bam} $STEP1DIR/{params.donor}.collate | \
        samtools fixmate \
        -m \
        -u \
        -@ {threads} \
        - - | \
        samtools view \
        -u \
        -f 1 \
        -F 12 \
        - | \
        samtools sort \
        -n \
        -u \
        -@ {threads} \
        -T $STEP1DIR/{params.donor}.sort - | \
        samtools fastq \
        -@ {threads} \
        -T CB \
        -1 $STEP1DIR/{params.donor}.r1.fq.gz \
        -2 $STEP1DIR/{params.donor}.r2.fq.gz \
        -0 /dev/null \
        -s /dev/null \
        -n \
        -
        bwa-mem2 mem \
        -C \
        -t {threads} \
        {params.genome_index} \
        $STEP1DIR/{params.donor}.r1.fq.gz \
        $STEP1DIR/{params.donor}.r2.fq.gz | \
        samtools view \
        -@ {threads} \
        -b \
        -q 10 \
        - | \
        samtools sort \
        -@ {threads} \
        -T $STEP1DIR/{params.donor} \
        --write-index \
        -o $STEP1DIR/{params.donor}.bam
       
        collect_qc "$STEP1DIR/{params.donor}.bam" "1_remap_input"
        
        # ========================================================================
        # STEP 2: Identify reads overlapping SNPs
        # ========================================================================

        STEP2DIR="{params.workdir}/02.intersect"
        
        python {params.wasp_dir}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir $STEP2DIR \
        --snp_tab {input.snp_tab} \
        --snp_index {input.snp_index} \
        --haplotype {input.haplotype} \
        --samples {params.donor} \
        $STEP1DIR/{params.donor}.bam
        
        collect_qc "$STEP2DIR/{params.donor}.keep.bam" "2a_keep_no_snps"
        collect_qc "$STEP2DIR/{params.donor}.to.remap.bam" "2b_to_remap"
        
        # Add CB tags to FASTQ for remapping
        python workflow/scripts/add_cb_to_fastq.py \
        $STEP2DIR/{params.donor}.to.remap.bam \
        $STEP2DIR/{params.donor}.remap.fq1.gz \
        $STEP2DIR/{params.donor}.remap.fq2.gz \
        $STEP2DIR/{params.donor}.remap.cb
        
        # ========================================================================
        # STEP 3: Remap reads with swapped alleles
        # ========================================================================

        STEP3DIR="{params.workdir}/03.remap"
        mkdir -p $STEP3DIR
        
        bwa-mem2 mem \
        -C \
        -t {threads} \
        {params.genome_index} \
        $STEP2DIR/{params.donor}.remap.cb.fq1.gz \
        $STEP2DIR/{params.donor}.remap.cb.fq2.gz | \
        samtools view \
        -@ {threads} \
        -b \
        -q 10 - | \
        samtools sort \
        -T $STEP3DIR/{params.donor} \
        -@ {threads} \
        --write-index \
        -o $STEP3DIR/{params.donor}.bam

        collect_qc "$STEP3DIR/{params.donor}.bam" "3_remapped"
        
        # ========================================================================
        # STEP 4: Filter reads showing mapping bias
        # ========================================================================

        STEP4DIR="{params.workdir}/04.filter"
        mkdir -p $STEP4DIR
        
        python {params.wasp_dir}/mapping/filter_remapped_reads.py \
        $STEP2DIR/{params.donor}.to.remap.bam \
        $STEP3DIR/{params.donor}.bam \
        $STEP4DIR/{params.donor}.bam {log}
        
        collect_qc "$STEP4DIR/{params.donor}.bam" "4_bias_filtered"
        
        # ========================================================================
        # STEP 5: Merge keep + bias-filtered reads
        # ========================================================================

        STEP5DIR="{params.workdir}/05.merge"
        mkdir -p $STEP5DIR
        
        samtools merge -@ {threads} \
        -o $STEP5DIR/{params.donor}.notsort.bam \
        $STEP4DIR/{params.donor}.bam \
        $STEP2DIR/{params.donor}.keep.bam
        
        samtools sort \
        -T $STEP5DIR/{params.donor} \
        -@ {threads} \
        --write-index \
        -o $STEP5DIR/{params.donor}.bam \
        $STEP5DIR/{params.donor}.notsort.bam
        
        collect_qc "$STEP5DIR/{params.donor}.bam" "5_merged"
        
        # ========================================================================
        # STEP 6: Remove duplicates
        # ========================================================================

        STEP6DIR="{params.workdir}/06.rmdup"
        mkdir -p $STEP6DIR
        
        python workflow/scripts/rmdup.py \
        $STEP5DIR/{params.donor}.bam \
        $STEP6DIR/{params.donor}.bam 2> {output.rmdup_stats}
        
        samtools sort -@ {threads} \
        -o {output.final_bam}##idx##{output.final_bai} \
        -T $STEP6DIR/{params.donor} \
        --write-index \
        $STEP6DIR/{params.donor}.bam

        collect_qc "{output.final_bam}" "6_final_dedup"
        '''


rule copy_split_bam_no_wasp:
    '''
    Skip WASP processing - just copy/index split BAMs to processed folder.
    '''
    input:
        bam=config['output_folder'] + '/split-bams/{celltype}/{donor}.bam',
        qc_done=config['output_folder'] + '/qc/flagstat-split-bams/{celltype}/{donor}.qc.tsv'
    output:
        bam=config['output_folder'] + '/processed-bams/{celltype}/{donor}.bam',
        bai=config['output_folder'] + '/processed-bams/{celltype}/{donor}.bam.bai'
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
        time_min=lambda wildcards, attempt: attempt * 20
    params:
        bai=config['output_folder'] + '/split-bams/{celltype}/{donor}.bai'
    shell:
        '''
        cp {input.bam} {output.bam}
        cp {params.bai} {output.bai}
        '''


def get_processed_bam_targets(wildcards):
    '''
    Get final processed BAM targets based on WASP config.
    '''
    # Ensure checkpoint has run
    checkpoints.aggregate_split_bams.get(**wildcards).output.done
    
    # Read filtered metadata
    filtered_df = pd.read_table(config['output_folder'] + '/metadata-filtered/metadata.filtered.tsv')
    
    # Get unique celltype/donor combinations
    combinations = filtered_df[[config["celltype_column"], "donor_id"]].drop_duplicates()
    
    # Build paths for processed BAMs
    bam_paths = []
    for _, row in combinations.iterrows():
        path = f"{config['output_folder']}/processed-bams/{row[config['celltype_column']]}/{row['donor_id']}.bam"
        bam_paths.append(path)
    
    return bam_paths
