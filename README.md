# Rasqual pipeline
Rasqual pipeline using snakemake


## Initial Setup

```bash
git clone --recursive https://github.com/davidebolo1993/rasqual_smk_pipeline
cd rasqual_smk_pipeline

#create environment (conda/mamba)

conda create -n rasqual_smk_pipeline_env \
	-c bioconda \
	-c conda forge \
	snakemake=7.32.4 \
	apptainer=1.3.2 \
	cookiecutter=2.6.0 \
	hdf5
	
conda activate $PWD/rasqual_smk_pipeline_env

#build scar
mkdir -p scar/build && cd scar/build
cmake ..
make
cd -

#build WASP
cd WASP/snp2h5
make
cd -
```

## Profiles 

### SLURM

An example slurm profile is included in `config/ht_slurm`

### LSF

An example lsf profile is included in `config/sanger_lsf`

## Run

### Adjust `config/config.yaml`

Adjust `config/config.yaml` to match your needs.

### SLURM

```bash

snakemake \
	--profile config/ht_slurm \
	--rerun-triggers mtime
```

### LSF

```bash
mkdir -p logs/lsf
snakemake \
    --cluster-config config/sanger_lsf/lsf.yaml \
    --cluster "bsub -q {cluster.queue} -G {cluster.group} -M {resources.mem_mb} -W {resources.time_min} -R \"select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]\" -eo logs/lsf/{rule}.{wildcards}.err -oo logs/lsf/{rule}.{wildcards}.out -n {threads}" \
    --default-resources mem_mb=1000 time_min=10 \
    -j 1 \
    --restart-times 5 \
    --rerun-triggers 'mtime' \
    --use-conda \
    --printshellcmds \
    --rerun-incomplete
```

