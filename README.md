# Rasqual pipeline
Rasqual pipeline using snakemake


## Initial Setup

```bash
git clone --recursive https://github.com/davidebolo1993/rasqual_smk_pipeline
cd rasqual_smk_pipeline

#create environment (conda/mamba/micromamba)
mamba create -n rasqual_smk_pipeline_env \
	-c bioconda \
	-c conda forge \
	snakemake=7.32.4 \
	apptainer=1.3.2 \
	cookiecutter=2.6.0 \
	hdf5
	
mamba activate $PWD/rasqual_smk_pipeline_env

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

```bash
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template
```

### LSF

```bash
#Add here
```

## Run 

### SLURM

```bash

snakemake \
	--profile config/slurm \
	--rerun-triggers mtime \
	--use-conda
```

### LSF

```bash
#Add here
```

