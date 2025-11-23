# Proximity coverage
Repository of the proximity coverage study

## Get ready

Place the snakemake profile config in the working directory.

```sh
mkdir my_working_dir
cd my_working_dir
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/config.yaml
```

## Test 1: dogs

### Get the code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.yaml
```

### Prepare the config files

```sh
nano 0_preprocessing.yaml
```

### Execute the pipelines

#### 0 - Preprocessing

```sh
screen -r hmmo
module load snakemake/9.9.0
snakemake -s 0_preprocessing.smk --cores 8 --profile .
```