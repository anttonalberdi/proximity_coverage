# Proximity coverage
Repository of the proximity coverage study

## Get ready

Place the snakemake profile config in the working directory.

```sh
mkdir /projects/alberdilab/people/jpl786/2026_proximity_coverage/dogs
cd /projects/alberdilab/people/jpl786/2026_proximity_coverage/dogs
wget https://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/config.yaml
```

## Test 1: dogs

### Get the code

```sh
wget https://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/0_preprocessing.smk
wget https://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/0_preprocessing.yaml
```

### Prepare the config files

```sh
nano 0_preprocessing.yaml
```

### Execute the pipelines

#### 0 - Preprocessing

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 0_preprocessing.smk --cores 8 --profile .
```

#### 1 - Single coverage

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 1_single_coverage.smk --cores 8 --profile .
```

#### 2 - All coverage

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 1_single_coverage.smk --cores 8 --profile .
```

#### 4 - Random coverage

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 1_single_coverage.smk --cores 8 --profile .
```