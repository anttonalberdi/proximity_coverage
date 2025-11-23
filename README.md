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

Runs the initial quality-filtering, host removal and individual assemblies.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 0_preprocessing.smk --cores 8 --profile .
```

#### 1 - Single coverage

Runs binning relying on a single coverage sample's coverage profile.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 1_single_coverage.smk --cores 8 --profile .
```

#### 2 - All coverage

Runs binning relying on all samples' coverage profiles.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 2_all_coverage.smk --cores 8 --profile .
```

#### 3 - Proximity coverage

Runs binning using an increasing number of coverage profiles, sorted by metagenomic similarity to the focal sample.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 3_proximity_coverage.smk --cores 8 --profile .
```

#### 4 - Random coverage

Runs binning using an increasing number of randomly selected coverage profiles.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 4_random_coverage.smk --cores 8 --profile .
```