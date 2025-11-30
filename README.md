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
wget http://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/1_single_coverage.smk
wget https://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/2_all_coverage.smk
wget https://raw.githubusercontent.com/anttonalberdi/proximity_coverage/refs/heads/main/3_proximity_coverage.smk

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
screen -r proximity3
module load snakemake/9.9.0
snakemake -s 2_all_coverage.smk --cores 8 --profile .
```

#### 3 - Proximity coverage

Runs binning using an increasing number of coverage profiles, sorted by metagenomic similarity to the focal sample.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 3_proximity_coverage.smk --cores 24 --profile .
```

#### 4 - Random coverage

Runs binning using an increasing number of randomly selected coverage profiles.

```sh
screen -r proximity1
module load snakemake/9.9.0
snakemake -s 4_random_coverage.smk --cores 24 --profile .
```

#### 5 - List and cluster genomes

List genomes under comparison and run dRep.

```sh
useconda
module load drep/3.6.2 fastani/1.33 mash/2.3

SAMPLES=($(ls 1_single_coverage/bowtie2/*.bam | sed 's/.*\///; s/\.bam$//'))

for ID in "${SAMPLES[@]}"; do
    echo "Processing $ID"

    python list_dereplicated_genomes.py \
        --assembly "$ID" \
        --binner metabat2 \
        --output "5_clustering/$ID"

    sbatch --partition=cpuqueue \
        --job-name="drep_$ID" \
        --qos=normal \
        --cpus-per-task=4 \
        --mem=8G \
        --time=0:15:00 \
        --wrap="dRep compare 5_clustering/$ID -g 5_clustering/${ID}/${ID}_metabat2_paths.txt"
done
```

#### 6 - Compile cluster data

Compile data required for statistical analyses and visualisations in R.

```sh
mkdir -p 6_final
ls 5_clustering > 6_final/clusters.csv
cat 5_clustering/*/*_genomeInfo.csv > 6_final/genomeInfo.csv
cat 5_clustering/*/data_tables/Ndb.csv > 6_final/aniDist.csv
```












#### 7 - Compute cluster statistics

Compute cluster statistics.

```sh
ID="EHI01280"
python compute_cluster_stats.py -m 10 -n "5_clustering/$ID/data_tables/Ndb.csv" -g "5_clustering/$ID/${ID}_metabat2_genomeInfo.csv" -o "5_clustering/$ID"
cat 5_clustering/$ID/binning_metrics.csv

ID="EHI01283"
python compute_cluster_stats.py -m 10 -n "5_clustering/$ID/data_tables/Ndb.csv" -g "5_clustering/$ID/${ID}_metabat2_genomeInfo.csv" -o "5_clustering/$ID"
cat 5_clustering/$ID/binning_metrics.csv

ID="EHI01321"
python compute_cluster_stats.py -m 10 -n "5_clustering/$ID/data_tables/Ndb.csv" -g "5_clustering/$ID/${ID}_metabat2_genomeInfo.csv" -o "5_clustering/$ID"
cat 5_clustering/$ID/binning_metrics.csv
```