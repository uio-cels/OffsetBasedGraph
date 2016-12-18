

### Reproduces gene experiment in Coordinates and Intervals on Graph-based Reference Genomes

# Get gene dataset
if [ ! -f knownGene.txt ]; then
    echo "Fetching gene dataset ..."
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
    gunzip knownGene.txt.gz
fi

python3 gene_experiment.py grch38.chrom.sizes grch38_alternative_loci.txt knownGene.txt