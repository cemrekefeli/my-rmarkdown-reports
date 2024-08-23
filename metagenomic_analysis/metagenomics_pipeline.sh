#!/bin/bash/env Rscript

echo "Metagenomics Pipeline"
echo "Created by Cemre Kefeli"

# Creating directory

mkdir sra_data

# Changing directory

cd sra_data

# Array of SRA accession numbers

accessions=(
SRR5181558
SRR5181559
SRR5181560
SRR5181561
SRR5181562
SRR5181563
SRR5181564
SRR5181565
SRR5181566
SRR5181567
SRR5181568
SRR5181569
SRR5181570
SRR5181571
SRR5181572
SRR5181573
SRR5181574
SRR5181575
SRR5181576
SRR5181577
SRR5181578
SRR5181579
SRR5181580
SRR5181581
)

# Save accessions to a text file

printf "%s\n" "${accessions[@]}" > accessions.txt

echo "Downloading fastq files - Bytes 784.60 Mb - 24 samples "

# Download SRA files

prefetch --option-file accessions.txt


# Convert SRA files to FASTQ format

while read -r line; do
    fastq-dump $line
done < accessions.txt

# Downloading Qiime2 (only if not already downloaded)
if [ ! -f "qiime2-amplicon-2024.2-py38-linux-conda.yml" ]; then
    wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
fi

# Create and activate the conda environment (only if it doesn't exist)
conda_env="qiime2-amplicon-2024.2"
if ! conda env list | grep -q "\<$conda_env\>"; then
    conda env create -n $conda_env --file qiime2-amplicon-2024.2-py38-linux-conda.yml
fi
conda activate $conda_env


# Running bashfile with source 
# This is needed because I create and activate conda environment
echo "Run the bash script with 'source metagenomics_pipeline.sh'"

# Define the file name
file_name="sample_files.tsv"

# Write header to the file
echo -e "sample-id\tabsolute-filepath" > "$file_name"

# Loop through each ID in the file
while IFS= read -r sample_id; do
    file_path="$PWD/${sample_id}.fastq"
    echo -e "${sample_id}\t${file_path}" >> "$file_name"
done < "accessions.txt"

echo "'$file_name' created."

SAMPLE=sample_files.tsv
OUT_FNAME=allreads
FMT=SingleEndFastqManifestPhred33V2

echo "Creating sample metadata"

cat > sample-metadata.tsv <<EOF
"sample-id"	"AvgSpotLen"	"Bases"	"BioSample"	"Bytes"	"Experiment"	"Library	Name"	"Sample	Name"	"SRA	Study"
"SRR5181558"    496     98644215        "SAMN06167590"  58291658        "SRX2497680"    "Sample8_repl3" "Site8" "SRP096885"
"SRR5181559"    496     114827099       "SAMN06167590"  67763180        "SRX2497681"    "Sample8_repl2" "Site8" "SRP096885"
"SRR5181560"    497     66352808        "SAMN06167590"  39078443        "SRX2497682"    "Sample8_repl1" "Site8" "SRP096885"
"SRR5181561"    488     35386291        "SAMN06167589"  21776576        "SRX2497683"    "Sample7_repl3" "Site7" "SRP096885"
"SRR5181562"    487     42579908        "SAMN06167589"  26307651        "SRX2497684"    "Sample7_repl2" "Site7" "SRP096885"
"SRR5181563"    485     69248953        "SAMN06167589"  38837522        "SRX2497685"    "Sample7_repl1" "Site7" "SRP096885"
"SRR5181564"    496     49928519        "SAMN06167588"  30779767        "SRX2497686"    "Sample6_repl3" "Site6" "SRP096885"
"SRR5181565"    498     40460822        "SAMN06167588"  25011727        "SRX2497687"    "Sample6_repl2" "Site6" "SRP096885"
"SRR5181566"    488     32148174        "SAMN06167588"  18336056        "SRX2497688"    "Sample6_repl1" "Site6" "SRP096885"
"SRR5181567"    500     64543452        "SAMN06167587"	36913089        "SRX2497689"    "Sample5_repl3" "Site5" "SRP096885"
"SRR5181568"    497     9366938	"SAMN06167587" 5887417 "SRX2497690"    "Sample5_repl2" "Site5" "SRP096885"
"SRR5181569"    487     49754063        "SAMN06167587"  30827297        "SRX2497691"    "Sample5_repl1" "Site5" "SRP096885"
"SRR5181570"	496	99318646	"SAMN06167586"	58731846	"SRX2497692"	"Sample4_repl3"	"Site4"	"SRP096885"
"SRR5181571"	496	76160167	"SAMN06167586"	45201511	"SRX2497693"	"Sample4_repl2"	"Site4"	"SRP096885"
"SRR5181572"	498	113261379	"SAMN06167586"	67250802	"SRX2497694"	"Sample4_repl1"	"Site4"	"SRP096885"
"SRR5181573"	499	23123019	"SAMN06167585"	14333462	"SRX2497695"	"Sample3_repl3"	"Site3"	"SRP096885"
"SRR5181574"	492	24342780	"SAMN06167585"	15213020	"SRX2497696"	"Sample3_repl2"	"Site3"	"SRP096885"
"SRR5181575"	490	62849051	"SAMN06167585"	35962677	"SRX2497697"	"Sample3_repl1"	"Site3"	"SRP096885"
"SRR5181576"	496	14726026	"SAMN06167584"	9214184	"SRX2497698"	"Sample2_repl3"	"Site2"	"SRP096885"
"SRR5181577"	496	32332587	"SAMN06167584"	20121336	"SRX2497699"	"Sample2_repl2"	"Site2"	"SRP096885"
"SRR5181578"	490	63542472	"SAMN06167584"	36301607	"SRX2497700"	"Sample2_repl1"	"Site2"	"SRP096885"
"SRR5181579"	497	78720444	"SAMN06167583"	46878848	"SRX2497701"	"Sample1_repl3"	"Site1"	"SRP096885"
"SRR5181580"	495	46699183	"SAMN06167583"	27944549	"SRX2497702"	"Sample1_repl2"	"Site1"	"SRP096885"
"SRR5181581"	496	76811107	"SAMN06167583"	45743669	"SRX2497703"	"Sample1_repl1"	"Site1"	"SRP096885"
EOF

echo "Importing sequencing data into QIIME 2"

# Importing sequencing data into QIIME 2

qiime tools import \
   --type 'SampleData[SequencesWithQuality]' \
   --input-path $SAMPLE \
    --output-path $OUT_FNAME.qza \
   --input-format=$FMT

qiime demux summarize \
  --i-data $OUT_FNAME.qza \
  --o-visualization $OUT_FNAME.qzv

echo "Denoising the data"

TRUNC_LEN=250
INPUT=allreads.qza

qiime dada2 denoise-single --i-demultiplexed-seqs allreads.qza \
  --p-trim-left 0 --p-trunc-len $TRUNC_LEN --p-n-threads 8 \
  --p-hashed-feature-ids \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza

# Create feature table

echo "Create feature table"

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


# Extracting otu table

qiime tools export --input-path table.qza --output-path table

# Annotating OTU's

echo "Downloading trained classifier"

wget https://data.qiime2.org/2018.11/common/gg-13-8-99-nb-classifier.qza

echo "Performing taxonomic classification"

qiime feature-classifier classify-sklearn \
    --i-reads rep-seqs.qza \
    --i-classifier /home/desert_tlt/data/gg-13-8-99-nb-classifier.qza \
    --o-classification taxonomy.qza

echo "Exporting taxonomy information"

qiime tools export \  
    --input-path taxonomy.qza \
    --output-path taxonomy




