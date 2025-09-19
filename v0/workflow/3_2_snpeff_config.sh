#!/bin/bash
#PBS -l select=1:ncpus=2:mem=5GB
#PBS -q bix
#PBS -l walltime=01:00:00
#PBS -o /home/dylang/dong_reseq/logs/3_2_snpeff_config.out
#PBS -e /home/dylang/dong_reseq/logs/3_2_snpeff_config.err
#PBS -m abe
#PBS -M dylang@sun.ac.za

#Load modules
module load app/snpeff/current
module load app/gffread
module load app/samtools

#Paths & variables
BASE_REF="/home/dylang/PN_ref"
FASTA_REF="${BASE_REF}/T2T_ref.fasta"
ANN_REF="${BASE_REF}/T2T_v5.gff"
SNPEFF_DIR="/apps2/snpeff/snpEff/snpEff"
ORGANISM="Vitis vinifera"
VERSION="T2T"
ACCESSION="GCF_030704535.1"
DB_DIR="${SNPEFF_DIR}/data/${ACCESSION}"
INPUT_DIR="/home/dylang/dong_reseq/final_filtering" #--> This points to the filtered output dir of bcftools/vcftools


#Check whether existing database exists
if [[ -d "${DB_DIR}" ]] && [[ -f "${DB_DIR}/snpEffectPredictor.bin" ]]; then
    echo "Database already exists; skipping build..."
else
    echo "Creating snpeff database..."

    #Validate required input files
    if [[ -f "${FASTA_REF}" ]] && [[ -f "${ANN_REF}" ]]; then
    echo "Found ${FASTA_REF} and ${ANN_REF}"

        #Check indexing
      if [ ! -f "${FASTA_REF}.fai" ]; then
        echo "Indexing reference"
        samtools faidx "${FASTA_REF}" || { echo "Failed to generate index for ${FASTA_REF}"; exit 1; }
      else
        echo "Found existing index"
      fi

        #Create CDS fasta
        if ! gffread -W -x "${BASE_REF}/cds.fa" -g "${FASTA_REF}" "${ANN_REF}"; then
            echo "Failed to extract CDS sequences"
            exit 1
        fi

        #Create protein fasta
        if ! gffread -y "${BASE_REF}/protein.fa" -g "${FASTA_REF}" "${ANN_REF}"; then
            echo "Failed to extract protein sequences"
            exit 1
        fi

        #Copy to snpeff data dir
        mkdir "${DB_DIR}"

        cp "${FASTA_REF}" "${DB_DIR}/sequences.fa"
        cp "${FASTA_REF}.fai" "${DB_DIR}/sequences.fa.fai"
        cp "${ANN_REF}" "${DB_DIR}/genes.gff"
        cp "${BASE_REF}/cds.fa" "${DB_DIR}/cds.fa"
        cp "${BASE_REF}/protein.fa" "${DB_DIR}/protein.fa"

        echo "Copied input files to ${SNPEFF_DIR}"

        #Update config file
        CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"
        if ! grep -q "${ACCESSION}" "${CONFIG_FILE}"; then
            echo "Adding entry to ${CONFIG_FILE}"

            sed -i "/^# Non-standard Databases$/a # ${ORGANISM} : version ${VERSION}\n${ACCESSION}.genome : ${ORGANISM}" "${CONFIG_FILE}"
        else
            echo "Config entry already exists"
        fi

        #Run snpeff database build
        cd "${SNPEFF_DIR}"
        java -jar snpEff.jar build -gff3 -v "${ACCESSION}" || \
            { echo "snpEff database build failed"; exit 1; }

        echo "Successful database creation for ${ACCESSION}"
    else
        echo "Failed to find ${FASTA_REF} and ${ANN_REF}"
        exit 1
    fi
fi

#Run snpeff annotation
echo "Running snpEff annotation"

java -jar -Xmx4g snpEff.jar -v "${ACCESSION}" -stats snpeff_stats.html \
    "${INPUT_DIR}/final_filt_maxmiss.vcf.gz" | gzip > "${INPUT_DIR}/final_filt_maxmiss_ann.vcf.gz" || \
    { echo "snpEff annotation failed"; exit 1; }

echo "Successfully annotated ${INPUT_DIR}/final_filt_maxmiss.vcf.gz"

## This needs to be implemented as well
## SNPSIFT

#Filter PCG variants --> snpsift subworkflow
zcat "${INPUT_DIR}/final_filt_maxmiss_ann.vcf.gz" \
    | java -jar SnpSift.jar filter \
    "(ANN[*].BIOTYPE = 'protein_coding')" \
    | gzip > "${INPUT_DIR}/final_filt_high_mod_pcg.ann.vcf.gz" || \
    { echo "SnpSift filter failed"; exit 1;}

#echo "Successfully extracted ${input_dir}/final_filt_high_mod_pcg.ann.vcf.gz"
