#!/bin/bash
#PBS -l select=1:ncpus=4
#PBS -N download_qc_align
#PBS -q bix
#PBS -M dylang@sun.ac.za
#PBS -m abe
#PBS -e /home/dylang/dong_reseq/logs/1_1_download_qc_align.err
#PBS -o /home/dylang/dong_reseq/logs/1_1_download_qc_align.out
#PBS -l walltime=120:00:00

#Module load
module load app/parallel/parallel

#Paths for downloading and processing
base_reseq_dir="/home/dylang/dong_reseq"
accessions_file="${base_reseq_dir}/meta_files/dong_accessions.txt"
downloaded_accessions_file="${base_reseq_dir}/meta_files/downloaded_dong_accessions.txt"
snp_chip_accessions_file="${base_reseq_dir}/meta_files/dong_accessions_snp_chip_overlap_final.txt"
excluded_accessions_file="${base_reseq_dir}/meta_files/excluded_dong_accessions.txt"
meta_file="${base_reseq_dir}/meta_files/dong_meta_Vv_Vs_final.txt"
qc="${base_reseq_dir}/workflow/1_2_qc_script.job"

#Parallel download function --> updated link (>15 MB/s)
download_fast() {
    accession="$1"
    #url="https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub$accession"
     url="ftp://download.big.ac.cn/gsa2$accession"
    #url="https://download.cncb.ac.cn/gsa2$accession"
     wget "$url" -P "${reseq_dir}" --continue
}

export -f download_fast

#Number accessions to download and process at a time (based on entry list in dong_accessions)
total_files=100
split_size=10
start=3901
end=4000
bins=$((total_files / split_size))

#Split accessions file into bins of split_size sequences each
for i in $(seq 0 $((bins - 1))); do
    bin_start=$((start + i * split_size))
    bin_end=$((start + (i + 1) * split_size - 1))
    if [ $bin_end -gt $end ]; then
        bin_end=$end
    fi

    reseq_dir="${base_reseq_dir}/reseq_${bin_start}_${bin_end}"
    mkdir -p "${reseq_dir}"

    bin_file="${reseq_dir}/$(basename "$reseq_dir").txt"
    sed -n "${bin_start},${bin_end}p" "${accessions_file}" > "${bin_file}"
    
    #Filter out snp_chip downloaded accessions
    temp_bin_file="${reseq_dir}/temp_$(basename "{$bin_file}")"
    grep -vxFf "${snp_chip_accessions_file}" "${bin_file}" > "${temp_bin_file}"
    mv "${temp_bin_file}" "${bin_file}"
    
    if [ ! -s "${bin_file}" ]; then
    	echo "All sequences already downloaded. Skipping ${bin_start} to ${bin_end}"
    	continue
    fi
    
    #Increase file descriptor limit
    ulimit -n 4096
    
    #Parallel download accessions in current bin 
    cd "${reseq_dir}" 
    #Could increase -j according to requirements
    parallel -j 4 download_fast :::: "${bin_file}"
    echo "Successfully downloaded accessions for ${bin_start} to ${bin_end}"


    #Compare paired-end file sizes; append to downloaded_dong_accessions.txt
    for f1_file in *_f1.fq.gz; do
        r2_file="${f1_file/_f1/_r2}"
        if [ ! -f "${r2_file}" ]; then
            echo "r2 file not found for ${f1_file}; exiting"
            exit 1
        fi
        echo "{r2_file} found for ${f1_file}; continue"

		#Calculate file size differences for each read-pair
        f1_size=$(stat -c "%s" "${f1_file}")
        r2_size=$(stat -c "%s" "${r2_file}")
        size_diff=$((f1_size - r2_size))
        abs_size_diff=${size_diff#-}

		#Print sizes and size differences for each read-pair
        echo "f1 file: ${f1_file}, size: ${f1_size}"
        echo "r2 file: ${r2_file}, size: ${r2_size}"
        echo "Size difference: ${abs_size_diff} bytes"
        
        #Check size conditions; reattempt download
        if [ $abs_size_diff -le 1073741824 ] \
        	&& [ $f1_size -ge 536870912 ] \
        	&& [ $r2_size -ge 536870912 ] \
        	&& [ $f1_size -le  10737418240 ] \
        	&& [ $r2_size -le  10737418240 ]; then
        	echo "Size conditions are met, adding accessions to downloaded_accessions.txt"
            f1_filename=$(basename "${f1_file}")
            r2_filename=$(basename "${r2_file}")
        	echo "${f1_filename} ${r2_filename}" >> "${downloaded_accessions_file}"
        else
            echo "Size conditions are not met, removing and re-downloading files"
            rm "${f1_file}" "${r2_file}"
            accession_number=$(basename "${f1_file}" | cut -d "_" -f 1)
            #wget "https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA006917/${accession_number}/$(basename "${f1_file}")" -P "${reseq_dir}"
            #wget "https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA006917/${accession_number}/$(basename "${r2_file}")" -P "${reseq_dir}"
            wget "ftp://download.big.ac.cn/gsa2/CRA006917/${accession_number}/$(basename "${f1_file}")" -P "${reseq_dir}"
            wget "ftp://download.big.ac.cn/gsa2/CRA006917/${accession_number}/$(basename "${r2_file}")" -P "${reseq_dir}"
            f1_size=$(stat -c "%s" "${f1_file}")
            r2_size=$(stat -c "%s" "${r2_file}")
            size_diff=$((f1_size - r2_size))
            abs_size_diff=${size_diff#-}

			#Recheck file size conditions for each read-pair after redownload
            if [ $abs_size_diff -le 1073741824 ] \
        		&& [ $f1_size -ge 536870912 ] \
        		&& [ $r2_size -ge 536870912 ] \
        		&& [ $f1_size -le  10737418240 ] \
        		&& [ $r2_size -le  10737418240 ]; then
            	echo "Size conditions are met after re-downloading, adding to downloaded_accessions"                         
            	f1_filename=$(basename "${f1_file}")
            	r2_filename=$(basename "${r2_file}")
            	echo "${f1_filename} ${r2_filename}" >> "${downloaded_accessions_file}"
            else
                echo "Sizes are still different after re-downloading; adding to excluded_accessions"
                echo "${accession_number}" >> "${excluded_accessions_file}"
                echo "Removing ${f1_file} and ${r2_file}"
                rm "${f1_file}" "${r2_file}"
            fi
        fi
        echo "-------------------------------------"
    done
    
    #Rename files --> find corresponding accession number in metafile and rename accordingly
    while IFS= read -r sequence_file; do
        accession_number=$(echo "${sequence_file}" | cut -d "_" -f 1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
        echo "Accession number: ${accession_number}"
        database_id=$(grep -F "${accession_number}" "${meta_file}" | awk -F'\t' '{print $5}')
        echo "Database_ID: ${database_id}"
        if [ -n "${database_id}" ]; then
            original_file_f1=$(find "${reseq_dir}" -name "${accession_number}_f1.fq.gz")
            original_file_r2=$(find "${reseq_dir}" -name "${accession_number}_r2.fq.gz")
            if [ -n "${original_file_f1}" ] && [ -n "${original_file_r2}" ]; then
                new_name_f1="${reseq_dir}/${database_id}_f1.fq.gz"
                new_name_r2="${reseq_dir}/${database_id}_r2.fq.gz"
                mv "${original_file_f1}" "${new_name_f1}"
                mv "${original_file_r2}" "${new_name_r2}"
                echo "Renamed ${original_file_f1} to ${new_name_f1}"
                echo "Renamed ${original_file_r2} to ${new_name_r2}"
            else
                echo "Original sequence files not found for ${accession_number}"
            fi
        else
            echo "Database_ID not found for ${accession_number}"
        fi
    done < "${downloaded_accessions_file}"

    #Submit QC script for current bin
    ssh hpc2 qsub -v bin_file="${bin_file}",reseq_dir="${reseq_dir}" "$qc"
done
