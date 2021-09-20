#!/bin/bash
# Limit background jobs to no more that $maxproc at once.
maxproc=1

# Declare an array of string with type
declare -a l=("UB-2845-11"
"UB-2845-12"
"UB-2845-13"
"UB-2845-14"
"UB-2845-21"
"UB-2845-22"
"UB-2845-23"
"UB-2845-24"
"UB-2845-31"
"UB-2845-32"
"UB-2845-33"
"UB-2845-34")

declare -a s=("S11"
"S12"
"S13"
"S14"
"S21"
"S22"
"S23"
"S24"
"S31"
"S32"
"S33"
"S34")

#S=11

for i in ${!l[@]}; do
    while [ $(jobs | wc -l) -ge "$maxproc" ]
    do
        sleep 1
    done
    #S=$((S + 1))
    #echo $S
    #echo starting new job $i with ongoing=$(jobs | wc -l)
    printf "Starting new job %s for %s with ongoing=$(jobs | wc -l)\n" "${s[i]}" "${l[i]}" 
    cellranger count --id=${s[i]} \
                   --transcriptome=/mnt/WD2/BIN21P006_ML/refdata-gex-GRCh38-2020-A \
                   --fastqs=/mnt/WD2/BIN21P006_ML/UB-2845/210427_A00181_0279_BHC5TNDRXY \
                   --sample=${l[i]} \
                   --expect-cells=5000 \
                   --localcores=48 \
                   --localmem=128
done

