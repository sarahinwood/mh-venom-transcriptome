#!/usr/bin/env bash

bbduk="/Volumes/BiocArchive/archive/deardenlab/asw-transcriptome/bin/bbmap/bbduk.sh"

outdir="output/bbduk_trim"

head_r1="data/head_r1.fastq.gz"
head_r2="data/head_r2.fastq.gz"

abdo_r1="data/abdo_r1.fastq.gz"
abdo_r2="data/abdo_r2.fastq.gz"

sting_r1="data/sting_r1.fastq.gz"
sting_r2="data/sting_r2.fastq.gz"

"${bbduk}" in="${head_r1}" in2="${head_r2}" \
    out=output/bbduk_trim/head_r1.fq.gz \
    out2=output/bbduk_trim/head_r2.fq.gz \
    ref=/Volumes/BiocArchive/archive/deardenlab/asw-transcriptome/bin/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=15 \
    &> "${outdir}"/head.log.txt &

"${bbduk}" in="${abdo_r1}" in2="${abdo_r2}" \
    out=output/bbduk_trim/abdo_r1.fq.gz \
    out2=output/bbduk_trim/abdo_r2.fq.gz \
    ref=/Volumes/BiocArchive/archive/deardenlab/asw-transcriptome/bin/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=15 \
    &> "${outdir}"/abdo.log.txt &

"${bbduk}" in="${sting_r1}" in2="${sting_r2}" \
    out=output/bbduk_trim/sting_r1.fq.gz \
    out2=output/bbduk_trim/sting_r2.fq.gz \
    ref=/Volumes/BiocArchive/archive/deardenlab/asw-transcriptome/bin/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=15 \
    &> "${outdir}"/sting.log.txt &

wait

echo "done"

exit 0
