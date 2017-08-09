#!/usr/bin/env bash

set -eu

bbmerge="bin/bbmap/bbmerge.sh"

head_r1="output/bbduk_trim/head_r1.fq.gz"
head_r2="output/bbduk_trim/head_r2.fq.gz"

abdo_r1="output/bbduk_trim/abdo_r1.fq.gz"
abdo_r2="output/bbduk_trim/abdo_r2.fq.gz"

sting_r1="output/bbduk_trim/sting_r1.fq.gz"
sting_r2="output/bbduk_trim/sting_r2.fq.gz"

"${bbmerge}" in="${head_r1}" in2="${head_r2}" \
	out=output/bbmerge/head_merged.fq.gz \
	outu1=output/bbmerge/head_r1_unmerged.fq.gz \
	outu2=output/bbmerge/head_r2_unmerged.fq.gz \
	ihist=output/bbmerge/head_ihist.txt \
	verystrict=t \
	adapters=bin/bbmap/resources/adapters.fa \
	&> output/bbmerge/head.log.txt &

"${bbmerge}" in="${abdo_r1}" in2="${abdo_r2}" \
	out=output/bbmerge/abdo_merged.fq.gz \
	outu1=output/bbmerge/abdo_r1_unmerged.fq.gz \
	outu2=output/bbmerge/abdo_r2_unmerged.fq.gz \
	ihist=output/bbmerge/abdo_ihist.txt \
	verystrict=t \
	adapters=bin/bbmap/resources/adapters.fa \
	&> output/bbmerge/abdo.log.txt &

"${bbmerge}" in="${sting_r1}" in2="${sting_r2}" \
	out=output/bbmerge/sting_merged.fq.gz \
	outu1=output/bbmerge/sting_r1_unmerged.fq.gz \
	outu2=output/bbmerge/sting_r2_unmerged.fq.gz \
	ihist=output/bbmerge/sting_ihist.txt \
	verystrict=t \
	adapters=bin/bbmap/resources/adapters.fa \
	&> output/bbmerge/sting.log.txt &

wait

echo "done"

exit 0