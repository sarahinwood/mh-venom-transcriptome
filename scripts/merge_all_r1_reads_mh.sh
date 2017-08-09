#!/usr/bin/env bash

cat output/bbmerge/head_r1_unmerged.fq.gz output/bbmerge/head_merged.fq.gz > output/bbmerge/head_all_r1.fq.gz

cat output/bbmerge/abdo_r1_unmerged.fq.gz output/bbmerge/abdo_merged.fq.gz > output/bbmerge/abdo_all_r1.fq.gz

cat output/bbmerge/sting_r1_unmerged.fq.gz output/bbmerge/sting_merged.fq.gz > output/bbmerge/sting_all_r1.fq.gz