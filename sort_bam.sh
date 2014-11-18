#!/bin/bash

for file in alignments/*
do
    samtools sort $file ${file}_sorted.bam
done

for file in alignments/*_sorted.bam
do
    mv $file alignments/sorted/.
done
