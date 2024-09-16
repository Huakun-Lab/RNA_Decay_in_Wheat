#!/usr/bin/env python
# -*- coding: utf-8 -*-


def count_introns(gff_file):
    transcript_introns = {}

    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                data = line.strip().split('\t')
                feature_type = data[2]
                attributes = data[8]

                if feature_type == 'mRNA':
                    transcript_id = get_transcript_id(attributes)
                    transcript_introns[transcript_id] = 0

                if feature_type == 'CDS':
                    transcript_id = get_transcript_id(attributes)
                    transcript_introns[transcript_id] += 1

    return transcript_introns


def get_transcript_id(attributes):
    id_field = attributes.split(';')[0]
    transcript_id = id_field.split('=')[1]
    return transcript_id


gff_file_path = '/Reference/Genome/gff/TRITD_HC_LC_addUTR.gff3'
transcript_introns = count_introns(gff_file_path)

for transcript_id, cds_count in transcript_introns.items():
    intron_count = max(0, cds_count - 1)
    print(f"{transcript_id} : {intron_count}")