#!/usr/bin/python3
import os
import sys
import re
from argparse import ArgumentParser

from Bio.SwissProt import parse
import gffutils


def uniprot_to_geneid(name):  # some issues with uniprot files that contain entries from multiple genomes
    # ORF names (a.k.a. sequencing names or contig names or temporary ORFNames).
    # A name temporarily attributed by a sequencing project to an open reading frame.
    # This name is generally based on a cosmid numbering system.
    # Examples: MtCY277.28c, SYGP-ORF50, SpBC2F12.04, C06E1.1, CG10954.
    if "ORFNames=" in name:
        start_position = name.find("ORFNames=") + 9
        gene_id = name[start_position:].split(' ', 1)[0]
        return gene_id
    if "OrderedLocusNames=" in name:
        start_position = name.find("OrderedLocusNames=") + 18
        gene_id = name[start_position:].split(' ', 1)[0]
        return gene_id
    return "nogeneid"


def geneid_to_genome_no(gene_name):
    lastn = ''
    for n in re.findall(r'\d+', gene_name):  # find last number in string
        lastn = n
    try:
        gene_no = lastn
        genome_id = gene_name[0:len(gene_name) - len(lastn)]
        return genome_id, gene_no
    except IndexError:
        return "", ""


def parse_uniprot_cross_ref(uniprot_cr_tag, resource_abbr, indices=None):
    cross_ref_list = []
    for cr in uniprot_cr_tag:
        if cr[0].lower() == resource_abbr.lower():
            try:
                if indices is None:
                    cross_ref_list.append(cr[1:])
                else:
                    cross_ref_list.append([r for idx, r in enumerate(cr) if idx in indices])
            except IndexError:
                continue
    return cross_ref_list


def get_pfam(uniprot_cr_tag):
    pfam_list = parse_uniprot_cross_ref(uniprot_cr_tag, "Pfam")
    if len(pfam_list) == 0:
        pfam_str = "#"
    else:
        pfam_str = ""
        for pf in pfam_list:
            try:
                pfam_str += pf[0].lower()
            except IndexError:
                continue
    return pfam_str


def get_arabidopsis_gene_id(uniprot_cr_tag):
    araport_list = parse_uniprot_cross_ref(uniprot_cr_tag, "Araport", [1])
    tair_list = parse_uniprot_cross_ref(uniprot_cr_tag, "TAIR", [2])
    ensemblplants_list = parse_uniprot_cross_ref(uniprot_cr_tag, "EnsemblPlants", [3])
    gene_id_list = set([cr_id for cr_list in araport_list + tair_list + ensemblplants_list for cr_id in cr_list])
    try:
        return sorted(gene_id_list)[0]
    except IndexError:
        return None


def parse_gff_for_info(gff_path, feature_type='gene', id_list=None):
    gff_info_dict = {}
    with open(gff_path, 'r') as gff_handle:
        for line in gff_handle:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            gff_line_info = gffutils.feature.feature_from_line(line)
            if gff_line_info.featuretype == feature_type:
                gff_seqid = gff_line_info.seqid
                gff_pos = str(gff_line_info.start) + '-' + str(gff_line_info.end)
                gff_strand = gff_line_info.strand
                gff_attributes = gff_line_info.attributes
                # gff_attr_id = gff_attributes['ID']
                feature_id_prefix_regex = re.compile('^' + feature_type + '\W+')
                try:
                    feature_id = re.sub(feature_id_prefix_regex, '', gff_attributes['ID'][0])
                except IndexError:
                    feature_id = ''
                if id_list is not None and feature_id not in id_list:
                    continue
                gff_info_dict.setdefault(feature_id, (gff_seqid, gff_pos, gff_strand))
    return gff_info_dict


def main(file_path, output_fasta, species=None, gff_path=None):
    if gff_path is not None and os.path.isfile(gff_path):
        gff_info_dict = parse_gff_for_info(gff_path)
    else:
        gff_info_dict = None
    with open(file_path, 'r') as input_handle, open(output_fasta, 'w') as output_handle:
        records = parse(input_handle)
        for record in records:
            strand = "#"
            chromosome = "#"
            # Orginally the field after strand is 'cog'
            # cog = "#"

            gene_name = uniprot_to_geneid(record.gene_name)
            gene_name = re.sub(r'^\W+', '', gene_name)
            gene_name = re.sub(r'\W+$', '', gene_name)

            if species is not None and species == "arabidopsis":
                gene_name = get_arabidopsis_gene_id(record.cross_references)
                gene_pos = "#"
            else:
                gene_name, gene_pos = geneid_to_genome_no(gene_name)
            if gene_name == "":
                continue
            if gff_info_dict is not None:
                try:
                    gff_seqid, gff_pos, strand = gff_info_dict[gene_name]
                    chromosome = gff_seqid
                    gene_pos = gff_pos
                except KeyError:
                    print("No GFF info for", gene_name)

            pfam = str(get_pfam(record.cross_references))
            taxid = str(record.taxonomy_id[0]).replace(" ", "_")
            try:
                output_handle.write(">" + '|'.join([gene_name, gene_pos, strand, chromosome, pfam, taxid,
                                                    str(record.accessions[0])]) + '\n')
                output_handle.write(record.sequence + '\n')
            except TypeError:
                print("Gene Name not found", record.gene_name, gene_pos, strand, chromosome, pfam, taxid) 


if __name__ == '__main__':
    parser = ArgumentParser(usage="uniprot_to_fa_annotated.py uniprot_text_format_file.txt "
                                  "uniprot_text_sequence_output.fasta")
    parser.add_argument('-i', '--input', dest='uniprot_txt', required=True)
    parser.add_argument('-o', '--output', dest='output_fasta')
    parser.add_argument('-g', '--gff', dest='gff_path')
    parser.add_argument('-s', '--species', dest='species')
    args = parser.parse_args()

    if args.output_fasta is None:
        output_path = args.uniprot_txt + '.fa'
    else:
        output_path = args.output_fasta
    # print(args)
    main(args.uniprot_txt, output_path, args.species, args.gff_path)

