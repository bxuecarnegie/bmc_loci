#!/usr/bin/python3

import sys
import os
from argparse import ArgumentParser
from subprocess import call

from Bio.SeqIO import parse as seq_parse
from Bio.SearchIO import parse as search_parse


script_path = os.path.dirname(os.path.realpath(__file__))

ref_hmm = os.path.join("reference_data/combined_hmms.hmm")
ep_ref_hmm = os.path.join("reference_data/ep_hmms.hmm")
org_table = os.path.join("reference_data/genome_vs_org_table.txt")

org_dict = {}
with open(org_table) as file:
    for org in file:
        line = org.strip()
        name = line.split()[0]
        org = " ".join(line.split()[1:])
        org_dict[name] = org


def run_hmm(input_fasta, hmm_library, output_path, hmm_cmd, additional_cmd="", log_path="/dev/null"):
    try:
        hmm_call = " ".join([hmm_cmd, additional_cmd, "--tblout", output_path, hmm_library, input_fasta, ">", log_path])
        print(hmm_call)
        retcode = call(hmm_call, shell=True)
        if retcode < 0:
            print("Child was terminated by signal", -retcode, file=sys.stderr)
        else:
            print("Child returned", retcode, file=sys.stderr)
        return retcode
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
        return -1 * sys.maxsize


def read_hmmscan_tab_output(hmm_output, hmm_cmd, parse_type='hmmer3-tab', evalue_cutoff=None, bitscore_cuttoff=None):
    hmm_top_hit_dict = {}
    for res in search_parse(hmm_output, parse_type):
        query_hmm = res.id
        if len(res.hits) > 0:
            for i in range(0, len(res.hits)):
                aln_start = 0
                aln_end = 0
                hit_hmm = res.hits[i].id
                if hmm_cmd == 'hmmscan':
                    query_id = query_hmm
                    hit_id = hit_hmm
                    hit_score = res.hits[i].hsps[0].evalue  # evalue
                    if evalue_cutoff is not None and hit_score >= evalue_cutoff:
                        continue
                    if parse_type == 'hmmer3-text':
                        aln_start = res.hits[i].hsps[0].query_start + 1
                        aln_end = res.hits[i].hsps[0].query_end
                elif hmm_cmd == 'hmmsearch':
                    query_id = hit_hmm
                    hit_id = query_hmm
                    hit_score = res.hits[i].bitscore
                    if bitscore_cuttoff is not None and hit_score <= bitscore_cuttoff:
                        continue
                    if parse_type == 'hmmer3-text':
                        aln_start = res.hits[i].hsps[0].hit_start + 1
                        aln_end = res.hits[i].hsps[0].hit_end
                else:
                    continue
                cur_hit_id, cur_hit_score, _, _ = \
                    hmm_top_hit_dict.setdefault(query_id, (hit_id, hit_score, aln_start, aln_end))
                if hmm_cmd == 'hmmscan' and hit_score < cur_hit_score:
                    hmm_top_hit_dict[query_id] = (hit_id, hit_score, aln_start, aln_end)
                elif hmm_cmd == 'hmmsearch' and hit_score > cur_hit_score:
                    hmm_top_hit_dict[query_id] = (hit_id, hit_score, aln_start, aln_end)
    return hmm_top_hit_dict


def fa_header_to_fields(fasta_header):
    # read information from fasta file
    fasta_info = fasta_header.split("|")
    try:
        gene_name = fasta_info[0]
        gene_pos = fasta_info[1]
        directionality = fasta_info[2]
        chromosome = fasta_info[3]
        pfam_list = fasta_info[4]
        taxid = fasta_info[5]
        uniprot = fasta_info[6]
        return gene_name, gene_pos, directionality, chromosome, pfam_list, taxid, uniprot
    except IndexError:
        return None


def get_signature_for_fasta_argument_helper(gene_name, argument, default_string=""):
    if argument is None:
        gene_argument = default_string
    elif type(argument) is str:
        gene_argument = argument
    elif type(argument) is dict:
        try:
            gene_argument = argument[gene_name]
        except KeyError:
            gene_argument = default_string
    else:
        gene_argument = default_string
    return gene_argument


def get_signature_for_fasta(fasta_file, hmm_hits, ep_hit_dict, type_assignment=None, locus_name=None):
    gene_list = []  # collection information for locus visualization
    sig_list = []  # one line space separated locus information summary with all genes and HMM matches
    total_protein_length = 0
    with open(fasta_file, "r") as handle:
        for record in seq_parse(handle, "fasta"):
            protein_length = len(str(record.seq))
            total_protein_length = total_protein_length + protein_length
            if fa_header_to_fields(record.id) is None:
                continue
            gene_name, gene_pos, directionality, chromosome, pfam_list, taxid, uniprot = fa_header_to_fields(record.id)
            gene_type_assignment = \
                get_signature_for_fasta_argument_helper(gene_name, type_assignment, default_string="unassigned")
            gene_locus_name = get_signature_for_fasta_argument_helper(gene_name, locus_name, default_string=gene_name)
            if pfam_list == "pf00936" or pfam_list == "pf03319":  # don't show pfams for shell proteins
                pfam_list = ""
            try:
                protein_info = hmm_hits[record.id][0].split(".")[0] + "(" + pfam_list + ")"
                bmctype = hmm_hits[record.id][0].split(".")[1]  # extract bmctype from HMM name string
                sig_list.append((gene_locus_name, gene_type_assignment, chromosome, gene_pos,
                                 "-" + hmm_hits[record.id][0].split(".")[0] + "-" + bmctype))
            except KeyError:
                protein_info = "-" + pfam_list
                bmctype = "-"
                if pfam_list == "":
                    pfam_list = "#"
                sig_list.append((gene_locus_name, gene_type_assignment, chromosome, gene_pos, "-" + pfam_list + "-#"))
            try:
                ep_match_from = ep_hit_dict[record.id][2]
                ep_match_to = ep_hit_dict[record.id][3]
                if ep_match_from < 25:
                    ep_match_string = "EP_NTERM"
                elif (protein_length - ep_match_to) < 25:
                    ep_match_string = "EP_CTERM"
                else:
                    ep_match_string = ""
                ep_str = ep_match_string
            except KeyError:
                ep_str = ""
            try:
                organism = "".join(org_dict[gene_name])
            except KeyError:
                organism = "unknown"
            gene_list.append((gene_locus_name, gene_type_assignment, organism, taxid, chromosome, gene_name, gene_pos,
                              uniprot, protein_length, protein_info, bmctype, directionality, ep_str))
    return gene_list, sig_list, total_protein_length


def main(fasta_path, hmm_hits_path, ep_hits_path, type_assignment=None, locus_name=None):
    if os.path.isfile(hmm_hits_path):
        hmm_hits = read_hmmscan_tab_output(hmm_hits_path, 'hmmsearch')
    else:
        ret = run_hmm(fasta_path, ref_hmm, hmm_hits_path, "hmmsearch", additional_cmd="-T 40")
        if ret != 0:
            raise SystemError
        hmm_hits = read_hmmscan_tab_output(hmm_hits_path, 'hmmsearch')
    if os.path.isfile(ep_hits_path):
        ep_hit_dict = read_hmmscan_tab_output(ep_hits_path, "hmmscan", parse_type='hmmer3-text', evalue_cutoff=5E-5)
    else:
        ret = run_hmm(fasta_path, ep_ref_hmm, ep_hits_path + '.tab', "hmmscan", log_path=ep_hits_path)
        if ret != 0:
            raise SystemError
        ep_hit_dict = read_hmmscan_tab_output(ep_hits_path, "hmmscan", parse_type='hmmer3-text', evalue_cutoff=5E-5)
    gene_list, sig_list, total_protein_length = \
        get_signature_for_fasta(fasta_path, hmm_hits, ep_hit_dict, type_assignment, locus_name)
    sig_list.sort(key=lambda x: x[2])
    gene_list.sort(key=lambda x: x[4])
    return gene_list, sig_list, total_protein_length


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest="fasta_path", required=True)
    parser.add_argument('-o', '--output', dest="output_folder")
    parser.add_argument('-m', '--hmm', dest="hmm_hits_path")
    parser.add_argument('-e', '--ep', dest="ep_hits_path")
    # BMC Type
    parser.add_argument('-t', '--type', dest="type_assignment")
    # Locus Name
    parser.add_argument('-l', '--locus', dest="locus_name")
    parser.add_argument('-f', '--format', dest="format", choices=["new", "og"], default="og")

    args = parser.parse_args()

    fasta_name, fasta_ext = os.path.splitext(os.path.basename(args.fasta_path))
    if args.output_folder is None:
        output_folder = os.path.dirname(args.fasta_path)
    else:
        output_folder = args.output_folder
    if args.hmm_hits_path is None:
        hmm_hits_path = os.path.join(output_folder, fasta_name + '_tbl.tmp')
    else:
        hmm_hits_path = args.hmm_hits_path
    if args.ep_hits_path is None:
        ep_hits_path = os.path.join(output_folder, fasta_name + '_ep_log')
    else:
        ep_hits_path = args.ep_hits_path

    gene_list, sig_list, total_protein_length = main(args.fasta_path, hmm_hits_path, ep_hits_path,
                                                     type_assignment=args.type_assignment, locus_name=args.locus_name)

    gene_list_output = os.path.join(output_folder, fasta_name + '.gene_list.tsv')
    signature_list_output = os.path.join(output_folder, fasta_name + '.sig.txt')
    with open(signature_list_output, 'w') as sop:
        sig_dict = {}
        for idx, sig in enumerate(sig_list):
            gene_locus_name, gene_type_assignment, chromosome, gene_pos, sig_attr = sig
            key = gene_locus_name + ' ' + gene_type_assignment
            if args.format == "og":
                try:
                    sig_dict[key].append(gene_pos.split('-')[0] + sig_attr)
                except KeyError:
                    sig_dict.setdefault(key, [gene_pos.split('-')[0] + sig_attr])
            else:
                key += chromosome
                try:
                    sig_dict[key].append(gene_pos + sig_attr)
                except KeyError:
                    sig_dict.setdefault(key, [gene_pos + sig_attr])
        for gene in sig_dict:
            sop.write(gene + ' ' + ' '.join(sig_dict[gene]) + '\n')
    with open(gene_list_output, 'w') as gop:
        for g in gene_list:
            if args.format == "og":
                gop.write("\t".join([str(i) for i in g[5:]]) + '\n')
            else:
                gop.write("\t".join([str(i) for i in g]) + '\n')