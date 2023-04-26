#!/usr/bin/env python3
"""Module that splits the input bam file into 3 bam files based on one input region."""


import argparse
import os


import pandas as pd
import pysam

def get_primary_lists(bam, region):
    """
    Generates to QNAME lists of reads whose primary alignment is in/out of the region

    Parameters
    ----------
    bam: pysam.AlignmentFile
    The input bam
    region: str
    The region's name in the fasta file

    Returns
    -------
    in_region_list: list
    a list of read qnames whose first alignment is in the region
    out_region_list: list
    a list of read qnames whose first alignment is out of the region
    """
    total_reads = bam.mapped + bam.unmapped
    temp_df = pd.DataFrame({"read_qname":["na"] * total_reads, "read_status":["na"] * total_reads})

    counter = 0
    for read in bam.fetch(until_eof=True):
        temp_df.iloc[counter, 0] = read.qname
        if not read.is_unmapped and not read.is_secondary:
            if read.reference_name == region:
                temp_df.iloc[counter, 1] = "in"
            else:
                temp_df.iloc[counter, 1] = "out"
        else:
            temp_df.iloc[counter, 1] = "not_interested"
        counter += 1

    in_region_list = list(temp_df[temp_df["read_status"] == "in"]["read_qname"])
    out_region_list = list(temp_df[temp_df["read_status"] == "out"]["read_qname"])

    return in_region_list, out_region_list


def write_split_bams(bam, bam_name, in_list, out_list, region):
    """
    Generate 3 bam files, where all alignments are in the region, resp out of the region, and mixed.

    Parameters
    ----------
    bam: pysam.AlignmentFile
    the input bam
    bam_name: str
    the input file name (for naming the output bams)
    in_list: list
    a list of all the reads qnames whose primary alignments are in the region
    out_list: list
    a list of all the reads qnames whose primary alignments are out of the region
    region: str
    the region's name in the fasta file
    """
    bam_nameind = pysam.IndexedReads(bam)
    bam_nameind.build()
    in_bam_name = bam_name.split(".")[0] + "_in_" + region + ".bam"
    out_bam_name = bam_name.split(".")[0] + "_out_" + region + ".bam"
    mixed_bam_name = bam_name.split(".")[0] + "_mixed_" + region + ".bam"
    # in_bam_sorted = in_bam_name.split(".")[0] + "_sorted" + ".bam"
    # out_bam_sorted = in_bam_name.split(".")[0] + "_sorted" + ".bam"
    # mixed_bam_sorted = in_bam_name.split(".")[0] + "_sorted" + ".bam"

    in_bam = pysam.AlignmentFile(in_bam_name, "wb", template=bam)
    out_bam = pysam.AlignmentFile(out_bam_name, "wb", template=bam)
    mixed_bam = pysam.AlignmentFile(mixed_bam_name, "wb", template=bam)

    write_match_bam(bam_nameind, in_bam, mixed_bam, in_list, region, True)
    write_match_bam(bam_nameind, out_bam, mixed_bam, out_list, region, False)
    # import pdb; pdb.set_trace()
    # pysam.sort(in_bam_name, "-o", in_bam_sorted)
    # pysam.sort(out_bam_name, "-o", out_bam_sorted)
    # pysam.sort(mixed_bam_name, "-o", mixed_bam_sorted)
    # pysam.index(in_bam_sorted)
    # pysam.index(out_bam_sorted)
    # pysam.index(mixed_bam_sorted)
    # os.remove(in_bam_name)
    # os.remove(out_bam_name)
    # os.remove(mixed_bam_name)


def write_match_bam(source_bam, out_bam, mixed_bam, match_list, region, check_in):
    """
    Write reads alignements to a bam file if all alignments are in resp out of the region.

    Parameters
    ----------
    source_bam: pysam.AlignmentFile
    the name indexed original bam file
    out_bam: pysam.AlignmentFile
    the output bam
    mixed_bam: pysam.AlignmentFile
    the output bam if a read is aligned in and out of the region
    match_list: list
    a list of all the reads qnames whose primary alignments are in resp out of the region
    region: str
    the region's name in the fasta file
    check_in: bool
    Wether the reads should have all alignements in resp out of the region
    """
    for qname in match_list:
       try:
           source_bam.find(qname)
       except KeyError:
           pass
       else:
           match_iter = source_bam.find(qname)
           all_ok = True
           for match in match_iter:
               if match.is_secondary:
                   if check_in ^ (match.reference_name != region):
                       continue
                   else:
                       all_ok = False
                       break
           match_store = source_bam.find(qname)
           if all_ok:
               for match in match_store:
                   out_bam.write(match)
           else:
               for match in match_store:
                   mixed_bam.write(match)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Splits a bam file into 3 based on one region.
        (1) Aligned reads to the region with all secondary alignments to the region
        (2) Aligned reads to other regions with all secondary alignments to other regions
        (3) The rest.""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "bam_file",
        nargs="?",
        help="the bam file to split"
    )
    # add bai file? no
    parser.add_argument(
        "region",
        nargs="?",
        help="The region (reference in pysam) as described in the fasta file"
    )
    args = parser.parse_args()
    source_bam = pysam.AlignmentFile(args.bam_file, "rb")
    in_region_list, out_region_list = get_primary_lists(source_bam, args.region)
    write_split_bams(source_bam, args.bam_file, in_region_list, out_region_list, args.region)
