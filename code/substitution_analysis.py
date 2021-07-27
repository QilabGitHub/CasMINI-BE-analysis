import glob
import os
import pandas as pd
from Bio.Seq import Seq
import math
from collections import defaultdict
import re


MAXWINGLENGTH = 0
EDGE_CUTOFF = 5
OUTPUT_DIR = f"output_{MAXWINGLENGTH}bp"

MIN_AVERAGE_READ_QUALITY_CUTOFF = 30
CONVERSION_NUC_FROM = "A"
CONVERSION_NUC_TO = "G"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

docker_stem = f"docker run -v {os.path.dirname(os.path.realpath(__file__))}:/DATA -w /DATA -i pinellolab/crispresso2"
boilerplate = f"--min_average_read_quality {MIN_AVERAGE_READ_QUALITY_CUTOFF} -o {OUTPUT_DIR}/ --exclude_bp_from_left {EDGE_CUTOFF} --exclude_bp_from_right {EDGE_CUTOFF} --base_editor_output --conversion_nuc_from {CONVERSION_NUC_FROM} --conversion_nuc_to {CONVERSION_NUC_TO} --no_rerun"


def rc(guide):
    return str(Seq(guide).reverse_complement())


def regenerate_log(sample_metadata):
    all_data = []

    for dat in sample_metadata:
        sample = dat["sample"].strip().upper()
        amplicon = dat["amplicon"].strip().upper()
        guide = dat["guide"].strip().upper()

        # Unzip the full allele frequencies tables
        fname = f"{OUTPUT_DIR}/CRISPResso_on_{sample}/Alleles_frequency_table.txt"
        if not os.path.exists(fname): os.system(f"unzip {OUTPUT_DIR}/CRISPResso_on_{sample}/Alleles_frequency_table.zip -d {OUTPUT_DIR}/CRISPResso_on_{sample}")
        if not os.path.exists(fname): continue

        sample_data = defaultdict(int)
        tot_reads = 0
        # Now process the full allele frequencies tables
        with open(fname, "r") as ropen:
            for idx, line in enumerate(ropen):
                if idx == 0: continue
                sl = line.split()

                new_sequence = sl[0]
                ref_sequence = sl[1]
                num_reads = int(sl[-2])

                if rc(guide) in amplicon:
                    new_sequence = rc(new_sequence)
                    ref_sequence = rc(ref_sequence)

                tot_reads += num_reads

                # Delete all insertion positions
                zipped = list(zip(ref_sequence, new_sequence))
                zipped = [(r, n) for r, n in zipped if r != "-"]

                ref_reassembled = ''.join([r for r, n in zipped])
                guide_idx = ref_reassembled.find(guide)
                left_idx = guide_idx - 4
                left_idx = guide_idx
                right_idx = guide_idx + len(guide)

                zipped = zipped[left_idx: right_idx]
                ref_reassembled = ''.join([r for r, n in zipped])
                new_reassembled = ''.join([n for r, n in zipped])

                dat = set()
                for idx, token in enumerate(list(zip(ref_reassembled, new_reassembled))):
                    ref, new = token
                    
                    if ref == CONVERSION_NUC_FROM and new == CONVERSION_NUC_TO:
                        dat.add(f"{ref}{idx + 1}->{new}")
                        dat.add(f"{ref}->{new}")
                
                for substitution_type in dat:
                    sample_data[substitution_type] += num_reads


        final_dat = {"sample": sample, "total_reads": tot_reads}
        for modification_type, read_count in sample_data.items():
            final_dat[modification_type] = read_count
        all_data.append(final_dat)

    if len(all_data) == 0: return
    df = pd.DataFrame.from_records(all_data)
    df = df.sort_values(by="sample")
    new_columns = ["sample", "total_reads", f"{CONVERSION_NUC_FROM}->{CONVERSION_NUC_TO}"] + [f"{CONVERSION_NUC_FROM}{idx}->{CONVERSION_NUC_TO}" for idx in range(1, 29)]
    new_columns = [col for col in new_columns if col in df.columns]
    df = df[new_columns]
    df.to_excel(f"{OUTPUT_DIR}/substitution_readcounts_summary.xlsx", index=False)


def main():
    sample_metadata = pd.read_csv("experiment_overview.csv").to_dict("records")

    for dat in sample_metadata:
        regenerate_log(sample_metadata)

        sample = dat["sample"].strip().upper()
        amplicon = dat["amplicon"].strip().upper()
        guide = dat["guide"].strip().upper()

        if len(amplicon) == 0 or len(guide) == 0: continue

        print(f"Analyzing sample {sample}")

        if not (guide in amplicon or rc(guide) in amplicon):
            print(f"ERROR: Couldn't find guide in amplicon for sample {sample}")
            continue

        if guide in amplicon:
            left = amplicon.index(guide)
            left = left - MAXWINGLENGTH

            right = amplicon.index(guide) + len(guide)
            right = right + MAXWINGLENGTH

            if EDGE_CUTOFF <= left and right <= len(amplicon) - EDGE_CUTOFF:
                for fname in glob.glob(f"raw_data/{sample}_*/*.fastq.gz"):
                    command = f"{docker_stem} CRISPResso --fastq_r1 {fname} --amplicon_seq {amplicon} --guide_seq {guide} -w {12 + MAXWINGLENGTH} -wc -12 --plot_window_size {12 + MAXWINGLENGTH} -n {sample} {boilerplate}"
                    print(command)
                    os.system(command)
            else:
                left = max(left, EDGE_CUTOFF)
                right = min(right, len(amplicon) - EDGE_CUTOFF)

                length = right - left
                wing_length = int(math.floor(length / 2))
                center = left + wing_length
                center = center - (amplicon.index(guide) + len(guide))

                for fname in glob.glob(f"raw_data/{sample}_*/*.fastq.gz"):
                    command = f"{docker_stem} CRISPResso --fastq_r1 {fname} --amplicon_seq {amplicon} --guide_seq {guide} -w {wing_length} -wc {center} --plot_window_size {wing_length} -n {sample} {boilerplate}"
                    print(command)
                    os.system(command)
        elif rc(guide) in amplicon:
            rc_guide = rc(guide)

            left = amplicon.index(rc_guide)
            left = left - MAXWINGLENGTH

            right = amplicon.index(rc_guide) + len(rc_guide)
            right = right + MAXWINGLENGTH

            if EDGE_CUTOFF <= left and right <= len(amplicon) - EDGE_CUTOFF:
                for fname in glob.glob(f"raw_data/{sample}_*/*.fastq.gz"):
                    command = f"{docker_stem} CRISPResso --fastq_r1 {fname} --amplicon_seq {amplicon} --guide_seq {guide} -w {12 + MAXWINGLENGTH} -wc -12 --plot_window_size {12 + MAXWINGLENGTH} -n {sample} {boilerplate}"
                    print(command)
                    os.system(command)
            else:
                left = max(left, EDGE_CUTOFF)
                right = min(right, len(amplicon) - EDGE_CUTOFF)

                length = right - left
                wing_length = int(math.floor(length / 2))
                center = left + wing_length
                center = amplicon.index(rc_guide) - center

                for fname in glob.glob(f"raw_data/{sample}_*/*.fastq.gz"):
                    command = f"{docker_stem} CRISPResso --fastq_r1 {fname} --amplicon_seq {amplicon} --guide_seq {guide} -w {wing_length} -wc {center} --plot_window_size {wing_length} -n {sample} {boilerplate}"
                    print(command)
                    os.system(command)


if __name__=="__main__":
    main()
