import glob
import os
import pandas as pd
from Bio.Seq import Seq
import math
from collections import defaultdict
import re

from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (20,6)

MAXWINGLENGTH = 15
EDGE_CUTOFF = 5
OUTPUT_DIR = f"output_{MAXWINGLENGTH}bp"

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

docker_stem = f"docker run -v {os.path.dirname(os.path.realpath(__file__))}:/sample_metadata -w /sample_metadata -i pinellolab/crispresso2"
boilerplate = f"--min_average_read_quality 30 -o {OUTPUT_DIR}/ --exclude_bp_from_left {EDGE_CUTOFF} --exclude_bp_from_right {EDGE_CUTOFF} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --no_rerun"


def rc(guide):
    return str(Seq(guide).reverse_complement())


def longest_gap(seq):
    to_return = 0
    cache = 0

    for char in seq:
        if char == "-":
            cache += 1
            to_return = max(cache, to_return)
        else:
            cache = 0

    return to_return


def deletion_pos_plot():
    sample_metadata = pd.read_csv("experiment_overview.csv").to_dict("records")
    for dat in sample_metadata:
        sample = dat["sample"].strip()
        amplicon = dat["amplicon"].strip().upper()
        guide = dat["guide"].strip().upper()

        fnames = glob.glob(f"{OUTPUT_DIR}/*_{sample}/Alleles_frequency_table_around_sgRNA_*.txt")
        if len(fnames) == 0: continue
        fname = fnames[0]

        tot_reads = 0
        deletions = defaultdict(int)

        with open(fname, "r") as ropen:
            for idx, line in enumerate(ropen):
                if idx == 0: continue
                sl = line.split()

                new_sequence = sl[0].strip().upper()
                ref_sequence = sl[1].strip().upper()
                num_reads = int(sl[-2])

                if "N" in new_sequence: continue
                tot_reads += num_reads

                if guide in rc(ref_sequence.replace("-", "")):
                    new_sequence = rc(new_sequence)
                    ref_sequence = rc(ref_sequence)

                zipped = list(zip(ref_sequence, new_sequence))
                zipped = [(ref, new) for ref, new in zipped if ref != "-"]

                zipped_ref = ''.join([ref for ref, new in zipped])
                guide_pos = zipped_ref.find(guide)
                if guide_pos == -1: continue

                for idx, val in enumerate(zipped):
                    ref_char, new_char = val
                    if new_char == "-":
                        deletions[idx - guide_pos + 5] += num_reads


        deletions = {k: v / tot_reads * 100 for k, v in deletions.items()}

        if not os.path.exists("plots/deleted_positions/"): os.makedirs("plots/deleted_positions/")
        plt.figure()
        plt.bar([str(k) for k in sorted(deletions.keys())], [deletions[k] for k in sorted(deletions.keys())])
        plt.xlabel("Base position")
        plt.ylabel(r"% of total reads with a deletion at this position")
        plt.savefig(f"plots/deleted_positions/{sample}_deleted_positions.png")
        plt.savefig(f"plots/deleted_positions/{sample}_deleted_positions.pdf")


def indel_size_plot():
    sample_metadata = pd.read_csv("experiment_overview.csv").to_dict("records")
    sample_names = [dat["sample"] for dat in sample_metadata]

    for sample_name in sample_names:
        fnames = glob.glob(f"{OUTPUT_DIR}/*_{sample_name}/Alleles_frequency_table_around_sgRNA_*.txt")
        if len(fnames) == 0: continue
        fname = fnames[0]

        tot_reads = 0

        insertions = defaultdict(int)
        deletions = defaultdict(int)
        indels = defaultdict(int)
        
        with open(fname, "r") as ropen:
            for idx, line in enumerate(ropen):
                if idx == 0: continue
                sl = line.split()

                new_sequence = sl[0]
                ref_sequence = sl[1]
                num_reads = int(sl[-2])

                if "N" in new_sequence: continue

                tot_reads += num_reads

                insertion_size = longest_gap(ref_sequence)
                deletion_size = -1 * longest_gap(new_sequence)

                if insertion_size != 0:
                    insertions[insertion_size] += num_reads
                    indels[insertion_size] += num_reads
                if deletion_size != 0:
                    deletions[deletion_size] += num_reads
                    indels[deletion_size] += num_reads

        insertions = {k: v / tot_reads * 100 for k, v in insertions.items()}
        deletions = {k: v / tot_reads * 100 for k, v in deletions.items()}
        indels = {k: v / tot_reads * 100 for k, v in indels.items()}

        if not os.path.exists("plots/indel_sizes/"): os.makedirs("plots/indel_sizes/")
        plt.figure()
        plt.bar([str(k) for k in sorted(indels.keys())], [indels[k] for k in sorted(indels.keys())])
        plt.xlabel("Insertion/Deletion size")
        plt.ylabel(r"% of total reads with an indel of this length")
        plt.savefig(f"plots/indel_sizes/{sample_name}_indel_sizes.png")
        plt.savefig(f"plots/indel_sizes/{sample_name}_indel_sizes.pdf")


def regenerate_log(sample_metadata):
    all_sample_data = []
    sample_names = [dat["sample"] for dat in sample_metadata]

    for sample_name in sample_names:
        fnames = glob.glob(f"{OUTPUT_DIR}/*_{sample_name}/Alleles_frequency_table_around_sgRNA_*.txt")
        if len(fnames) == 0: continue
        fname = fnames[0]

        single_sample_data = defaultdict(int)
        tot_reads = 0
        with open(fname, "r") as ropen:
            for idx, line in enumerate(ropen):
                if idx == 0: continue
                sl = line.split()

                new_sequence = sl[0]
                ref_sequence = sl[1]
                num_reads = int(sl[-2])

                if "N" in new_sequence: continue

                tot_reads += num_reads

                dat = set()
                for ref, new in zip(ref_sequence, new_sequence):
                    if ref != new and ref != "-" and new != "-":
                        substitution_type = f"{ref}->{new}"
                        dat.add(substitution_type)
                
                for substitution_type in dat:
                    single_sample_data[substitution_type] += num_reads

                if "-" in new_sequence:
                    single_sample_data["deletions"] += num_reads
                if "-" in ref_sequence:
                    single_sample_data["insertions"] += num_reads
                if "-" in new_sequence or "-" in ref_sequence:
                    single_sample_data["insertions_or_deletions"] += num_reads
                if new_sequence != ref_sequence:
                    single_sample_data["any_modification"] += num_reads


        final_dat = {"sample": sample_name}
        for modification_type, read_count in single_sample_data.items():
            final_dat[modification_type] = read_count / tot_reads * 100
        all_sample_data.append(final_dat)

    if len(all_sample_data) == 0: return
    df = pd.DataFrame.from_records(all_sample_data)
    df = df.sort_values(by="sample")
    new_columns = ["sample"] + [f"{ref}->{new}" for ref in ["A", "T", "C", "G"] for new in ["A", "T", "C", "G"]] + ["insertions", "deletions", "insertions_or_deletions", "any_modification"]
    new_columns = [col for col in new_columns if col in df.columns]
    df = df[new_columns]
    df.to_excel(f"{OUTPUT_DIR}/modification_summary.xlsx", index=False)


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
    
    indel_size_plot()
    deletion_pos_plot()


if __name__=="__main__":
    main()
