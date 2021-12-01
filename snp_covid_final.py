#PIPELINE GALO VERSION JON

#!/usr/bin/env python

#   Standard library imports
import os
import sys
import re
import argparse
import subprocess
import pandas as pd
import numpy as np

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'snp_covid19.py', description = 'Pipeline to do variant calling with Sars-CoV-2 samples')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest="input_dir", type=str, required=True, help="Required. Input dir with the fastq files")
    input_group.add_argument('-r', '--reference_genome', metavar="reference", type=str, required=True, help='File to map against')
    input_group.add_argument('-p', '--primers_bed', dest="primers", required=True, help='File with the primers')
    input_group.add_argument('-a', '--annotation', dest="annotation", required=True, help='File with the annotation files')
    input_group.add_argument('-R', '--run_id', dest="run_id", required=True, help='Run ID')
    input_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')

    output_group = parser.add_argument_group('Output', 'Output parameters')

    output_group.add_argument('-o', '--output', type=str, required=True, help="Required.Output directory to extract all results")

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=16, help='Threads to use')

    arguments = parser.parse_args()

    return arguments

args = get_arguments()

def extract_read_list(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    190615 - Limit only parent folder, not subdirectories
    """
    input_dir = os.path.abspath(input_dir)
    all_fasta = []
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        if root == input_dir: # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_fasta = re.match(r'.*\.f(ast)*[aq](\.gz)*',filename)
                if is_fasta:
                    all_fasta.append(filename)
    all_fasta = sorted(all_fasta)
    if len(all_fasta) % 2 == 0:
        for index, fasta_file in enumerate(all_fasta):
            if index % 2 == 0:
                r1_list.append(fasta_file)
            elif index % 1 == 0:
                r2_list.append(fasta_file)          
    else:
        print('ERROR: The number of fastq sequence are not paired')
        
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)
    
    return r1_list, r2_list

def extract_sample(R1_file, R2_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R1 = os.path.basename(R1_file)
    basename_R2 = os.path.basename(R2_file)

    sample_name_R = os.path.commonprefix([basename_R1, basename_R2])
  
    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    
    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip("_")
    else:
        sample_name = sample_name_R

    return sample_name

def file_to_list(file_name):
    list_F = []
    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, "r") as f:
        for line in f:
            list_F.append(line.strip())
    return list_F

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

out_consensus_dir = os.path.join(args.output, "Consensus")
out_variants_dir = os.path.join(args.output, "Variants")
out_reports_dir = os.path.join(args.output, "Reports")
out_log_dir = os.path.join(args.output, "Log")

out_bam_dir = os.path.join(args.output, "Bam")
out_trim_dir = os.path.join(args.output, "Trimmed")
out_phyl_dir = os.path.join(args.output, "Phylogeny")

out_line_dir = os.path.join(args.output, "Lineages")
out_compare_dir = os.path.join(args.output, "Compare")

out_emptyvariants_dir = os.path.join(args.output, "Empty_variants")
out_onlysnps_variants = os.path.join(args.output, "Variants_snps")

#Obtain all R1 and R2 from folder
r1, r2 = extract_read_list(args.input_dir)

#Check if there are samples to filter
sample_list_F = []
if args.sample_list == None:
    print("\n" + "No samples to filter")
    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)
        sample_list_F.append(sample)
else:
    print("samples will be filtered")
    sample_list_F = file_to_list(args.sample_list)
print("\n%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
    if sample in sample_list_F:

        sample_number = str(sample_list_F.index(sample) + 1)
        sample_total = str(len(sample_list_F))

        
        
        args.r1_file = r1_file
        args.r2_file = r2_file


        print("\n" + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")")
        # previous trimming

        check_create_dir(out_trim_dir)

        r1 = os.path.abspath(args.r1_file)
        r2 = os.path.abspath(args.r2_file)

        r1_trim = sample + "_R1.trim.fastq.gz"
        out_r1_trim = os.path.join(out_trim_dir, r1_trim)
        r2_trim = sample + "_R2.trim.fastq.gz"
        out_r2_trim = os.path.join(out_trim_dir, r2_trim)

        json = sample + "_fastp_trim.json"
        json_path = os.path.join(out_reports_dir, json)
        html = sample + "_fastp_trim.html"
        html_path = os.path.join(out_reports_dir, html)

        adapter = "/home/laura/COVID/adaptors/IlluminaAdaptors.fasta"

        cmd_trim = ["fastp --in1 " + r1 + " --in2 " + r2 + " --out1 " + out_r1_trim + " --out2 " + out_r2_trim + " --detect_adapter_for_pe --adapter_fasta " + adapter + " --cut_tail --cut_window_size 10 --cut_mean_quality 30 --thread 4 --max_len1 250 --max_len2 250 --json " + json_path + " --html " + html_path]
        subprocess.run(cmd_trim, shell=True)

        # MAPPING

        #bwa mem  -Y -M -R '@RG\tID:.\tSM:.' -t ${THREADS} ${REF_GENOME} ${FQ1} ${FQ2} 2>>${RUN_ID}.log | samtools sort | samtools view -F 4 -b -@ ${THREADS} - > ${SAMPLE}.sort.bam

        check_create_dir(out_bam_dir)
        #os.mkdir(out_log_dir)

        out_map_name = sample + ".sort.bam"
        out_map_file = os.path.join(out_bam_dir, out_map_name)

        reference = os.path.abspath(args.reference_genome)
        r1_trimmed = sample + "_R1.trim.fastq.gz"
        out_r1_trimmed = os.path.join(out_trim_dir, r1_trimmed)
        r2_trimmed = sample + "_R2.trim.fastq.gz"
        out_r2_trimmed = os.path.join(out_trim_dir, r2_trimmed)
        #r1 = os.path.abspath(args.r1_file)
        #r2 = os.path.abspath(args.r2_file)
        
        run = args.run_id
        run_file = run + ".log"
        run_path = os.path.join(args.output, run_file)

        cmd_map = ['bwa mem -Y -M -R "@RG\\tID:.\\tSM:." -t 4 ' + reference + " " + out_r1_trimmed + " " + out_r2_trimmed + ' | samtools sort | samtools view -F 4 -b -@ 4 -o ' + out_map_file]
        subprocess.run(cmd_map, shell=True)

        # PRME MASKIN RIMMING

        #ivar trim -i ${SAMPLE}.sort.bam -p ${SAMPLE}.trim -m ${MINTRIM_LENGTH} -q ${MINQUALTRIM} -s ${WINDOW} -b ${PRIMERS_BED} &>>${RUN_ID}.log

        in_trim = sample + ".sort.bam"
        input_trim = os.path.join(out_bam_dir, in_trim)

        out_trim = sample + ".trim"
        output_trim = os.path.join(out_bam_dir, out_trim)

        primer = os.path.abspath(args.primers)

        cmd_trim = ["/home/laura/Downloads/ivar-1.2/src/ivar trim -i " + input_trim + " -p " + output_trim + " -m 30 -q 20 -s 4 -b " + primer + " -e"]
        subprocess.run(cmd_trim, shell=True)

        # resort bam, index and depth

        #samtools sort -o ${SAMPLE}.trim.sort.bam ${SAMPLE}.trim.bam 1>&2 2>>${RUN_ID}.log
        #samtools index ${SAMPLE}.trim.sort.bam &>>${RUN_ID}.log
        #samtools depth -aa --reference ${REF_GENOME} ${SAMPLE}.trim.sort.bam > ${SAMPLE}.depth 2>>${RUN_ID}.log

        in_sort = sample + ".trim.bam"
        input_sort = os.path.join(out_bam_dir, in_sort)

        out_sort = sample + ".trim.sort.bam"
        output_sort = os.path.join(out_bam_dir, out_sort)

        cmd_sort = ["samtools sort -o " + output_sort + " " + input_sort]
        subprocess.run(cmd_sort, shell=True)

        in_index = sample + ".trim.sort.bam"
        input_index = os.path.join(out_bam_dir, in_index)

        cmd_index = ["samtools index " + input_index]
        subprocess.run(cmd_index, shell=True)

        check_create_dir(out_reports_dir)

        in_depth = sample + ".trim.sort.bam"
        input_depth = os.path.join(out_bam_dir, in_depth)
        out_depth = sample + ".depth"
        output_depth = os.path.join(out_reports_dir, out_depth)
        reference = os.path.abspath(args.reference_genome)

        cmd_depth = ["samtools depth -aa --reference " + reference + " " + input_depth + " > " + output_depth]
        subprocess.run(cmd_depth, shell=True)

        ### IVAR CNSESUS

        #samtools mpileup -aa -A -d 0 -B -Q 0 ${SAMPLE}.trim.sort.bam | ivar consensus -p ${SAMPLE} -q ${MINQUAL_CONS} -t ${MINFREQ_CONS} -m ${MINDEPTH_CONS} -n ${AMBIGUOUS_CHAR} &>>${RUN_ID}.log

        check_create_dir(out_consensus_dir)

        in_consensus = sample + ".trim.sort.bam"
        input_consensus = os.path.join(out_bam_dir, in_consensus)

        cmd_consensus = ["samtools mpileup -aa -A -d 0 -B -Q 0 " + input_consensus + " | /home/laura/Downloads/ivar-1.2/src/ivar consensus -p " + sample + " -q 20 -t 0.8 -m 30 -n N"]
        subprocess.run(cmd_consensus, shell=True)

        out_consensus = sample + ".fa"
        output_consensus = os.path.join(args.output, out_consensus)
        cmd_mv_one = ["mv " + output_consensus + " " + out_consensus_dir]
        subprocess.run(cmd_mv_one, shell=True)

        out_qual_consensus = sample + ".qual.txt"
        output_qual_consensus = os.path.join(args.output, out_qual_consensus)
        cmd_mv_two = ["mv " + output_qual_consensus + " " + out_consensus_dir]
        subprocess.run(cmd_mv_two, shell=True)

        ### ivar variant calling 

        check_create_dir(out_variants_dir)

        #samtools mpileup -aa -A -d 0 -B -Q 0 --reference ${REF_GENOME}  ${SAMPLE}.trim.sort.bam | ivar variants -p ${SAMPLE} -q ${MINQUAL_VAR} -t ${MINFREQ_VAR}
        
        in_variant = sample + ".trim.sort.bam"
        input_variant = os.path.join(out_bam_dir, in_variant)

        reference = os.path.abspath(args.reference_genome)
        annotation = os.path.abspath(args.annotation)
        
        cmd_variant = ["samtools mpileup -aa -A -d 0 -B -Q 0 --reference " + reference + " " + input_variant + " | /home/laura/Downloads/ivar-1.2/src/ivar variants -p " + sample + " -q 20 -t 0.05 -m 20 -r " + reference + " -g " + annotation]
        subprocess.run(cmd_variant, shell=True)

        out_var = sample + ".tsv"
        output_var = os.path.join(args.output, out_var)
        cmd_mv_three = ["mv " + output_var + " " + out_variants_dir]
        subprocess.run(cmd_mv_three, shell=True)

        ### calculate depth (.depth is done before)

        depth = sample + ".depth"
        in_depth = os.path.join(out_reports_dir, depth)

        wdepth = sample + ".wdepth"
        out_depth = os.path.join(out_reports_dir, wdepth)

        primer = os.path.abspath(args.primers)

        cmd_depth_one = ["python /home/laura/COVID/coverage_files/COVID-window-coverage.py -i " + in_depth + " -p " + primer + " -w 100 -H -s " + sample + " | sort -nk2,2 > " + out_depth]
        subprocess.run(cmd_depth_one, shell=True)

        first_variant = sample + ".trim.sort.bam"
        first_input_variant = os.path.join(out_bam_dir, first_variant)

        bamstat = sample + ".bamstats"
        bamstat_path = os.path.join(out_reports_dir, bamstat)
        cmd_flagstat = ["samtools flagstat " + first_input_variant + " > " + bamstat_path]
        subprocess.run(cmd_flagstat, shell=True)

        depth_two = sample + ".wdepth"
        out_depth_two = os.path.join(out_reports_dir, depth_two)

        consensus = sample + ".fa"
        consensus_out = os.path.join(out_consensus_dir, consensus)

        report_qc = sample + ".depthQC"
        out_report_qc = os.path.join(out_reports_dir, report_qc)

        global_summary = "QC_summary." + args.run_id + ".txt"
        global_summary_path = os.path.join(out_reports_dir, global_summary)

        cmd_depth_two = ["python /home/laura/COVID/coverage_files/COVID-depth-report.py -d " + out_depth_two + " -c " + consensus_out + " -H -s " + sample + " --min-depth 400 -o " + out_report_qc + " >> " + global_summary_path]
        subprocess.run(cmd_depth_two, shell=True)

        # This file does not exist anymore, but we keep it. Nevertheless it superposes the final one, we must solve it

        depth_two = args.run_id + ".depthQC.txt"
        depth_two_out = os.path.join(args.output, depth_two)

        cmd_one = ["cat " + out_report_qc +  " |  awk 'NR == 1 || $1 !=\"Sample\"' >> " +  depth_two_out]
        subprocess.run(cmd_one, shell=True)

        move = ["mv " + depth_two_out + " " + out_reports_dir]
        subprocess.run(move, shell=True)

## Move empty variants to another folder just to know, because it does not affect to the compare and other steps, as we initially thought

check_create_dir(out_emptyvariants_dir)

for r,d,f in os.walk(out_variants_dir):
    for file in f:
        if file.endswith(".tsv"):
            whole_name = os.path.join(r, file)
            read = pd.read_csv(whole_name, sep="\t")
            if read.shape == (0, 19):
                print("This file is empty and will be removed from the variants analysis")
                cmd_copy = ["cp " + whole_name + " " + out_emptyvariants_dir]
                subprocess.run(cmd_copy, shell=True)
            else:
                pass

## Group coverage stats

def calculate_cov_stats(file_cov):
    df = pd.read_csv(file_cov, sep="\t", names=["#CHROM", "POS", "COV" ])
    unmmaped_pos = len(df.POS[df.COV == 0].tolist())
    pos_0_10 = len(df.POS[(df.COV > 0) & (df.COV <= 10)].tolist())
    pos_10_20 = len(df.POS[(df.COV > 10) & (df.COV <= 20)].tolist())
    pos_high20 = len(df.POS[(df.COV > 20)].tolist())
    pos_high30 = len(df.POS[(df.COV > 30)].tolist())
    pos_high500 = len(df.POS[(df.COV >= 500)].tolist())
    pos_high800 = len(df.POS[(df.COV >= 800)].tolist())
    pos_high1200 = len(df.POS[(df.COV >= 1200)].tolist())
    total_pos = df.shape[0]
    unmmaped_prop = "%.2f" % ((unmmaped_pos/total_pos)*100)
    prop_0_10 = "%.2f" % ((pos_0_10/total_pos)*100)
    prop_10_20 = "%.2f" % ((pos_10_20/total_pos)*100)
    prop_high20 = "%.2f" % ((pos_high20/total_pos)*100)
    prop_high30 = "%.2f" % ((pos_high30/total_pos)*100)
    prop_high500 = "%.2f" % ((pos_high500/total_pos)*100)
    prop_high800 = "%.2f" % ((pos_high800/total_pos)*100)
    prop_high1200 = "%.2f" % ((pos_high1200/total_pos)*100)
    
    mean_cov = "%.2f" % (df.COV.mean())
    
    return mean_cov, unmmaped_prop, prop_0_10, prop_10_20, prop_high20, prop_high30, prop_high500, prop_high800, prop_high1200

directory_path = os.path.abspath(out_reports_dir)

output_file_name = args.run_id + ".coverage.tab"
output_file = os.path.join(directory_path,output_file_name)

with open(output_file, "w") as outfile:
        outfile.write("#SAMPLE" + "\t" + "MEAN_COV" + "\t" + "UNMMAPED_PROP" + "\t" + "COV1-10X" + "\t" + "COV10-20X" + "\t" + "COV>20X" + "\t" + "COV>30X" + "\t" + "COV>500X" + "\t" + "COV>800X" + "\t" + "COV>1200" + "\n")
        for root, _, files in os.walk(directory_path):
            for name in files:
                filename = os.path.join(root, name)
                file_name_cov = os.path.basename(filename)
                sample = file_name_cov.split(".")[0]
                if filename.endswith(".depth") and (os.path.getsize(filename) > 0):
                    coverage_stats = calculate_cov_stats(filename)
                    mean_cov = coverage_stats[0]
                    unmmaped_prop = coverage_stats[1]
                    outfile.write(sample + "\t" + ("\t").join(coverage_stats) + "\n")

## Mutation

mut_file = args.run_id + ".mutation.txt"
mutation_file = os.path.join(out_variants_dir, mut_file)

with open(mutation_file, "w") as fout:
    for r, d, f in os.walk(out_variants_dir):
        for file in f:
            if file.endswith(".tsv"):
                name = file.split(".")[0]
                path = os.path.join(r, file)
                prueba = pd.read_csv(path, sep="\t")
                for index, row in prueba.iterrows():
                    if (str(row["REF"]) + str(row["POS"]) + str(row["ALT"])) == "A23403G":
                        fout.write("There's D614G mutation in " + name + "\n")
                    else:
                        pass


## Lineage

# this works, but im not sure with the positions. i've 3 papers reference, but i must contrast.
# also, redirect the results to a txt

# legacy
"""
check_create_dir(out_line_dir)

res = args.run_id + ".lineages.txt"
res_file = os.path.join(out_line_dir, res)

with open(res_file, "w") as fout:
    for r, d, f in os.walk(out_variants_dir):
        for file in f:
            if file.endswith(".tsv"):
                name = file.split(".")[0]
                path = os.path.join(r, file)
                prueba = pd.read_csv(path, sep="\t")
                pos_prueba = list(prueba["POS"])
                if 28144 and 8782 and 14408 and 23403 and 26144 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje V")
                    fout.write(name + " tiene SNPs marcadores de linaje V" + "\n")
                elif 14408 and 23403 and 28144 and 8782 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje S") 
                    fout.write(name + " tiene SNPs marcadores de linaje S" + "\n")
                elif 14408 and 23403 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje G")
                    fout.write(name + " tiene SNPs marcadores de linaje G" + "\n")
                else:
                    print(name + " no tiene SNPs marcadores de linaje")
                    fout.write(name + " no tiene SNPs marcadores de linaje" + "\n")
"""

check_create_dir(out_line_dir)

res = args.run_id + ".lineages.txt"
res_file = os.path.join(out_line_dir, res)

with open(res_file, "w") as fout:
    for r, d, f in os.walk(out_variants_dir):
        for file in f:
            if file.endswith(".tsv"):
                name = file.split(".")[0]
                path = os.path.join(r, file)
                prueba = pd.read_csv(path, sep="\t")
                pos_prueba = list(prueba["POS"])
                if 3037 and 14408 and 23403 in pos_prueba: #quitamos 8782 a pesar de que pensamos que debe estar
                    print(name + " tiene SNPs marcadores de linaje 20A")
                    fout.write(name + " tiene SNPs marcadores de linaje 20A" + "\n")
                elif 28881 and 28882 and 28883 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje 20B") 
                    fout.write(name + " tiene SNPs marcadores de linaje 20B" + "\n")
                elif 1059 and 25563 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje 20C")
                    fout.write(name + " tiene SNPs marcadores de linaje 20C" + "\n")
                elif 8782 and 28144 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje 19A")
                    fout.write(name + " tiene SNPs marcadores de linaje 19A" + "\n")
                #elif 8782 in pos_prueba:
                    #print(name + " tiene SNPs marcadores de linaje 19A")
                    #fout.write(name + " tiene SNPs marcadores de linaje 19A" + "\n")
                elif 11083 and 26144 in pos_prueba:
                    print(name + " tiene SNPs marcadores de linaje 19B")
                    fout.write(name + " tiene SNPs marcadores de linaje 19B" + "\n")
                else:
                    print(name + " no tiene SNPs marcadores de linaje conocidos")
                    fout.write(name + " no tiene SNPs marcadores de linaje" + "\n")


## Compare

check_create_dir(out_compare_dir)

# 3 parts: presence matrix, pairwise and dendogram

#1 - Presence matrix

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb

final_ddbb = blank_database()
sample_filter_list = []
all_samples = 0
new_samples = 0
for r,d,f in os.walk(out_variants_dir):
    for file in f:
        if file.endswith(".tsv"):
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            file_name = file.split(".")[0]
            print(file_name)
            file_path = os.path.join(r, file)
            print(file_path)
            read_file = pd.read_csv(file_path, sep="\t") #read the file
            file_pass = read_file[read_file["PASS"] == True] ###
            final_file = file_pass[~file_pass.ALT.str.startswith("-")&~file_pass.ALT.str.startswith("+")] ###
            for position in final_file["POS"].unique(): #access to the snps of each sample
                if position not in final_ddbb["Position"].values:
                    positions_added.append(int(position))
                    
                    new_row = len(final_ddbb.index)
                    final_ddbb.loc[new_row, "Position"] = int(position)
                    final_ddbb.loc[new_row, "Samples"] = file_name
                    final_ddbb.loc[new_row, "N"] = int(1)
                    final_ddbb.loc[new_row, file_name] = str(1)
                else:
                    positions_shared.append(int(position))
                    
                    index_position = final_ddbb.index[final_ddbb["Position"] == int(position)][0]
                    
                    number_samples_with_position = final_ddbb.loc[index_position, "N"]
                    names_samples_with_position = final_ddbb.loc[index_position, "Samples"]
                    new_names_samples = names_samples_with_position + "," + file_name
                    
                    final_ddbb.loc[index_position, "N"] = number_samples_with_position + 1
                    final_ddbb.loc[index_position, "Samples"] = new_names_samples
                    final_ddbb.loc[index_position, file_name] = str(1)
                
final_ddbb = final_ddbb.fillna(0)      
final_ddbb["Position"] = final_ddbb["Position"].astype(int)
final_ddbb["N"] = final_ddbb["N"].astype(int)

file_presence = args.run_id + ".presence.tsv"
file_presence_def = os.path.join(out_compare_dir, file_presence)

final_ddbb.to_csv(file_presence_def, sep="\t", index=False)

#2 - Pairwise

from sklearn.metrics import pairwise_distances, accuracy_score

def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = accuracy_score(df[sample1], df[sample2]) #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns:
            for sample2 in dataframe.iloc[:,3:].columns:
                if sample1 != sample2:
                    snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                    line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                    f.write(line_distance)

def import_to_pandas(file_table, header=False, sep='\t'):
    if header == False:
        #exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        #Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)
    
    return dataframe

file_presence_again = args.run_id + ".presence.tsv"
file_presence_again_def = os.path.join(out_compare_dir, file_presence_again) # instead of using previous one, we upload again, as in jupyter there was an error if we don't upload de output file

presence_ddbb = import_to_pandas(file_presence_again_def, header=True)

file_pairwise = args.run_id + ".snps.pairwise.tsv"
file_pairwise_path = os.path.join(out_compare_dir, file_pairwise)

snp_distance_pairwise(presence_ddbb, file_pairwise_path)

#3 - Dendrogram

import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt

def dendogram_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average') #method='single'

    plt.rcParams['lines.linewidth'] = 8 #Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10 #Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30}) #Increase x tick label size
    #plt.tick_params(labelsize=30)
    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending', show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    
    plt.savefig(output_file, format="png")

png_file = args.run_id + ".snp.dendrogram.png"
png_dend_file = os.path.join(out_compare_dir, png_file)

dendogram_dataframe(presence_ddbb, png_dend_file)

## Edit variant files to quit indels and put the final file in another folder

check_create_dir(out_onlysnps_variants)

for r, d, f in os.walk(out_variants_dir):
    for file in f:
        new = os.path.join(r, file)
        if new.endswith("tsv"):
            one = new.split("/")[-1]
            two = one.split(".")[0]
            #print(two)
            read_file = pd.read_csv(new, sep="\t")
            file_pass = read_file[read_file["PASS"] == True] ###
            final_file = file_pass[~file_pass.ALT.str.startswith("-")&~file_pass.ALT.str.startswith("+")]
            final_name = two + ".tsv"
            path_dest = out_onlysnps_variants
            final_path = os.path.join(path_dest, final_name)
            final_file.to_csv(final_path, sep="\t", index=False)

# We move it to the bottom because this step lasts more, so in the meantime we could analize the rest of the results
## Phylogeny and trees 

# legacy

"""
for root, dirs, file in os.walk(out_consensus_dir):
    for f in file:
        if f.endswith(".fa"):
            new = os.path.join(root, f)
            with open(new, "r+") as f:
                for line in f:
                    if not line.startswith(">"):
                        new_line = line.rstrip()
                        #print(len(new_line))
                        content_n = "N" * 29903
                        #print(len(content_n))
                        if new_line == content_n:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
"""

"""
for root, dirs, file in os.walk(out_consensus_dir):
    for f in file:
        if f.endswith(".fa"):
            new = os.path.join(root, f)
            with open(new, "r+") as f:
                for line in f:
                    if not line.startswith(">"):
                        new_line = line.rstrip()
                        #print(len(new_line))
                        content_n = "N" * 29903
                        #print(len(content_n))
                        content_na = "N" * 29868 + "A" * 35
                        #print(len(content_n_a))
                        content_na_dif = "N" * 29870 + "A" * 33
                        content_nan = "N" * 29870 + "A" * 30 + "N" * 3
                        if new_line == content_n:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
                        elif new_line == content_na:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
                        elif new_line == content_nan:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
                        elif new_line == content_na_dif:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
"""
"""
for root, dirs, file in os.walk(out_consensus_dir):
    for f in file:
        if f.endswith(".fa"):
            new = os.path.join(root, f)
            with open(new, "r+") as f:
                for line in f:
                    if not line.startswith(">"):
                        new_line = line.rstrip()
                        #print(len(new_line))
                        content_n = "N" * 29866
                        #print(len(content_n))
                        if content_n in new_line:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)
"""
for root, dirs, file in os.walk(out_consensus_dir):
    for f in file:
        if f.endswith(".fa"):
            new = os.path.join(root, f)
            with open(new, "r+") as f:
                for line in f:
                    if not line.startswith(">"):
                        new_line = line.rstrip()
                        count = 0
                        letter_count = 0
                        #print(len(new_line))
                        for char in new_line:
                            if char.isalpha():
                                count += 1
                        for e in new_line:
                            if e == "N":
                                letter_count += 1
                        p = float(letter_count) / float(count) * 100
                        #print(len(content_n))
                        if p > 50:
                            name = new.split(".")[0]
                            print(name) # to see which is the name; it does work in juyter but not here
                            name_two = name.split("/")[-1]
                            print(name_two + " has not enough coverage and will be removed from the analysis")
                            f.truncate(0)

# cat all the consensus files in just one and quit the consensus word from the fasta name in order not to haver problemas eventually


pre_converge = os.path.join(out_consensus_dir, "*.fa")

converge = args.run_id + ".fa"
converge_file = os.path.join(out_consensus_dir, converge)

cmd_join = ["cat " + pre_converge + " | sed 's/Consensus_//g' > " + converge_file]
subprocess.run(cmd_join, shell=True)

        # run clustalo to align the multifasta file (if not, iqtree will raise an error, even if you choose fasta format)

check_create_dir(out_phyl_dir)

converged = args.run_id + ".fa"
converged_file = os.path.join(out_consensus_dir, converged)

phy = args.run_id + ".phy"
phy_file = os.path.join(out_phyl_dir, phy)

cmd_clustal = ["clustalo -i " + converged_file + " -o " + phy_file + " --outfmt=phy --threads 8"]
subprocess.run(cmd_clustal, shell=True)

        # run iqtree to obtain tree file. with this file treefile we upload to itol

phy_in = args.run_id + ".phy"
phy_file_in = os.path.join(out_phyl_dir, phy)

cmd_iq = ["iqtree -s " + phy_file_in + " -nt AUTO"]
subprocess.run(cmd_iq, shell=True)

## PANGOLIN LINEAGE

"""
pangolin_on = ["conda activate pangolin"]
subprocess.run(pangolin_on, shell=True)

converged = args.run_id + ".fa"
converged_file = os.path.join(out_consensus_dir, converged)

pangolin_file = args.run_id + ".lineage.report.csv"
final_pangolin_file = os.path.join(out_line_dir, pangolin_file)

run_pangolin = ["pangolin " + converged_file + " -o " + final_pangolin_file]
subprocess.run(run_pangolin, shell=True)

pangolin_off = ["conda deactivate"]
subprocess.run(pangolin_off, shell=True)
"""

# As it is very difficult to activate a conda environment from a python script, we will create and independent 
# BASH script to activate and run the pangolin program, and we call it from here.
# If everything goes ok, we will have a LIneage 

cmd_pangolin = ["bash -i pangolin.sh"]
subprocess.run(cmd_pangolin, shell=True)

