#!/user/bin/env python
#Created by Megan Neveau
#Last edited 6/9/17
import sys
import os.path
import argparse
import tempfile
import time
from subprocess import Popen, PIPE
from collections import defaultdict
import FisherExact

#check to ensure files have corresponding index files
def check_for_index(file):
    if os.path.isfile(file):
        print("Index file exists")
    else:
        print(file)
        print("Index file does not exist")
        sys.exit()

#execute bam-readcount
def bam_readcount(bam_file):
    bam_readcount_cmd = ['/opt/bam-readcount', '-f', ref_fasta, '-l', snp_file]
    print(ref_fasta, snp_file)
    bam_file_name = bam_file.split('/')
    output_file = os.path.join(tempdir, bam_file_name[-1] + '_rc.tsv')
    print(output_file)
    bam_readcount_cmd.append(bam_file)
    print(bam_readcount_cmd)
    execution = Popen(bam_readcount_cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = execution.communicate()
    if execution.returncode == 0:
        with open(output_file, 'wb') as output_fh:
            output_fh.write(stdout)
            print("Output written")
        output_fh.close()
        #quits if output file is empty
        if os.stat(output_file).st_size == 0:
            print("No output from bam-readcount")
            sys.exit()
    else:
        print("no output written")
        sys.exit(stderr)
    return output_file 

def parse_rc(rc_output_file, i):
    parse_file = os.path.join(tempdir, "parse_file_" + str(i))
    print(parse_file)
    with open(parse_file, 'w') as parse_file_fh:
        with open(rc_output_file, 'r') as rc_output:
            for line in rc_output:
                #replace : with /t
                edited_line = line.replace(":","\t")
                #split at each tab
                field = edited_line.split("\t")
                #write the array elements you want to keep onto parse_file 
                #(chr pos refbase A # C # G # T #)
                parse_file_fh.write(field[0] + "\t" + field[1] + "\t" + field [2] + "\t" + field[18] + "\t" + field[19] + "\t" + field[32] + "\t" + field[33] + "\t" + field[46] + "\t" + field[47] + "\t" + field[60] + "\t" + field[61] + "\t" + field[74] + "\t" + field[75] + " \n")
        rc_output.close()
        parse_file_fh.close()
    return parse_file

#make genotype calls with fisher exact test
def geno_calc(parsed_file, i):
    deep_enough = 0
    output_file = os.path.join(tempdir, "geno_calc" + str(i)) 
    with open(parsed_file, 'r') as data_file:
        with open(output_file, 'w') as out_f:
            for line in data_file:
                counter, letter, alpha, num = 0, 3, [None, None, None, None], [0,0,0,0]
                field = line.rstrip('\n').split("\t")
                depth = int(field[4]) + int(field[6]) + int(field[8]) + int(field[10])
                if depth >= 10:
                    deep_enough += 1

                    for val in field[4:11:2]:
                        if int(val) > 5:
                            alpha[counter] = field[letter]
                            num[counter] = val
                            counter += 1
                        letter += 2
                    if int(num[0]) >= 12 and int(num[1]) == 0:
                        field.append(alpha[0] + "\n")
                    elif counter == 2 and (int(num[0]) + int(num[1]) >= 12):
                                total = int(num[0]) + int(num[1])
                                half = total / 2
                                p_hetero = FisherExact.fisher_exact([[num[0], num[1]], [half, half]])
                                p_homo_1 = FisherExact.fisher_exact([[num[0], num[1]], [total, 0]])
                                p_homo_2 = FisherExact.fisher_exact([[num[0], num[1]], [0, total]])
                                if p_hetero > p_homo_1 and p_hetero > p_homo_2:
                                    field.append(alpha[0] + "/" + alpha[1] + " \n")
                                elif p_homo_1 > p_homo_2 and p_homo_1 > p_hetero:
                                    field.append(alpha[0] + "\n")
                                else:
                                    field.append(alpha[1] + "\n")   
                    else:
                        field.append("NA\n")
                    line = "\t".join(field)
                    out_f.write(line)
                else:
                    field.append("NA\n")
                    line = "\t".join(field)
                    out_f.write(line)
            #error out if no data with 10x coverage
            if deep_enough == 0:
                print("No site has greater than 10x coverage")
                sys.exit()
        out_f.close()
    data_file.close()
    return output_file

#open the output file from R and get the genotype
def get_genotypes(r_output_file,  n):
    total, snp = [], []
    valid_genotypes = 0
    with open (r_output_file, 'r') as r_file:
        for line in r_file:
            field = line.split('\t')
            if field[13] != "NA":
                total.append(field[0] + '_' + field[1])
                snp.append(field[0] + '_' + field[1] + '_' + field[13])
                valid_genotypes = valid_genotypes + 1
            #if full genotype output is wanted, creates dictionary to do so
            if output_geno != None:
                genotypes_dict[field[0] + "\t" + field[1] + "\t" + field[2]]["samp" + str(n)] = field[13]
        #quit if no valid genotypes
        if valid_genotypes <= 1:
            print("Not enough information to continue")
            sys.exit()
    r_file.close()
    return total, snp

#Program start:
start_time = time.time()
#makes sure that all of the required arguments are present
parser = argparse.ArgumentParser()
parser.add_argument("bam_file_1")
parser.add_argument("bam_file_2")
parser.add_argument("ref_fasta")
parser.add_argument("snp_file")
parser.add_argument("output_file")
parser.add_argument("--output_geno")
args = parser.parse_args()

bam_1 = args.bam_file_1 
bam_2 = args.bam_file_2
ref_fasta = args.ref_fasta
snp_file = args.snp_file
output = args.output_file
output_geno = args.output_geno

#checks for index for bam & reference FASTA files
check_for_index(bam_1 + ".bai") #will need + ".bai" & fasta will need + ".fa"
check_for_index(bam_2 + ".bai")
check_for_index(ref_fasta + ".fai")

#create temp directory for munging
tempdir = tempfile.mkdtemp()
print(tempdir)

#execute bam-readcount on both bam files, returns rc .tsv files
rc_1 = bam_readcount(bam_1)
rc_2 = bam_readcount(bam_2)

#parse the rc return file for the data fields that we want
parsed_1 = parse_rc(rc_1, 1)
parsed_2 = parse_rc(rc_2, 2)

#run R on the fields to determine genotypes
geno_out_1 = geno_calc(parsed_1, 1)
geno_out_2 = geno_calc(parsed_2, 2)

#get genotype data from R results
genotypes_dict = defaultdict(dict)
total_1, snp_1 = get_genotypes(geno_out_1, 1)
total_2, snp_2 = get_genotypes(geno_out_2, 2)

#write full genotype file in .bed format if specified in args
if output_geno != None:
    with open(output_geno, 'w') as geno_file:
        #write header
        geno_file.write("Chr" + "\t" + "Start" + "\t" + "Stop" + "\t" + "Genotype" + "\t" + "Sample1" + "\t" + "Sample2" + "\n")
        for key in genotypes_dict:
            field = key.split('\t')
            #makes sure each has 2 samples
            if 'samp1' in genotypes_dict[key] and 'samp2' in genotypes_dict[key]:
                geno_file.write(field[0] + "\t" +  field[1] + "\t"  + field[1] + "\t" + genotypes_dict[key]['samp1'] + "\t" + genotypes_dict[key]['samp2'] + "\n")
    geno_file.close()

#calculates matches
samples_wo_comparison = 0
all_samples = total_1 + total_2
total_samples = len(all_samples)
for samp in all_samples:
    if all_samples.count(samp) != 2:
        samples_wo_comparison += 1
total_samples_fixed = total_samples - samples_wo_comparison
total_matches_possible = total_samples_fixed / 2

total_uniq = len(list(set(snp_1 + snp_2)))

total_matches = total_samples_fixed - total_uniq

#write calculations to output file
with open(output, 'w') as output_file:
    output_file.write("Concordance Summary: \n")
    output_file.write("%7s\t%22s\t%15s\n"%("Matches", "Total Matches Possible","Percent Matched"))
    output_file.write("%7s\t%22s\t%15s\n"%(total_matches, total_matches_possible, (float(total_matches) / float(total_matches_possible)) * 100))

    print("%7s\t%22s\t%15s"%("Matches", "Total Matches Possible","Percent Matched"))
    print("%7d\t%22d\t%15f"%(total_matches, total_matches_possible, float(total_matches) / float(total_matches_possible) * 100))
output_file.close()

print("--- %s seconds ---"%(time.time() - start_time))
