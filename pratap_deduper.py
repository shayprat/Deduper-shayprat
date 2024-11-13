#!/usr/bin/env python

import  argparse
import re

def get_args(): 
    parser = argparse.ArgumentParser(description="Python 3.12 program for Reference-based PCR duplicate removal of single-end reads", add_help=True)
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to deduplicated sam file", type=str, required=True)
    parser.add_argument("-u", "--umi_txt", help="designates file containing the list of UMIs", type=str, required=True)
    #parser.add_argument("-h", "--help", help="./pratap_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>", type=str)
    return parser.parse_args()

args = get_args()

input_sam = args.file 
out_sam = args.outfile
umi_txt = args.umi_txt

#initialize a set of valid/known UMIs

valid_umi_set = set()
with open(umi_txt, "r") as umis:
    for line in umis:
        umi = line.strip("\n")
        valid_umi_set.add(umi)

#Functions

def get_orientation(bitwise_flag: int) -> str:
    '''This function takes in a bitwise flag (second column in a SAM file) and tertermins whether the sequence is on the plus or minus strand'''
    if bitwise_flag & 16 == 16:
        return "minus"
    else:
        return "plus"
    

def get_SAM_info(line:str):
    field = line.split("\t")
    umi = field[0].split(":")[-1].strip()
    bitwise_flag = int(field[1].strip())
    chr = field[2].strip()
    cigar = field[5].strip()
    #pos is the 1-based leftmost mapping, NOT the 5prime start
    pos = int(field[3].strip())
    strand = get_orientation(bitwise_flag)
    

    return umi, chr, cigar, pos, strand


def get_corrected_position(cigar:str, strand:str, pos:int):
    broken_cigar = re.findall(r'(\d+)([A-Z]{1})', cigar)
    #print(broken_cigar)
    if strand == "plus":
        if broken_cigar[0][1] == "S":    
            pos -= int(broken_cigar[0][0])
    else:
        # if soft clipping on left side, skip
        if broken_cigar[0][1] == "S":
            for i in range(1,len(broken_cigar)):
                if broken_cigar[i][1] == "M":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "N":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "D":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "S":
                    pos += int(broken_cigar[-1][0])
        # If no soft clipping on left side
        else:
            for i in range(len(broken_cigar)):
                if broken_cigar[i][1] == "M":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "N":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "D":
                    pos += int(broken_cigar[i][0])
                if broken_cigar[i][1] == "S":
                    pos += int(broken_cigar[-1][0])
        pos -= 1
    return pos       

#I think these are all the functions I need .... lets write the big one


#initiate counters
header_lines = 0
invalid_umi = 0
pcr_dupes = 0
biological_dupes = 0
unique_reads = set()
#make sure SAM file is sorted before running this!
current_chr = str("1")
chrom_dict = {}
chrom_counter = 0 


with open (input_sam, "r") as input, open(out_sam, "w") as output:
    for line in input:
        if line == "":  
            #end of file
            break
        if line.startswith("@"): #header line
            header_lines += 1
            output.write(line)
        else: 
            #it's a read
            umi, chr, cigar, pos, strand = get_SAM_info(line)

            if umi not in valid_umi_set:
                invalid_umi += 1
                continue

            #print(f"{cigar}, {strand}, {pos}")

            corrected_pos = get_corrected_position(cigar, strand, pos)

            line_info = (chr, corrected_pos, strand, umi)

            if chr != current_chr: #we're on a new chromosome, reset variables and write the line to output
                chrom_dict[current_chr] = chrom_counter 
                current_chr = chr
                chrom_counter = 0
                unique_reads.clear()
                unique_reads.add(line_info)
                biological_dupes += 1
                chrom_counter += 1
                output.write(line)

             #we're on the same chromosome so we need to check if its a PCR duplicate
            else:
                if line_info not in unique_reads:
                    unique_reads.add(line_info)
                    biological_dupes += 1
                    chrom_counter += 1
                    output.write(line)

                else:
                    pcr_dupes += 1
            
for chrom in chrom_dict:
    print(f"{chrom}\t{chrom_dict[chrom]}")

print(f"Number of Header Lines = {header_lines}")
print(f"Number of PCR Duplicates = {pcr_dupes}")
print(f"Number of Biological Duplicates = {biological_dupes}")  
print(f"Number of Invalid UMIs = {invalid_umi}")  






        



            

    

        

        


