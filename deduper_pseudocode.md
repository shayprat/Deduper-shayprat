Write up a strategy for writing a Reference Based PCR Duplicate Removal tool
i.e. given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read)

**Defining the problem**
If a library is amplified prior to sequencing, we could be introducing PCR duplicates to the system. For DNA/WGS applications, duplicates aren't an issue - you just need aligments that map to the reference genome. 

However, for RNA/WES/WTS applications, duplicates hinder relative abundance measurements (i.e. are there more copies of a gene sequence due to that gene being actually transcribed more or is it because that sequences had bias in PCR due to GC/AT richness, difference in Tms, secondary structures etc)

To help identify PCR duplicates, the library prep (pre amplifciation step)introduces UMIs, unique molecular identifiers, which are short random barcodes that serve as an ID tag. 

PCR Duplicate criteria (all must be TRUE):
1. the UMIs are the same (look at QName for UMI - column 1 of SAM file)
2. present on the same strand. Need to check sequence strandedness i.e. plus or minus strand (bitwise flag - column 2 of SAM file)
    NOTE: if minus strand, need the reverse compliment of the UMI!!
3. Are located on the same chromosome (RName - column 3 of SAM file)
4. Same start position on the 5' end (Pos - column 4 of SAM file)
    NOTE: take into account minus strand and CIGAR string (softclipping!!)

**write examples**


**Pseudocode**

```
Split SAM by chromosome and split by UMI (or use samtoools to sort by chromosome) outside of .py script

In .py script:

import  argparse

def get_args(): 
    parser = argparse.ArgumentParser(description="Python 3.12 program for Reference-based PCR duplicate Removal of single-end reads")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to deduplicated sam file", type=str, required=True)
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs", type=str, required=True)
    parser.add_argument("-h", "--help", help="./pratap_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>", type=str, required=False)
    return parser.parse_args()

in = get_args().file 
out = get_args().outfile
umi = get_args().umi

#initialize set of known/valid UMIs
valid_umi = ()

open list of valid UMI file (i.e STL96.txt)
    while True:
        read line and add to valid_umi
    if line == "" 
        break


intialize empty set to keep track of where you are
current_info (type=set) = (UMI, Chrom, corrected start position)

open input SAM file for reading, and output SAM file for writing
   While True loop?
    if line == ""
        break

    if line starts with "@", it's a header 
        write to output file
        continue to next line

    else next line isn't a header
        new_umi = extract from line.split(\t),  QName (COL 1)
        new_chrom = extract from line.split(\t), RName (COL 3)
        strand = get_orientation(line)
        position = get_corrected_position(line)

        if strandedness = minus, generate rev_comp of umi (not sure if this is needed)
            new_umi = bioinfo.rev_comp(new_umi)

        else (it's a positive strand)
        continue

        if new_umi not in valid_umi
        discard? and continue to next line and restart process

        create a tuple of [new_umi, new_chrome, position]
        if tuple is in current_info set, it's a duplicate
            discard and continue
        else it's unique
            write to output file

        
        
        #not going this route, lets save everything to a set (current_info)
        #new_chrom != current_chrom and (or?) new_umi != current_umi
           # update/reset current_chrom to new_chrom
           # update/reset current_umi to new_umi
           # new_position = get_correct_position(line)
```
    

**High-level functions**

1. Strandedness
```
def get_orientation(line: str) ->: str
    ```determines whether a sequences is on the plus or minus strand depending on BITWISE FLAG```
    flag = line.split()[1]
    if flag is [whatever the bit value is for minus, 16?]
        strand = "minus"
    else
        strand = "plus"
return strand
```

**Input** = NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
**Output** = plus

**Input** = NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	16	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
**Output** = minus
            
2. Corrected Position
```
def get_corrected_position(line: str) ->: int
    ```determines left-most start postion taking softclipping into account```
    sam_position = line.split()[3]
    cigar_string = line.split()[5]
    corrected_position = sam_position adjusted by integer before the "S"
    Take into account INDELS????
return corrected_position

```
**input** = NS500451:154:HWKTMBGXX:1:11101:13631:1222:TCGTAGGT	0	2	76944432	36	71M	*	0	0	GCCAGTACCAGTAGATTCGGTCAGATCGCTCAATTTTTACACCATTTTTATACCATTCACACTCAGGGTCT	/AEEEEEEEEEEEEEEEEEEEEEEAEEE<EEEEEEEEEEE6<EEEEEEEEEEEEEEEEEEAEEEEE<AEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU

**output** = 76944432

**input** = NS500451:154:HWKTMBGXX:1:11101:13631:1222:TCGTAGGT	0	2	76944432	36	3S71M	*	0	0	GCCAGTACCAGTAGATTCGGTCAGATCGCTCAATTTTTACACCATTTTTATACCATTCACACTCAGGGTCT	/AEEEEEEEEEEEEEEEEEEEEEEAEEE<EEEEEEEEEEE6<EEEEEEEEEEEEEEEEEEAEEEEE<AEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU

**output** = 76944429










