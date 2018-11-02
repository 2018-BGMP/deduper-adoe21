import argparse
import re
import numpy as np

def get_arguments():
	parser = argparse.ArgumentParser(description="Program to remove PCR duplicates", add_help="False")
	parser.add_argument("-s", "--file", help="SAM file to be checked for PCR duplicates", required=True, type=str)
	parser.add_argument("-u", "--umi", help="File containing known UMIs", required=False, type=str)
	parser.add_argument("-p", "--paired", help="Tell deduper if reads are paired end", required=False, action="store_true")
	return parser.parse_args()
	
args = get_arguments()

def isolate_UMI(string):
	'''Returns UMI from SAM qname'''
	UMI = re.match('.*?([A-Z]+)$', string).group(1)
	return UMI
	
def start_pos_check(start1,cigar1):
	'''Checks if start position of two SAM reads match'''
	firstletter1 = re.search('([A-Z])',cigar1).group(1)
	if firstletter1 == "S":
		firstnumber1 = re.search('(\d+)',cigar1).group(1)
		start1 = int(start1) - int(firstnumber1)
	return start1

def end_pos_check(start1,cigar1):
	'''Checks if the end position of two SAM reads match'''
	sum1 = 0
	if cigar1[-1] == "S":
		sum1 += int(re.search('(\d+)$',cigar1[:-1]).group(1))
	cigarlist1 = cigar1.split("M")
	cigarlist1 = cigarlist1[:-1]
	for elem in cigarlist1:
		if elem[-1].isnumeric() == True:
			sum1 += int(re.search('(\d+)$',elem).group(1))
	Dlist1 = cigar1.split("D")
	Dlist1 = Dlist1[:-1]
	for elem in Dlist1:
		if elem[-1].isnumeric() == True:
			sum1 += int(re.search('(\d+)$',elem).group(1))
	Nlist1 = cigar1.split("N")
	Nlist1 = Nlist1[:-1]
	for elem in Nlist1:
		if elem[-1].isnumeric() == True:
			sum1 += int(re.search('(\d+)$',elem).group(1))
	end_pos1 = int(start1) + int(sum1)
	return end_pos1
	
def check_chrom(string):
	'''Extracts chromosome from RNAME of SAM line'''
	chrom = int(string)
	return chrom
	
def strand_check(string):
	'''Extracts strandedness from bitwise flag of SAM line'''
	bitflag = string
	if ((int(bitflag) & 16) != 16):
		return "+"
	else:
		return "-"

	
chromdump = 0
UMI_set = set()
if args.umi != None:
	stop = "start"
	with open(args.umi, "r") as umi:
		while stop != "":
			item = umi.readline()
			stop = item
			if stop == "":
				continue
			UMI_set.add(item[:-1])

outfile = re.search('([a-zA-Z0-9]+\.sam)$',args.file).group(1)
			
with open(args.file,"r") as sam:
	with open(outfile + "_deduped","a") as out:
		for line in sam:
			breaker = False
			if line.startswith("@"):
				out.write(line)
				continue
			read = line.strip().split()
			array = np.array(read)
			qname = array[0]
			rname = array[2]
			flag = array[1]
			start = array[3]
			cigar = array[5]
			UMI = isolate_UMI(qname)
			chrom = check_chrom(rname)
			if chrom != chromdump:
				UMI_dict = {}
				chromdump = chrom
			strand = strand_check(flag)
			if strand == "+":
				position = start_pos_check(start,cigar)
			if strand == "-":
				position = end_pos_check(start,cigar)
			important = (chrom, strand, position)
			if len(UMI_set) < 1:
				if UMI in UMI_dict:
					for elem in UMI_dict[UMI]:
						if elem == important:
							breaker=True
							break
					if breaker:
						continue
					out.write(line)
					UMI_dict[UMI].add(important)
				if UMI not in UMI_dict:
					UMI_dict[UMI] = set()
					UMI_dict[UMI].add(important)
					out.write(line)
			if len(UMI_set) > 1:
				if UMI in UMI_dict and UMI in UMI_set:
					for elem in UMI_dict[UMI]:
						if elem == important:
							breaker=True
							break
					if breaker:
						continue
					out.write(line)
					UMI_dict[UMI].add(important)
				if UMI not in UMI_dict and UMI in UMI_set:
					UMI_dict[UMI] = set()
					UMI_dict[UMI].add(important)
					out.write(line)
				if UMI not in UMI_set:
					continue
				
