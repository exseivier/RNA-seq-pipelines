#!/usr/bin/env python
import sys
import pickle

def load_genome(filename):
	"""(STR) -> HASH[STR:STR]
	Requires a fasta file name of the reference genome sequence
	It takes the chromosomes of a genome and their sequences
	and stores them into a HASH variable "HASH[header:sequence]"
	"""
	# Counting total headers
	IN = open(filename, "r")
	total = 0
	for line in IN:
		if line[0] == ">":
			total += 1
	IN.close()
	# Declaring and Assigning GENOME HASH variable
	GENOME = {}
	counter = 0
	IN = open(filename, "r")
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			counter += 1
			print "\rParsing %s [%f" % (line, (float(counter)/float(total))*100) + " %]"
			header = line
		else:
			try:
				GENOME[header] += line
			except Exception:
				GENOME[header] = line
	IN.close()
	# Dumping HASH variable to a file
	if len(GENOME) != 0:
		OUT = open("%s.dat" % filename, "wb")
		pickle.dump(GENOME, OUT)
		OUT.close()
	# Returning HASH
	return GENOME

def rev_comp(dna):
	"""(STR) -> STR
	Requires a string of characters
	Returns the reverse complement of a dna string
	"""
	rc = ""
	d_nt = {
			"A":"T",
			"C":"G",
			"G":"T",
			"T":"A",
			"N":"N"
			}
	for nt in dna:
		rc = d_nt[nt] + rc
	return rc

def main():
	"""
	MAIN PROGRAM
	"""
	# Trying with Genome fasta.dat file to load HASH data structure. Checking if file exists
	try:
		DF = open("%s.dat" % sys.argv[1], "rb")
		genome = pickle.load(DF)
		DF.close()
	# If an error is rised, then HASH structure is created with load_genome()
	except IOError:
		genome = load_genome(sys.argv[1])
	# Printing data structure to stdout
	for head, seq in genome.iteritems():
		print "%s\n%s\n" % (head, seq[:70])
	
	# Declaring variables to parse gff entries
	gstart = 0
	gend = 0
	cds_start = 0
	cds_end = 0
	name = ""
	locus = ""
	orient = ""
	CDS = []
	# Open output file to write 5' and 3' UTR's
	DOWN_UTR = open("%s.3utr.fasta" % sys.argv[1], "w+")
	UP_UTR = open("%s.5utr.fasta" % sys.argv[1], "w+")
	# Parsing gff file entries
	GFF_IN = open(sys.argv[2], "r")
	for entry in GFF_IN:
		entry = entry.strip()
		if entry[0] == "#":
			continue
		DATA = entry.split("\t")
		if DATA[2] == "gene":
			if len(CDS) != 0:
				print "Processing gene "+name+"..."
				cds_start = min(CDS)
				cds_end = max(CDS)
				"""
				header = "%s_%s_%d_%d" % (locus, name, cds_start, cds_end)
				print ">"+header
				print gstart, cds_start, cds_end, gend
				"""
				#"""
				if orient == "+":
					header = "%s_%s_gs%d_ge%d_cs%d_ce%d" % (locus, name, gstart, gend, cds_start, cds_end)
					if cds_start - gstart != 0:
						UP_UTR.write("%s\n%s\n" % (header, genome[">"+locus][gstart:cds_start]))
					if gend - cds_end != 0:
						DOWN_UTR.write("%s\n%s\n" % (header, genome[">"+locus][cds_end:gend]))
				elif orient == "-":
					header = "%s_%s_gs%d_ge%d_cs%d_ce%d" % (locus, name, gstart, gend, cds_start, cds_end)
					if cds_start - gstart != 0:
						DOWN_UTR.write("%s\n%s\n" % (header, rev_comp(genome[">"+locus][gstart:cds_start])))
					if gend - cds_end != 0:
						UP_UTR.write("%s\n%s\n" % (header, rev_comp(genome[">"+locus][cds_end:gend])))
				#"""
			gstart = int(DATA[3])
			gend = int(DATA[4])
			orient = DATA[6]
			name = DATA[8].split(";")[0].split("=")[1]
			locus = DATA[0]

		elif DATA[2] == "CDS":
			CDS = []
			CDS.append(int(DATA[3]))
			CDS.append(int(DATA[4]))

	DOWN_UTR.close()
	UP_UTR.close()

if __name__ == "__main__":
	main()

