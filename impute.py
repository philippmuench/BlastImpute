#!/usr/bin/python
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import argparse
import os

def replaceN(str, evalue, quite = True, ):
	"This functions performs a blast search of str to replace a nt at given position pos"
	# export substring as fasta file
	with open("query.fasta", "w") as fasta_query_file:
		fasta_query_file.write(">query\n")
		fasta_query_file.write(str)
	# blast query file to database
	blastn_cline = NcbiblastnCommandline(query='query.fasta', db="db/proteobacteria_rep.fasta", 
		evalue=evalue, outfmt=5, out="query.xml")
	stdout, stderr = blastn_cline()
	result_handle = open("query.xml")
	blast_records = NCBIXML.parse(result_handle)
	pred_list = [] # list to store all prediction
	eval_list = [] # list to store evalue for each alignment that will be used for imputation
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					ambigous_pos_alignment = hsp.query.find('N')
					if hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1] != "-" and hsp.query[ambigous_pos_alignment:ambigous_pos_alignment+1] == 'N':
						prediction = hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1]
						query_char = hsp.query[ambigous_pos_alignment:ambigous_pos_alignment+1]
						pred_list.append(prediction)
						eval_list.append(hsp.expect)
						if not quite:
							print(hsp.query[100:200])
							print(hsp.match[100:200])
							print(hsp.sbjct[100:200])
					else:
						prediction = "X" # there is no match
						query_char = "N"
						hsp.expect = 1
	return(query_char, pred_list, eval_list)

def find_majority(k):
	#This functions returns the majority character of a list and its count
    myMap = {}
    maximum = ( '', 0 ) # (occurring element, occurrences)
    for n in k:
        if n in myMap: myMap[n] += 1
        else: myMap[n] = 1
        # Keep track of maximum on the go
        if myMap[n] > maximum[1]: maximum = (n,myMap[n])
    return maximum

# parse arguments
parser = argparse.ArgumentParser(description='Impute by blast')
parser.add_argument('--input', help='path to input fasta file')
parser.add_argument('--output', help='where to write the corrected sequence', default='output.fasta')
parser.add_argument('--window', type=int, help='window size', default=300)
parser.add_argument('--evalue', type=float, help='evalue threshold for blast', default=0.01)
parser.add_argument("--verbose", help="increase output verbosity",
                    action="store_true")
args = parser.parse_args()

# load fasta file
fasta_sequences = SeqIO.parse(open(args.input),'fasta')

# check output file
try:
    os.remove(args.output)
except OSError:
    pass

# iterate ofer fasta sequences and search for N
with open(args.output, "w") as fasta_output_file:
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		print("processing " + name)
		# loop over every N remaining in sequence and either replace to nucleotide or to X if there is no match
		while(sequence.find('N') > -1):
			ambigous_pos = sequence.find('N') # find the next N in the sequence, will be a -1 if no found
			# check if window would exeed the sequence
			if ambigous_pos - args.window/2 < 0: 
				left_window = ambigous_pos # reduce the window site if its too near on the beginning of the sequence
			else:
				left_window = args.window/2
			if ambigous_pos + args.window/2 > len(sequence):
				right_window = len(sequence) # reduce right window if window would exeed sequence length
			else: 
				right_window = args.window/2 
			subsequence = sequence[ambigous_pos - left_window:ambigous_pos + right_window]
			query, pred_list, eval_list = replaceN(subsequence, args.evalue,  quite = args.verbose)
			if eval_list:
				majority = find_majority(pred_list) # get the most occuring prediction from all alignments
				print('pos: ' + str(ambigous_pos) + ' ' + query + ' > ' + majority[0] + ' (mean alignment evalue: ' 
				+ str(sum(eval_list) / float(len(eval_list))) + ', ' + str(majority[1]) + ' alignments)')
			# correct sequence
			sequence_list = list(sequence)
			sequence_list[ambigous_pos] = majority[0]
			sequence = "".join(sequence_list)
		# write imputed sequence

		fasta_output_file.write('>' + name + '\n')
		fasta_output_file.write(sequence.replace('X', 'N') + '\n') # X was used to decode ambigous sequences that cannot be corrected

