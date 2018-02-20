#!/usr/bin/python
import warnings
from Bio import BiopythonDeprecationWarning
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
import argparse


def replaceN(str, evalue, quite = True, ):
	"This functions performs a blast search of str to replace a nt at given position pos"
	# export substring as fasta file
	with open("query.fasta", "w") as fasta_query_file:
		fasta_query_file.write(">query\n")
		fasta_query_file.write(str)
	# blast query file to database
	blastn_cline = NcbiblastnCommandline(query='query.fasta', db="db/proteobacteria_rep.fasta", evalue=0.001, outfmt=5, out="query.xml")
	stdout, stderr = blastn_cline()
	result_handle = open("query.xml")
	blast_records = NCBIXML.parse(result_handle)
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					ambigous_pos_alignment = hsp.query.find('N')
					if hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1] != "-" and hsp.query[ambigous_pos_alignment:ambigous_pos_alignment+1] == 'N':
						prediction = hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1]
						query_char = hsp.query[ambigous_pos_alignment:ambigous_pos_alignment+1]
						if not quite:
							print(hsp.query[100:200])
							print(hsp.match[100:200])
							print(hsp.sbjct[100:200])
					else:
						prediction = "X" # there is no match
						query_char = "N"
						hsp.expect = -1
	return(query_char, prediction, hsp.expect)


parser = argparse.ArgumentParser(description='Impute by blast')
parser.add_argument('--input', help='path to input fasta file')
parser.add_argument('--window', type=int, help='window size', default=300)
parser.add_argument('--evalue', type=float, help='evalue threshold for blast', default=0.01)

args = parser.parse_args()

fasta_sequences = SeqIO.parse(open(args.input),'fasta')

# iterate ofer fasta sequences and search for N
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	# loop over every N remaining in sequence and either replace to nucleotide or to X if there is no match
	while(sequence.find('N') > -1):
		ambigous_pos = sequence.find('N') # find the next N in the sequence, will be a -1 if no found
		#print("processing fasta header" + name)
		#print("ambigous sequence found at: " + str(ambigous_pos))
		#print("extracting neighborhood around position " + str(ambigous_pos))

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
		query, prediction, evalue = replaceN(subsequence, args.evalue,  quite = True)
		print('pos: ' + str(ambigous_pos) + ' ' + query + ' > ' + prediction + ' (alignment evalue: ' + str(evalue) + ')')
		# correct sequence
		sequence_list = list(sequence)
		sequence_list[ambigous_pos] = prediction
		sequence = "".join(sequence_list)
