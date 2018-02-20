#!/usr/bin/python
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

INPUT_FILE='test.fasta'
LEFT_WINDOW = 150
RIGHT_WINDOW = 150
E_VALUE_THRESH = 0.1

def replaceN (str):
	"This functions performs a blast search of str to replace a nt at given position pos"
	with open("query.fasta", "w") as fasta_query_file:
		fasta_query_file.write(">query\n")
		fasta_query_file.write(str)
	blastn_cline = NcbiblastnCommandline(query='query.fasta', db="db/proteobacteria_rep.fasta", evalue=0.001, outfmt=5, out="query.xml")
	stdout, stderr = blastn_cline()
	print(stderr)
	result_handle = open("test.xml")
	blast_records = NCBIXML.parse(result_handle)
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					ambigous_pos_alignment = hsp.query.find('N')
					if hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1] != "-":
						prediction = hsp.sbjct[ambigous_pos_alignment:ambigous_pos_alignment+1]
					else:
						prediction = "X" # there is no match
						hsp.expect = -1
	return(prediction, hsp.expect)

fasta_sequences = SeqIO.parse(open(INPUT_FILE),'fasta')

# iterate ofer fasta sequences and search for N
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	# loop over every N remaining in sequence and either replace to nucleotide or to X if there is no match
	while(sequence.find('N') > -1):
		ambigous_pos = sequence.find('N') # find the next N in the sequence, will be a -1 if no found
		print("processing fasta header" + name)
		print("ambigous sequence found at: " + str(ambigous_pos))
		print("extracting neighborhood around position " + str(ambigous_pos))
		subsequence = sequence[ambigous_pos - LEFT_WINDOW:ambigous_pos + RIGHT_WINDOW]
		prediction, evalue = replaceN(subsequence)
		print('pos: ' + str(ambigous_pos) + ' N > ' + prediction + ' (alignment evalue: ' + str(evalue) + ')')
		# correct sequence
		sequence_list = list(sequence)
		sequence_list[ambigous_pos] = prediction
		sequence = "".join(sequence_list)
