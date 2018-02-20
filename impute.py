from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIStandalone
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

input_file='test.fasta'
left_window = 100
right_window = 100
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
E_VALUE_THRESH = 0.04

for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	ambigous_pos = sequence.find('N') # find the next N in the sequence, will be a -1 if no found
	if ambigous_pos > 0: # only proceed if there is a ambigous char in the sequence
		print("processing fasta file" + name)
		print("ambigous sequence found at: " + str(ambigous_pos))
		print("extracting neighborhood around position " + str(ambigous_pos))
		subsequence = sequence[ambigous_pos - left_window:ambigous_pos + right_window]
		#query_record = SeqRecord(Seq(subsequence, IUPAC.protein),id="query", name="query", description="query sequence")
		with open("query.fasta", "w") as fasta_query_file:
			fasta_query_file.write(">query\n")
			fasta_query_file.write(subsequence)
		blastn_cline = NcbiblastnCommandline(query='query.fasta', db="db/proteobacteria_rep.fasta", evalue=0.001, outfmt=5, out="query.xml")
		stdout, stderr = blastn_cline()
		print(stderr)
		result_handle = open("test.xml")
		blast_records = NCBIXML.parse(result_handle)

		for blast_record in blast_records:
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					if hsp.expect < E_VALUE_THRESH:
						print('****Alignment****')
						print('sequence:', alignment.title)
						print('length:', alignment.length)
						print('e value:', hsp.expect)
						print(hsp.query[0:75] + '...')
						print(hsp.match[0:75] + '...')
						print(hsp.sbjct[0:75] + '...')