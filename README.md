# BlastImpute
Alignment based imputation

Requirements: 
```
Blast > 2.7.1+
Biopython
```

## prepare reference database

```bash
cd db
makeblastdb -in proteobacteria_rep.fasta -dbtype nucl 
```

## usuage

```bash
>python impute.py --input test.fasta --window 300
pos: 619 N > A (alignment evalue: 3.73174e-78)
pos: 677 N > G (alignment evalue: 2.86447e-84)

```

## method

- for every N it will create an alignment to reference database with the windowsize/2 on each direction
- alignments will get screened for the _opposite_ nucleotide. e.g. in this example the N (middle of the sequence) will be replaced by `A`.

```
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
```

