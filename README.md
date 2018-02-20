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
>python impute.py --input test.fasta --window 300 --verbose
pos: 619 N > A (mean alignment evalue: 3.73174e-79, 10 alignments)
pos: 677 N > G (mean alignment evalue: 2.86447e-85, 9 alignments)
pos: 941 N > G (mean alignment evalue: 4.6918e-22, 4 alignments)

```

without `--verbose` it will print parts of the alignment, you can also use `--evalue` parameter e.g. `--evalue 0.00001` to increase the threshold of the alignments that will be used for imputation

```bash
>python impute.py --input test.fasta --window 300
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| || ||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCCAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
TGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATTT
| ||||||||||||||||||||||||||||| ||||| ||||||||||| ||||| |||||| |||| |||||||||||||||||||||||||| |||||
TAGTAGTAGGACTCGTTATTAATCCAAATAGATTAATAGAGATAAGAGAAGCTAGGTTAAATTTATTGCAAATTAACGAAAATAAAAGCTATACAGATTT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
||||| ||||||||||| |||||||||||||||||||||||||||||||| |||||||||||| |||| |||||||| ||||||||||||||||||||||
TTGGTGGTAGGACTCGTCATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATTTATTGCAAATTAATGAAAATAAAAGCTATACGGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|| |||||||| ||||| |||||||||||||| ||||| ||||||||||| ||||| |||||| |||| |||||||||||||||||||||||||| ||||
TTAGTAGTAGGGCTCGTCATTAATCCAAATAGATTAATAGAGATAAGAGAAGCTAGGTTAAATTTATTGCAAATTAACGAAAATAAAAGCTATACAGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|| |||||||| ||||| |||||||||||||| ||||| ||||||||||| ||||| |||||| |||| ||||||||||| |||||||||||||| ||||
TTAGTAGTAGGGCTCGTAATTAATCCAAATAGATTAATAGAGATAAGAGAAGCTAGGTTAAATTTATTGCAAATTAACGACAATAAAAGCTATACAGATT
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|| ||||||||||| ||||||||||||||||||||||||||||||||||| |||||||||||| |||| |||||||||||||||||||||||||| ||||
TTAGTAGTAGGACTTGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATTTATTGCAAATTAACGAAAATAAAAGCTATACAGATT
TGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATTT
| |||||||| || || |||||||||||||| ||||| || ||||||||  ||||  | ||| |  |||||||||| ||||||| ||  ||||| |||||
TAGTAGTAGGGCTTGTAATTAATCCAAATAGATTAATAGAAATAAGAGAAACTAGGCTGAATTTGCTACAAATTAATGAAAATAGAAATTATACAGATTT
pos: 619 N > A (mean alignment evalue: 3.73174e-79, 10 alignments)
...
```


## method

- for every N it will create an alignment to reference database with the windowsize/2 on each direction
- alignments will get screened for the _opposite_ nucleotide. e.g. in this example the N (middle of the sequence) will be replaced by `A`.
- if there are multiple alignments the majority will be used for imputation
- the evalue is the mean evalue of all alignments

```
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGANGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
|||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
TTGGTAGTAGGACTCGTTATTAATCCAAATAGGTTAATTGAGATAAGAGAAGCTAGATTAAATCTATTACAAATTAACGAAAATAAAAGCTATACGGATT
```

