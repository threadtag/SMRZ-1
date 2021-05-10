## SMRZ-1
Scripts to find candidate SMRZ-1 ribozyme sequence in meta-genome

## Requirement
1. Perl
2. Ruby
3. Vienna RNA package

## usage exmample
```shell
perl locate_subseq.pl Bacteria.fna "\\w{8}T[AGT]CTACG\\w(\\w{1,})\\w[ACG]CGGGG\\w{8}" >bac_candidates.fasta
ruby smrz_candi.rb -t 37 -d bac_candidates.fasta >bac_found.fasta
```
