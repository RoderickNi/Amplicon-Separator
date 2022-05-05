# Amplicon-Separator
A pipline for  amplicons from different population separation  (Suitable for illumina sequencing platform)

- step1:    
Using *Overlapper.py* for overlapping splicing    
Please input pair-end sequencing result <fq1> & <fq2> and the path of <output.fa>
  
- step2:
Using *Separator.bash* for amplicons separating    
Please input <output.fa from step1> and <primers file>    
The formate of <primers file> was show in *Merge_DNA_primer.txt*
