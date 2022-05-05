# Amplicon-Separator
A pipline for  amplicons from different population separation  (Suitable for illumina sequencing platform)

## step1:  
- Dependencies:   
<b>pandaseq 2.11</b>    

Using *Overlapper.py* for overlapping splicing    
Please input pair-end sequencing result fastq1,fastq2 and the path of <output.fa>
  
## step2:    
- Dependencies:   
<b>seqtk 1.3</b>  
<b>fastx_toolkit 0.0.14</b>  

Using *Separator.bash* for amplicons separating    
Please input <output.fa from step1> and <primers file>    
The formate of <primers file> was show in *Merge_DNA_primer.txt*
