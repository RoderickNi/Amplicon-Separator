import os

def overlap_splicing(f,r,o,cpu=3):
    '''
    利用pandaseq拼接Pair-end reads
    '''
    os.system(f'pandaseq -f {f} -r {r} -w {o} -T {cpu} > {o}.log')
    

if __name__ == "__main__":
    
    input_fq_1=<Path>
    input_fq_2=<Path>
    output_fa=<Path>
    
    overlap_splicing(input_fq_1,input_fq_2,output_fa)