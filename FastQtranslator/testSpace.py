from fastQreader import FastQReader

file = FastQReader('illumina1.8.fastq')
#header, seq, qD, qScore
for header, seq, qD, qScore  in file.readFastQ():
    print('@'+header)
    print(seq)
    print('+'+qD)
    print(qScore)
    # print('header:', header)
    # print('seq:', seq)
    # print('qDes:',qD)
    # print('qScore:', qScore)
    # print()
