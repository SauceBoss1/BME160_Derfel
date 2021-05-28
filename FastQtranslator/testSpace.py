from fastQreader import FastQReader

file = FastQReader('Galaxy1.solexa.fastq')

for header, seq, qD, qScore in file.readFastQ():
    print('header:', header)
    print('seq:', seq)
    print('qDes:',qD)
    print('qScore:', qScore)
    print()
