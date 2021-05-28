from fastQreader import FastQReader
from solexaDictionary import solexaToPhred
#file = FastQReader('illumina1.3.fastq')
file = FastQReader('illumina1.8.fastq')    
#file = FastQReader('Galaxy1.solexa.fastq')
#header, seq, qD, qScore
for header, seq, qD, qScore  in file.readFastQ():
    print('@'+header)
    print(seq)
    print('+'+qD)
    print(qScore)
    #print()
    # print('header:', header)
    # print('seq:', seq)
    # print('qDes:',qD)
    # print('qScore:', qScore)
    # print()

#print(solexaToPhred)