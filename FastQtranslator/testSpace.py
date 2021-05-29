from fastQreader import FastQReader
from solexaDictionary import solexaToPhred
#file = FastQReader('illumina1.3.fastq')
file = FastQReader('illumina1.8.fastq')    
#file = FastQReader('Galaxy1.solexa.fastq')
#header, seq, qD, qScore

qScore = [-5, -4, -3, -2]

newStr = ''
for q in qScore:
    newStr += newStr.join(chr(solexaToPhred[q]+64))
print(newStr)

# for header, seq, qD, qScore  in file.readFastQ():
#     print('@'+header)
#     print(seq)
#     print('+'+qD)
#     print(qScore)
    #print()
    # print('header:', header)
    # print('seq:', seq)
    # print('qDes:',qD)
    # print('qScore:', qScore)
    # print()

#print(solexaToPhred)

