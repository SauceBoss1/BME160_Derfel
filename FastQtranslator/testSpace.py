from fastQreader import FastQReader
from solexaDictionary import solexaToPhred
import conversions
file = FastQReader('illumina1.3.fastq')
#file = FastQReader('illumina1.8.fastq')    
#file = FastQReader('Galaxy1.solexa.fastq')
#header, seq, qD, qScore

# qScore = [-5, -4, -3, -2]

# newStr = ''
# for q in qScore:
#     newStr += newStr.join(chr(solexaToPhred[q]+64))
# print(newStr)
clInOption = 'p64bIn'
clOutOp  = 'p64'
for header, seq, qD, qScore  in file.readFastQ():
    print('@'+header)
    curSeq = seq
    newQscore = qScore
    converter = conversions.PhredConversions(qScore, seq, clOutOp)
    if clInOption == 'p33in':
        converter.p33toPhred()
    elif clInOption == 'p64in':
        converter.p64bToPhred()
    elif clInOption == 'p64bIn':
        newQscore, curSeq = converter.p64bToPhred()
    elif clInOption == 'p64SolIn':
        newQscore = converter.pSolToPhred()
    print(curSeq)
    print('+'+qD)
    print(newQscore)

    # print('header:', header)
    # print('seq:', seq)
    # print('qDes:',qD)
    # print('qScore:', qScore)
    # print()

#print(solexaToPhred)

