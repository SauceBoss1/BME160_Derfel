'''
DON'T TURN IN THIS FILE!!! THIS IS ONLY A TESTING WORKSPACE
'''
from fastQreader import FastQReader
from solexaDictionary import solexaToPhred
import conversions
file = FastQReader('illumina1.3.fastq')
#file = FastQReader('illumina1.8.fastq')    
#file = FastQReader('Galaxy1.solexa.fastq')
#header, seq, qD, qScore
#####################################################################
# qScore = [-5, -4, -3, -2]

# newStr = ''
# for q in qScore:
#     newStr += newStr.join(chr(solexaToPhred[q]+64))
# print(newStr)
#####################################################################
# clInOption = 'p64bIn'
# clOutOp  = True
# for header, seq, qD, qScore  in file.readFastQ():
#     print('@'+header)
#     curSeq = seq
#     newQscore = qScore
#     converter = conversions.PhredConversions(qScore, seq, clOutOp)
#     if clInOption == 'p33in':
#         converter.p33toPhred()
#     elif clInOption == 'p64in':
#         converter.p64bToPhred()
#     elif clInOption == 'p64bIn':
#         newQscore, curSeq = converter.p64bToPhred()
#     elif clInOption == 'p64SolIn':
#         newQscore = converter.pSolToPhred()
#     print(curSeq)
#     print('+'+qD)
#     print(newQscore)

    # print('header:', header)
    # print('seq:', seq)
    # print('qDes:',qD)
    # print('qScore:', qScore)
    # print()
#####################################################################
#print(solexaToPhred)
#####################################################################

# import argparse

# parser = argparse.ArgumentParser(description = 'handles the fastQ command input',
#                                     add_help= True,
#                                     prefix_chars= '-',
#                                     usage = '%(prog)s [options] -option1[default] <input >output'
#                                     )

# parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False, help = 'PHRED33 input')
# parser.add_argument('-P64in','--PHRED64input', action = 'store', nargs= '?', const=True, default=False, help = 'PHRED64 input')
# parser.add_argument('-P64Bin','--PHRED64Binput', action= 'store', nargs='?',const=True, default=False, help= 'PHRED64 with B offset in quality values')
# parser.add_argument('-P64SOLin', '--PHRED64SOLinput', action='store',nargs='?',const='True',default=False, help = 'PHRED64 with SOLEXA interpretation of Q score')

# parser.add_argument('-P33out', '--PHRED33output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED33')
# parser.add_argument('-P64out', '--PHRED64output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED64')

# args = parser.parse_args()
#####################################################################

from commandLine import CommandLine
thisCL = CommandLine()

if thisCL.args.PHRED33input:
    print('phed33in')
elif thisCL.args.PHRED64input:
    print('phred64In')
elif thisCL.args.PHRED64Binput:
    print('phred64binput')
elif thisCL.args.PHRED64SOLinput:
    print('phred64SOLinput')

if thisCL.args.PHRED64output:
    print('p64 out')
else:
    print('p33 out')