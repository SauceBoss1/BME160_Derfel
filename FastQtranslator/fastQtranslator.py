#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)

from conversions import PhredConversions
from commandLine import CommandLine
from fastQreader import FastQReader
import sys
import time
start_time = time.time()

def main(inFile = ''):
    thisCL = CommandLine()
    fastQfile = FastQReader(inFile)

    if thisCL.args.PHRED64output:
        outputConversion = True
    else:
        outputConversion = False

    for header, seq, qDesc, qScore in fastQfile.readFastQ():
        qConv = PhredConversions(qScore, seq, clOut= outputConversion)
        newSeq = seq
        if thisCL.args.PHRED33input:
            newQscore = qConv.p33toPhred()
        elif thisCL.args.PHRED64input:
            newQscore = qConv.p64toPhred()
        elif thisCL.args.PHRED64Binput:
            newQscore, newSeq = qConv.p64bToPhred()
        elif thisCL.args.PHRED64SOLinput:
            newQscore = qConv.pSolToPhred()
        print('@'+header)
        print(newSeq)
        print('+'+qDesc)
        print(newQscore)
    print('run time (in seconds) : ', time.time()-start_time, file= sys.stderr)

if __name__ == '__main__':
    main()