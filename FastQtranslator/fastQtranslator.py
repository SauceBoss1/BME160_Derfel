#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)
'''
Translate one fastQ format to another fastQ format

written by Derfel Terciano & Sagarika Kannoly

Summary: 
    -Stream in a fastQ file, specify the type of fastQ format, and return the same fastQ file but in a different fastQ format that the user specifies

Notes:
    -use the -h flag in order to see all input and output arguments
    -program can only currently handle outputting to PHRED 64 and PHRED33
    -program run time will be displayed in sys.stderr
'''
from conversions import PhredConversions #remember to also have the solexaDictionary.py file in the same path
from commandLine import CommandLine
from fastQreader import FastQReader
import sys
import time
start_time = time.time() #starts the program runtime timer

def main(inFile = ''):
    '''Execute main translation sequence'''
    thisCL = CommandLine() #instiate/get the command arguments the user can use
    fastQfile = FastQReader(inFile) #instantiate the FastQReader

    if thisCL.args.PHRED64output: #figure out if the user has specified the p64 output, else fallback to p33 default
        outputConversion = True
    else:
        outputConversion = False

    for header, seq, qDesc, qScore in fastQfile.readFastQ(): #iterate through each fastQ read
        qConv = PhredConversions(qScore, seq, clOut= outputConversion) #for each read, initiate a phred conversion object
        newSeq = seq
        if thisCL.args.PHRED33input: #when P33in is specified, call the function that handles p33 conversions
            newQscore = qConv.p33toPhred()
        elif thisCL.args.PHRED64input: #when P64in is specified, call the function that handles p64 conversions
            newQscore = qConv.p64toPhred()
        elif thisCL.args.PHRED64Binput: #when P64Bin is specified, call the function that handles p64 B offset conversions
            newQscore, newSeq = qConv.p64bToPhred() #remember that for p64B, we need to change the sequence when we find 'B'
        elif thisCL.args.PHRED64SOLinput: #when P64SOLin is specified, call the function that handles p64 solexa conversions
            newQscore = qConv.pSolToPhred()

        #the lines below print the formatted output of the file conversion
        print('@'+header)
        print(newSeq)
        print('+'+qDesc)
        print(newQscore)

    print('run time (in seconds) : ', time.time()-start_time, file= sys.stderr) #print the final execution time of the program

if __name__ == '__main__':
    main()