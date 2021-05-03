#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''Doc String'''

import sequenceAnalysis
import sys
from importlib import reload
reload(sequenceAnalysis)



########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    ''' 
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (0,100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def fileFormatter(startCodons, stopCodons, minGene, biggestGeneOnly, inFile='', outFile=''):
    genome = sequenceAnalysis.FastAreader(inFile)
    if outFile=='':
        pass
    else:
        file = open(outFile, 'w')
        sys.stdout = file

    for header, seq in genome.readFasta():
        print(header)
        orfsFound = sequenceAnalysis.OrfFinder(seq, startCodons, stopCodons)
        for currentOrf in sorted(orfsFound.finalORFlist(minGene, biggestGeneOnly), key=lambda a:(a[2],-a[0]), reverse=True):
            print(f'{currentOrf[3]:2s} {currentOrf[0]:>5d}..{currentOrf[1]:>5d} {currentOrf[2]:2d}') #{seq[currentOrf[0]-1:currentOrf[1]]}
        print('\n')

    if outFile =='': pass
    else: file.close()

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
def main(inFile = '', options = None):
    '''
    Find some genes.  
    '''
    thisCommandLine = CommandLine(options)
    ###### replace the code between comments.
    #print (thisCommandLine.args)
    # thisCommandLine.args.longestGene is True if only the longest Gene is desired
    # thisCommandLine.args.start is a list of start codons
    # thisCommandLine.args.stop is a list of stop codons
    # thisCommandLine.args.minGene is the minimum Gene length to include
    #
    #######
    fileFormatter(thisCommandLine.args.start, thisCommandLine.args.stop, thisCommandLine.args.minGene, thisCommandLine.args.longestGene, inFile, '')
    
if __name__ == "__main__":
    main() # delete this stuff if running from commandline