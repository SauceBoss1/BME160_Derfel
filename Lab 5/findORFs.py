#!/usr/bin/env python3
# Name: Derfel Terciano (dtercian)
# Group Members: None
'''Doc String'''

#PSEUDO CODE FOR orfFinder METHOD IN sequenceAnalysis.orfFinder

#########################################################################################################
#lists that are a MUST: startCodons and stopCodons
#lists needed: startPosition and StopPosition
#
#for frames in frame 1, 2, 3 (this will be used to shift the search by 1 character to the right):
#   (remember) clear the stopPosition list
#   for codonPosition from range (frames to the length of the sequence)
#       codon = extract actual codon nucleotides using codon position       
#
#       if codon is a start codon:
#           append the current codon position
#
#
#       if codon is a stop codon:
#           append the current position to stopPosition
#           if startPos is not empty:
#               if length reqs have been met:
#                   save current ORF info
#               if we have shifted frames (i.e not frame 1) and the startPos[0] is the same as the frame:
#                   save current ORF info that starts at position 1
#               if multiple starts have been found (only if this option is enabled):
#                   save ORF info of all of those (from each start to current codonPos)
#               clear startPos list
#       if there no other ORFs in the frame (and) there are not starts (and) this is the only stop found:
#          save ORF info if the len requirements have been met (this is a dangling stop)
#           clear startPos list
#
#   if we're at the end, and there is still a start codon:
#       save ORF if length reqs are met from start pos to the end of the seq
#   if no ORFs have been found:
#       save entire length as an ORF since it could be a potential gene candidate
##########################################################################################################


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

    for header, seq in genome.readFasta():  #print a message if the minGene > seq length
        print(header)

        if len(seq) < minGene: #if the minGenelength is bigger than seq length, algorithm will not function properly
            print('The specified minimum gene length is bigger than the length of the sequence.\n\tChange the genelength for this program to work. \n')
        else:
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
