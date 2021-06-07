#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)
'''
Handle all command line options for the fastQ translator

Written by: Derfel Terciano & Sagarika Kannoly

Contains:
    -CommandLine() [object] - an object that handles command line options
'''
class CommandLine():
    '''
    Handle all command line arguments for fastQ translation

    Written By: Derfel Terciano & Sagarika Kannoly

    Options:
        Input handling:
            -args.PHRED33input
            -args.PHRED64input
            -args.PHRED64Binput
            -args.PHRED64SOLinput
        Output Handling:
            -args.PHRED64output
            -args.PHRED33output
    
    For usage use the -h flag in the command line
    '''

    def __init__(self, inOpts = None):
        '''
        Usage: 
            1) instantiate class
            2) use .args.argument as an input in your program
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'handles the fastQ command input',
                                                add_help= True,
                                                prefix_chars= '-',
                                                usage = '%(prog)s [options] -option1[default] <input >output'
                                                )

        self.parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False, help = 'PHRED33 input') #PHRED 33 input specifier/argument
        self.parser.add_argument('-P64in','--PHRED64input', action = 'store', nargs= '?', const=True, default=False, help = 'PHRED64 input') #PHRED 64 input specifier/argument
        self.parser.add_argument('-P64Bin','--PHRED64Binput', action= 'store', nargs='?',const=True, default=False, help= 'PHRED64 with B offset in quality values') #PHRED 64 B Offset input specifier/argument
        self.parser.add_argument('-P64SOLin', '--PHRED64SOLinput', action='store',nargs='?',const='True',default=False, help = 'PHRED64 with SOLEXA interpretation of Q score') #PHRED 64 Solexa Scoring input specifier/argument

        self.parser.add_argument('-P33out', '--PHRED33output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED33') #PHRED 33 output specifier/argument
        self.parser.add_argument('-P64out', '--PHRED64output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED64') #PHRED 33 output specifier/argument

        if inOpts is None: #handles command line input from stdin
            self.args = self.parser.parse_args() 
        else:
            self.args = self.parser.parse_args(inOpts) #parse only specific commands if specified