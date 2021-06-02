#!/usr/bin/env python3
# Collaborators: Sagarika Kannoly (skannoly) & Derfel Terciano (dtercian)

class CommandLine():

    def __init__(self, inOpts = None):
        import argparse
        self.parser = argparse.ArgumentParser(description = 'handles the fastQ command input',
                                                add_help= True,
                                                prefix_chars= '-',
                                                usage = '%(prog)s [options] -option1[default] <input >output'
                                                )

        self.parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False, help = 'PHRED33 input')
        self.parser.add_argument('-P64in','--PHRED64input', action = 'store', nargs= '?', const=True, default=False, help = 'PHRED64 input')
        self.parser.add_argument('-P64Bin','--PHRED64Binput', action= 'store', nargs='?',const=True, default=False, help= 'PHRED64 with B offset in quality values')
        self.parser.add_argument('-P64SOLin', '--PHRED64SOLinput', action='store',nargs='?',const='True',default=False, help = 'PHRED64 with SOLEXA interpretation of Q score')

        self.parser.add_argument('-P33out', '--PHRED33output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED33')
        self.parser.add_argument('-P64out', '--PHRED64output', action='store',nargs='?',const='True',default=False, help = 'output to PHRED64')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)