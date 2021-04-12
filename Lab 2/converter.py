#!/usr/bin/env python3 
# Name: Derfel Terciano (dtercian) 
# Group Members: None

'''
converter.py


This programs takes in a user input and outputs its corresponding amino acid
    Inputs can be RNA or DNA codons, Amino Acid abbreviations, or one-letter amino acid abbreviations
    Outputs are vary depending on the input. The following shows the general input and the expected output:
        -RNA or DNA codon => Three letter Amino Acid Abbreviation
        -Three letter Amino Acid Abbreviation => corresponding one letter abbreviation of AA
        -One letter Amino Acid Abbreviation => corresponding three letter Amino Acid Abbreviation
    
    Example:
    if I enter “ATG” (without quotes), then the program will output:
        ATG = MET
    if I enter “UAG” (without quotes), then the program will output:
        UAG = ---
    if I enter “E” (without quotes), then the program will output:
        E = GLU
    if I enter “Asp” (without quotes), then the program will output:
        ASP = D

Assumptions:
    I'm going to assume that only amino acid related inputs (like mentioned above) are the only inputs otherwise return unknown
'''


short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

#I use a seperate function to compare the user input in order to break the problem up into smaller pieces
def codonInterpreter (inputString):
    '''
    This function determins the correct dictionary output based on what input is given

    If a codon is incorrect or unrecogninzable 'unknown' is returned

    inputs: a string
    output: a formatted string
    '''
    inputString=inputString.upper() #makes input all uppercase

    # a series of if and elif statments that determines whether or not the input is in one of the four dictionaries
    if inputString in dnaCodonTable: #checks if input is a DNA codon
        return f'{inputString} = ' + f'{dnaCodonTable.get(inputString,"unknown").upper()}'
    elif inputString in RNA_codon_table: #checks if input is a RNA codon
        return f'{inputString} = ' + f'{RNA_codon_table.get(inputString,"unknown").upper()}'
    elif inputString in short_AA: #checks if input is a three-letter AA abbreviation
        return f'{inputString} = ' + f'{short_AA.get(inputString.upper(),"unknown")}'
    elif inputString in long_AA: #checks if input is a three-letter AA abbreviation
        return f'{inputString} = ' + f'{long_AA.get(inputString,"unknown")}'
    else:
        return f'{inputString} = unknown'

def main():
    '''Asks user for a string and returns the correct output'''
    print(codonInterpreter(input('Enter String: ')))

main()