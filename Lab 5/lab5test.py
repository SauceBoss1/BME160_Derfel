import sequenceAnalysis
from importlib import reload
reload(sequenceAnalysis)

# gene = sequenceAnalysis.FastAreader('tass2.fa')

# sequence =''

# for header,seq in gene.readFasta():
#     sequence = ''.join(seq)

#testGene = 'ATGATGTAA'
seq=''
with open('tass2.fa','r') as f:
    for line in f:
        if line.startswith ('>'):
            header = line[1:].rstrip()
            seq = ''
        else :
            seq += ''.join(line.rstrip().split()).upper()
def reverseStr(seq):
    tempString = list(seq)

    chars = { 'A' : 'T', 'T': 'A', 'C' : 'G', 'G' : 'C'}
    tempString = reversed([chars.get(base,base) for base in tempString])
    seq = ''.join(tempString)
    return seq

print(len(reverseStr(seq)))
orf = sequenceAnalysis.OrfFinder(reverseStr(seq))
orf1=orf.finalORFlist(100, False)
#print(orf1)
for orf in sorted(orf1, key=lambda a:(a[2]), reverse= True):
    print(orf)
#print(len(sequence))

#print(sequence[75659:75760])



