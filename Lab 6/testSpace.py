from findUnique import FindUnique #type: ignore
#testStr = 'CACUAAGA"LCUAUAUAGCACPAACCU∃UU6AGUUAGAGAUUGAGAGCCAU"UACUCUCCUUGGUGACCA'
#seq = FindUnique(testStr)
#print(seq.powerset())

curSet = {'ATG','GTG','NNN'}

print(curSet-{'ATG','GTG'})
genes = {'gene1': ['ATGGCG', {'ATG', 'GCG', 'gg','g', 'cc','c'}], 'gene2':['GGG',{'gg','g'}], 'gene3':['ccc',{'cc','c'}]}
un1 = set()
un2 = {'1','2','3'}
print(un1 - un2)

#print(genes['gene1'][1])

def uniqueFinder(curHeader):
    allSets = []
    for header, seq in genes.items():
        if curHeader != header:     
            allSets.append(genes[header][1])
    return genes[curHeader][1] - set().union(*allSets)

print('final: ',uniqueFinder('gene1'))

set1 = {'AAAACG','ACG'}

def essential(inSet):
    nonEssential = []
    uniqueSet = list(inSet)
    compareSet = list(uniqueSet.copy())
    for subset in uniqueSet:
        compareSet.remove(subset)
        for compareSub in compareSet:
            if compareSub in subset:
                nonEssential.append(subset)
                print('nonEssential: ', subset)
        compareSet.append(subset)
    return set(uniqueSet) - set(nonEssential)

print('essential', essential(set1))

seq = 'CCCACTGACTGACTGGGACTGGTT'
print(seq.find('GGG'))

test11 = { 0: 'aaa', 1: '2222'}
print(test11[0])

def findAllOccurances( subSeq, mainSeq):
    posFound = []
    i = mainSeq.find(subSeq)
    while i != -1:
        posFound.append(i)
        i = mainSeq.find(subSeq, i+1)
    return posFound
print(str(findAllOccurances('AAAAUA', 'UGGUACUU"GUUUAAAAUAAAAUAAAUGAUUUCGACPCAUUAGAUUAUGAUUUAAUUCAUAAUUACCAACCA')))

emptyObj = FindUnique()
print(sorted(emptyObj.powerset('ABCDEF')))