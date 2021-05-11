from findUnique import FindUnique #type: ignore
testStr = 'CACUAAGA"LCUAUAUAGCACPAACCU∃UU6AGUUAGAGAUUGAGAGCCAU"UACUCUCCUUGGUGACCA'
seq = FindUnique(testStr)
#print(seq.powerset())

curSet = {'ATG','GTG','NNN'}

print(curSet-{'ATG','GTG'})