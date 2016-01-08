##!/usr/bin/env python3

# Open Q file and cellProc
f = open('Q','r')
c = open('/home/vitors/workspace/tutorials/coupled8/constant/fuel/polyMesh/cellRegionAddressing','r')

# Read its contents to one variable
lidof = f.read()
lidoc = c.read()

# Find the nonUniform List, or
# everything between parenthesis
findex = lidof.find("(")+1
lindex = lidof.find(")")
listaf = lidof[findex:lindex]

findexc = lidoc.find("(")+1
lindexc = lidoc.find(")")
listac = lidoc[findexc:lindexc]

# Split and sort the list
listf = listaf.split()
listf.sort()

listc = listac.split()
listc.sort()

# Join list back to string
listac = '\n'.join(listc)
listaf = '\n'.join(listf)

# Insert \n at both begin and end
listaf = '\n'+listaf+'\n'
listac = '\n'+listac+'\n'

print(listac)
print(listaf)


