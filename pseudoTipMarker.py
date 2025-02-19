from ete3 import Tree
import pandas as pd
import numpy as np

#first read in the timeTree tree made from the list generated in Model_Codes.R
t = Tree('tree.nwk',format=1)
tipNames = []
for leaf in t:
    tipNames.append(leaf.name)


#next read in the comparative data
#we'll do this with pandas
d = pd.read_csv("Database.csv")

d = d.drop(0) # Drop the first row which is notation

spNames = d["Scientific_name"].tolist() #the method here coerces this from an array
#the names in the tree have underscores instead of spaces
#we need to change that to find matches in our comparative data
for i in range(0,len(spNames)-1):
    nm = spNames[i]
    nnm = nm.replace(' ', '_')
    spNames[i] = nnm

uniqueNames = set(spNames) # Find unique species
nameCounts = {}
for i in uniqueNames: # Count the occurance of each species
    x = spNames.count(i)
    nameCounts[i] = x
#nameCounts is a dictionary with each name as a key, with its count as the value
#print(spNames)

# The shortest real branch length is 3.688320, 0.000001 should be small enough for population branches.
#branches=[]
#for i in uniqueNames:
#	if i in tipNames:
#		x = t.get_leaves_by_name(i)[0]
#		branches.append(x.dist)
#branches = pd.Series(branches)
#branch_summary = branches.describe()

#now for multi-model species, we add new branches to the tree
branchLength = 0.000001 # a minimum branch length between populations
for i in uniqueNames:
	if i in tipNames:
		if nameCounts[i] > 1:
			x = t.get_leaves_by_name(i)[0]
			x.dist = x.dist-branchLength
			mom = x
	    
			counter = 1
			for j in range(1,nameCounts[i]+1):
				nm = i + str(counter)
				mom.add_child(name=nm, dist=branchLength)
				#print(nm)
				counter = counter + 1

fruitTipNames=[]            
fruitTipNames = []
for leaf in t:
    fruitTipNames.append(leaf.name)
print("Tip number: " + str(len(tipNames)) + " -> " + str(len(fruitTipNames)))


t.write(format=1, outfile="treePseudoTips.nwk")  

#for multi-model species, species names are not updated in the Database but in Model_Codes.R





