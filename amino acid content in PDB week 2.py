

#TASK: Protein quaternary structure, percentage of each amino acid and type of amino acid that is more exposed
#Plan: 
# 1. Download PDB file from RCSB.org and parse it getting list of names of residues from Chain branch. 
# 2. Create a dictionary:
#   keys - aminoacids names (objects from standart_aa_names); 
#   values - the amounts of each residue in a structure of interest.
# 3. print data as a percentage & abb. of aminoacid that is more exposed

from Bio.PDB import *
import numpy as np


#Download PDB file from RCSB.org
pdb=PDBList()
pdb.retrieve_pdb_file("7RKW", pdir=".", file_format="pdb")
 
# Parsing the structure and obtaining the names of all residues as a list (from the "Chain" branch)
parser=PDBParser(PERMISSIVE=True, QUIET=True)
structure=parser.get_structure("7RKW", "pdb7RKW.ent")
aminoacids=[]
for model in structure:
      for chain in model:
          for residue in chain.get_list():
             is_aa(residue)
             if True:
                  residue = str(residue.resname)
                  aminoacids.append(residue)
                  
                
# Making a list of standard_aa_names, which originaly have tuple type
standard_aa_names=list(standard_aa_names)

# Creating a dictionary: keys - aminoacids names; values - the amounts of each residue in a structure of interest
aa_dict=dict()
for Name in standard_aa_names:
    aa_dict[Name]=0

for Aminoacid in aminoacids:
    if Aminoacid in aa_dict:
        aa_dict[Aminoacid]+=1
    else: continue

# Calculating the value of 100% amino acids by finding the total number of amino acids
total=0
for key in aa_dict.keys():
    total+=aa_dict[key]

# Calculating % of each residue and residue which is most exposed 
for key in aa_dict.keys():
    aa_content = np.round((aa_dict[key]/total)*100,2)
    print(aa_content, "% of", key)
aa_of_max = max(aa_dict, key = lambda k: aa_dict[k])
print("An amino acid that is more exposed:", aa_of_max)
