
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import MutableSeq
# For input DNA sequence
dna=input('paste your DNA sequence: ')
dna=dna.upper()
dna=Seq(dna)

#ATGcount is needed to count start codons, so count an amount of possible coding sequences 
ATGcount=dna.count('ATG')
print(ATGcount, " ATG are found")
ii=0
#CDSs is a list of possible coding DNA sequences which start from ATG start codon to the end of a sequence
CDSs=[]
find_CDS=len(dna)
while ii<ATGcount:
    find_CDS=dna.rfind('ATG',end=find_CDS)
    CDS=str(dna[find_CDS:])        
    CDSs+=[CDS]
    #dna=dna[:find_CDS]
    
    ii=ii+1
    
print(CDSs)

#followed loop (str 30-42) is neeeded to modify sequences in CDSs making them mulpiple of 3
#(with responce to ByoPython warning, that if sequence is not multiple of 3 it can cause mistakes)
CDSs_full_codon=[]
for sequence in CDSs:
    lenght=float(len(sequence))
    
    if lenght%3!=0:
        sequence=MutableSeq(sequence)
        sequence=sequence+"N"
        
        if float(len(sequence))%3!=0:
            sequence=sequence+"N"
                     
    CDSs_full_codon+=[sequence]
print(CDSs_full_codon)    
# To get RNA and protein sequences of all of the CDSs from CDSs list
proteinforblast=[]
for dnaseq in CDSs_full_codon:
    dnaseq=str(dnaseq)
    dnaseq=Seq(dnaseq)
    rna=dnaseq.transcribe()
    protein=dnaseq.translate()
    print('\n==**===**=RNA=**===**==\n')
    print(rna)
    #
    print('\n>>>PROTEIN<<<\n')
    print(protein)
    print('\n')
    proteinforblast+=[protein]
print(proteinforblast)
#Protein sequences are transfered to NCBI BLAST; as a result we get a title of alignment and E value, alignments with E value < 0.0001 can be 
evalue=0.0001
ct=0
for protein_blast in proteinforblast:
    print('\n')
    print('>>>  >>  > For protein: ', str(protein_blast))
    result_handle = NCBIWWW.qblast("blastp", "pdb", protein_blast)
    blast_record=NCBIXML.read(result_handle)
    print: (ct, ' sequences are found:')
    for alignment in blast_record.alignments:
        print ('\n')
        
        print('For aligment: ', alignment.title)
        for hsp in alignment.hsps:
            ct+=1
            if hsp.expect < evalue:
                
                
                print ('\n')
                
                print('*********alignment**********')
                print('sequence: ', alignment.title)
                print('E value: ', hsp.expect)
                #print(hsp.query [0:75]+ '...')
                #print(hsp.match [0:75]+ '...')
                #print(hsp.sbjct [0:75]+ '...')
                print: (ct, ' sequences are found:')
            else:
                print('E value: ', hsp.expect)
                print: (ct, ' sequences are found:')