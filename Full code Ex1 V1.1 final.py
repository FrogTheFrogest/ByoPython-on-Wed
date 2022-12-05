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

#CDSs is a list of possible coding DNA sequences which start from ATG start codon to the end of a sequence
CDSs=[]
ii=0
find_CDS=len(dna)
while ii<ATGcount:
    find_CDS=dna.rfind('ATG',end=find_CDS)
    CDS=str(dna[find_CDS:])        
    CDSs+=[CDS]
    ii=ii+1
    
#followed loop is neeeded to modify sequences in CDSs making them mulpiple of 3
#(with responce to ByoPython warning, that if sequence is not multiple of 3 it can cause mistakes)
CDSs_full_codon=[]
for sequence in CDSs:
    lenght=float(len(sequence))
    if lenght%3!=0:
        sequence=MutableSeq(sequence)
        sequence=sequence+"N    
        if float(len(sequence))%3!=0:
            sequence=sequence+"N"                
    CDSs_full_codon+=[sequence]
    

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

#Protein sequences are transfered to NCBI BLAST for processing; as a result we get a title of alignment, E value, and a bit score of alignments if e-value < 0.0001
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
                print('Sequence: ', alignment.title)
                print('E value: ', hsp.expect)
                print("Bit score: ", hsp.bits)                          
