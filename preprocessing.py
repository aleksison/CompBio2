from ctypes import alignment
from Bio import AlignIO
from Bio import SeqIO

file = open('ALIGN.ann')
alignments = AlignIO.read(file, "stockholm")

aligns = []

for record in alignments:
    #print(f">{record.id}")
    #print(record.seq)
    align = SeqIO.SeqRecord(record.seq, id=record.id)
    #print(align)
    # align look like this:
    # ID: A0R1P9_MYCS2/28-74
    # Name: <unknown name>
    # Description: <unknown description>
    # Number of features: 0
    # Seq('----------------------ILSVA-E-R--L-L--A-------R----N--...---')
    
    aligns.append(align)

sth_file = "newgen.sth"

with open(sth_file, "w") as output_handle:
    SeqIO.write(aligns, output_handle, "stockholm")

# Step 2: eliminate higly gapped columns
# we want to extract the sequences
for record in aligns:
    print(record.seq)