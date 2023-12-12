from ctypes import alignment
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Step 2: eliminate higly gapped columns
fileshort = open('shortdoc.sth')
alignments = AlignIO.read(fileshort, "stockholm")

#extract sequences
i = 0
countgap = 0
for i in range(0,50): #107 length of the shortest sequence
  for record in alignments:
      #record.seq
      #print(align)
      #acces amino acid
      if record.seq[i] == '-':
          countgap += 1
  #if gaps represent at least 30% we get rid of the column
  if countgap >= 50*0.3:
      for record in alignments:
          # delete by slicing?
          record.seq = record.seq[:i] + record.seq[i+1:]
          # record.seq[i] delete

output_file = "aminoacids.sth"

with open(output_file, "w") as output_handle:
    SeqIO.write(alignments, output_handle, "stockholm")
