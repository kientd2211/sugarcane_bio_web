from Bio import SeqIO
records = SeqIO.parse("data/genomic.gbff", "genbank")
SeqIO.write(records, "data/genomic.fasta", "fasta")