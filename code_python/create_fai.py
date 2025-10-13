from Bio import SeqIO

with open("data/genomic.fasta", "r") as f:
    sequences = SeqIO.parse(f, "fasta")
    with open("data/genomic.fasta.fai", "w") as fai:
        for record in sequences:
            fai.write(f"{record.id}\t{len(record.seq)}\t0\t{len(record.seq)+1}\t{len(record.seq)+1}\n")