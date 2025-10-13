from Bio import SeqIO

# Đường dẫn đến file GBFF
gbff_file = "data/genomic.gbff"

# Đọc file GBFF
try:
    records = list(SeqIO.parse(gbff_file, "genbank"))
    print(f"Found {len(records)} records in {gbff_file}")

    # In thông tin cơ bản của từng record
    for record in records:
        print(f"\nRecord ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Sequence length: {len(record.seq)}")
        print(f"Number of features: {len(record.features)}")

        # In một số feature (ví dụ: gene, CDS)
        for feature in record.features[:5]:  # Giới hạn 5 feature đầu
            print(f"Feature type: {feature.type}, Location: {feature.location}")
            if feature.type == "gene":
                gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                print(f"Gene name: {gene_name}")

except Exception as e:
    print(f"Error reading file: {str(e)}")