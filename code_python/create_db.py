import gffutils

# Tạo database từ GBFF (có thể cần chuyển đổi định dạng)
db = gffutils.create_db("data/genomic.gbff", "data/genome.db", force=True, keep_temp=False, 
                       disable_infer_transcripts=True, disable_infer_genes=True)
print("Database created successfully")