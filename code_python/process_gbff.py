from Bio import SeqIO
import json
import os

# Danh sách các loài và thư mục dữ liệu
species_list = ["S_spontaneum", "S_officinarum", "S_bicolor", "O_sativa"]
data_dir = "data"
processed_data = {}

print("--- BẮT ĐẦU XỬ LÝ CHÚ GIẢI GBFF ---")

for species in species_list:
    gbff_file = os.path.join(data_dir, species, "genome.gbff")
    
    if not os.path.exists(gbff_file):
        print(f"CẢNH BÁO: Bỏ qua {species}. Không tìm thấy file GBFF.")
        continue

    print(f"\nĐang xử lý GBFF của loài: {species}...")
    species_data = {}
    try:
        # Tải dữ liệu GBFF bằng Biopython
        for record in SeqIO.parse(gbff_file, "genbank"):
            features_list = []
            
            # Trích xuất Features (gene, CDS, tRNA, rRNA...)
            for feature in record.features:
                if feature.type in ["gene", "CDS", "tRNA", "rRNA"]:
                    features_list.append({
                        "type": feature.type,
                        "location": str(feature.location),
                        "qualifiers": dict(feature.qualifiers) # Lưu các thông tin chi tiết (tên gen, sản phẩm,...)
                    })
            
            # Lưu dữ liệu theo Contig ID
            species_data[record.id] = {
                "description": record.description,
                "sequence_length": len(record.seq),
                "features": features_list
            }
        
        processed_data[species] = species_data
        print(f"   -> Hoàn thành xử lý {len(species_data)} Contigs/Chromosomes.")
        
    except Exception as e:
        print(f"   !!! LỖI khi xử lý GBFF của {species} !!!")
        print(f"   Error: {e}")

# Lưu toàn bộ dữ liệu đã xử lý vào một file JSON duy nhất
output_json_path = os.path.join(data_dir, "processed_genome_data.json")
with open(output_json_path, "w", encoding='utf-8') as f:
    json.dump(processed_data, f, indent=4, ensure_ascii=False)

print(f"\n--- HOÀN THÀNH. Dữ liệu đã lưu tại: {output_json_path} ---")