from flask import Flask, render_template, jsonify, request
import os
import subprocess
import json

app = Flask(__name__)

# -----------------------------------------------------------------
# 1. CẤU HÌNH VÀ TẢI DỮ LIỆU CHÚ GIẢI
# -----------------------------------------------------------------
PROCESSED_DATA_FILE = "data/processed_genome_data.json"
GENOME_DATA = {}
if os.path.exists(PROCESSED_DATA_FILE):
    with open(PROCESSED_DATA_FILE, "r", encoding='utf-8') as f:
        GENOME_DATA = json.load(f)
        print("Đã tải dữ liệu chú giải thành công.")
else:
    print(f"Lỗi: Không tìm thấy file {PROCESSED_DATA_FILE}. Vui lòng chạy process_gbff.py.")

# Danh sách các loài có sẵn (để truyền cho HTML)
AVAILABLE_SPECIES = list(GENOME_DATA.keys())

# Đường dẫn đến công cụ BLAST (CẦN KIỂM TRA LẠI)
BLASTN_PATH = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastn.exe"

# -----------------------------------------------------------------
# 2. ROUTES CHÍNH
# -----------------------------------------------------------------

@app.route('/')
def home():
    # Trang chủ và Tính năng Tìm kiếm Gen
    return render_template('index.html', species_list=AVAILABLE_SPECIES)

@app.route('/blast')
def blast():
    # Trang BLAST
    return render_template('blast.html', species_list=AVAILABLE_SPECIES)

# -----------------------------------------------------------------
# 3. API: TÌM KIẾM GEN/CONTIG
# -----------------------------------------------------------------
# Ví dụ: /api/gene_search?species=S_spontaneum&contig_id=JAKOGQ010001029.1
@app.route('/api/gene_search', methods=['GET'])
def gene_search():
    species = request.args.get('species')
    contig_id = request.args.get('contig_id')

    if species not in GENOME_DATA:
        return jsonify({"error": f"Loài '{species}' không có sẵn."}), 404
    
    result = GENOME_DATA[species].get(contig_id)
    
    if result:
        # Trả về kết quả từ dữ liệu JSON đã xử lý
        return jsonify({
            "species": species,
            "contig_id": contig_id,
            "description": result["description"],
            "sequence_length": result["sequence_length"],
            "features": result["features"] 
        })
    else:
        # Thử tìm kiếm theo tên/chú giải trong features (ví dụ đơn giản)
        found_contig_id = None
        for cid, data in GENOME_DATA[species].items():
            if contig_id.lower() in data["description"].lower():
                 found_contig_id = cid
                 break
            # Logic tìm kiếm phức tạp hơn có thể được thêm vào đây
        
        if found_contig_id:
             return jsonify(GENOME_DATA[species].get(found_contig_id))
        
        return jsonify({"error": f"Không tìm thấy Contig/Gen ID '{contig_id}' trong loài {species}."}), 404

# -----------------------------------------------------------------
# 4. API: PHÂN TÍCH BLAST
# -----------------------------------------------------------------
@app.route('/api/blast/search', methods=['POST'])
def blast_search():
    target_species = request.form.get('species') 
    sequence = request.form['sequence']
    evalue = request.form.get('evalue', '10')
    max_hits = request.form.get('max_hits', '10')

    if not os.path.exists(BLASTN_PATH):
        return jsonify({"error": "BLAST executable không tìm thấy. Vui lòng kiểm tra đường dẫn trong app.py"}), 500
        
    if target_species not in AVAILABLE_SPECIES:
        return jsonify({"error": "Loài mục tiêu không hợp lệ."}), 400
        
    # Xây dựng đường dẫn đến BLAST DB đã tạo
    blast_db = os.path.join("data", target_species, "genome_blast_db") 
    
    # Kiểm tra BLAST DB đã tồn tại chưa
    if not os.path.exists(blast_db + ".nhr"):
        # Lỗi này cũng sẽ bị bắt
        return jsonify({"error": f"BLAST database cho {target_species} không tìm thấy..."}), 500

    try:
        # 1. Tạo file truy vấn tạm thời
        query_file = "data/query.fasta"
        with open(query_file, "w") as f:
            lines = sequence.strip().split("\n")
            if not lines[0].startswith(">"):
                f.write(">query\n" + sequence.strip() + "\n")
            else:
                f.write("\n".join(lines) + "\n")
        
        # 2. Chạy blastn.exe
        result = subprocess.run(
            [BLASTN_PATH, 
             "-db", blast_db, 
             "-query", query_file, 
             "-outfmt", "6", 
             "-evalue", evalue, 
             "-max_target_seqs", max_hits],
            capture_output=True, text=True, check=True # **check=True RẤT QUAN TRỌNG**
        )
        
        # 3. Trả về kết quả (nếu thành công)
        return jsonify({"result": result.stdout.splitlines() if result.stdout else ["No hits found"]})

    except subprocess.CalledProcessError as e:
        # THAY ĐỔI: In lỗi ra console Flask và trả về thông báo chi tiết
        print(f"\n--- BLASTN EXECUTION ERROR ---\n")
        print(f"Error command: {e.cmd}")
        print(f"Error output (stderr): {e.stderr}")
        print(f"--------------------------------\n")
        return jsonify({"error": f"Lỗi thực thi BLASTN: {e.stderr.strip()}"}), 500
    except Exception as e:
        print(f"General error: {str(e)}")
        return jsonify({"error": f"Lỗi chung: {str(e)}"}), 500

# -----------------------------------------------------------------
# 5. CHẠY ỨNG DỤNG
# -----------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)