from flask import Flask, render_template, jsonify, request
import os
import subprocess
import json
import re 
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN as mt
from Bio.SeqUtils import gc_fraction 

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

MOCK_SPECIES = "S_spontaneum"
MOCK_CONTIG = "JAKOGQ010001029.1"
MOCK_DATA = {
    "description": "hypothetical protein, putative (Mock Data)",
    "sequence_length": 5241,
    "features": [
        {"type": "source", "location_start": 1, "location_end": 5241, "strand": "+", "product": None, "gene": None, "locus_tag": None},
        {"type": "CDS", "location_start": 100, "location_end": 1500, "strand": "+", "product": "Sucrose Synthase 1 (SuSy-1)", "gene": "SuSy1", "locus_tag": "SSpTAG00100"},
        {"type": "tRNA", "location_start": 2000, "location_end": 2075, "strand": "-", "product": "tRNA-Ala", "gene": None, "locus_tag": None},
        {"type": "CDS", "location_start": 3500, "location_end": 4500, "strand": "+", "product": "Cellulose Synthase A1", "gene": "CesA1", "locus_tag": "SSpTAG00101"}
    ]
}

# Hợp nhất dữ liệu mẫu với dữ liệu chính
if MOCK_SPECIES not in GENOME_DATA:
    GENOME_DATA[MOCK_SPECIES] = {}
GENOME_DATA[MOCK_SPECIES][MOCK_CONTIG] = MOCK_DATA


AVAILABLE_SPECIES = list(GENOME_DATA.keys())
BLASTN_PATH = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastn.exe"

# -----------------------------------------------------------------
# 2. ROUTES CHÍNH
# -----------------------------------------------------------------

@app.route('/')
def home():
    return render_template('index.html', species_list=AVAILABLE_SPECIES)

@app.route('/blast')
def blast():
    return render_template('blast.html', species_list=AVAILABLE_SPECIES)

@app.route('/primer')
def primer_design_page():
    return render_template('primer.html')

@app.route('/crispr_designer')
def crispr_designer():
    """Hiển thị module thiết kế CRISPR Tefor."""
    return render_template('crispr_tefor_designer.html')

@app.route('/gene_view')
def gene_view():
    species = request.args.get('species')
    contig_id = request.args.get('contig_id')
    
    if not species or not contig_id:
        return render_template('error.html', message="Thiếu thông tin loài hoặc Contig ID."), 400

    # Lấy dữ liệu chi tiết từ GENOME_DATA
    contig_data = GENOME_DATA.get(species, {}).get(contig_id)

    # KHẮC PHỤC: Kiểm tra chặt chẽ hơn. Nếu URL có Contig ID nhưng nó không tồn tại 
    # trong dictionary GENOME_DATA (kể cả sau khi thêm mock data), thì báo lỗi 404.
    if not contig_data:
        # Nếu không tìm thấy, Flask sẽ trả về 404 (Lỗi Not Found)
        return render_template('error.html', message=f"Không tìm thấy dữ liệu chi tiết cho Contig ID: {contig_id} trong loài {species}. Vui lòng kiểm tra lại ID."), 404

    # Truyền dữ liệu chi tiết cho template gene_view.html
    return render_template('gene_view.html', species=species, contig_id=contig_id, data=contig_data)


# -----------------------------------------------------------------
# 3. API: TÌM KIẾM GEN/CONTIG 
# -----------------------------------------------------------------
@app.route('/api/gene_search', methods=['GET'])
def gene_search():
    species = request.args.get('species')
    contig_id = request.args.get('contig_id')

    if species not in GENOME_DATA:
        return jsonify({"error": f"Loài '{species}' không có sẵn."}), 404
    
    result = GENOME_DATA[species].get(contig_id)
    
    if result:
        # THÊM: Cần đảm bảo rằng khi tìm thấy, nó trả về tất cả các trường cần thiết.
        return jsonify({
            "species": species,
            "contig_id": contig_id,
            "description": result["description"],
            "sequence_length": result["sequence_length"],
            "features": result["features"] 
        })
    else:
        found_contig_id = None
        for cid, data in GENOME_DATA[species].items():
            # Thử tìm theo mô tả
            if contig_id.lower() in data["description"].lower():
                 found_contig_id = cid
                 break
        
        if found_contig_id:
            result = GENOME_DATA[species].get(found_contig_id)
            return jsonify({
                "species": species,
                "contig_id": found_contig_id,
                "description": result["description"],
                "sequence_length": result["sequence_length"],
                "features": result["features"] 
            })
        
        return jsonify({"error": f"Không tìm thấy Contig/Gen ID '{contig_id}' trong loài {species}."}), 404

# -----------------------------------------------------------------
# 4. API: PHÂN TÍCH BLAST (Giữ nguyên)
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
        
    blast_db = os.path.join("data", target_species, "genome_blast_db") 
    
    if not os.path.exists(blast_db + ".nhr"):
        return jsonify({"error": f"BLAST database cho {target_species} không tìm thấy. Vui lòng chạy create_blast_dbs.py."}), 500

    try:
        query_file = "data/query.fasta"
        with open(query_file, "w") as f:
            lines = sequence.strip().split("\n")
            if not lines[0].startswith(">"):
                f.write(">query\n" + sequence.strip() + "\n")
            else:
                f.write("\n".join(lines) + "\n")
        
        result = subprocess.run(
            [BLASTN_PATH, 
             "-db", blast_db, 
             "-query", query_file, 
             "-outfmt", "6", 
             "-evalue", evalue, 
             "-max_target_seqs", max_hits],
            capture_output=True, text=True, check=True
        )
        
        return jsonify({"result": result.stdout.splitlines() if result.stdout else ["No hits found"]})



    except subprocess.CalledProcessError as e:
        print(f"\n--- BLASTN EXECUTION ERROR ---\n")
        print(f"Error command: {e.cmd}")
        print(f"Error output (stderr): {e.stderr}")
        print(f"--------------------------------\n")
        return jsonify({"error": f"Lỗi thực thi BLASTN: {e.stderr.strip()}"}), 500
    except Exception as e:
        print(f"General error: {str(e)}")
        return jsonify({"error": f"Lỗi chung: {str(e)}"}), 500

# -----------------------------------------------------------------
# 5. API: THIẾT KẾ PRIMER PCR
# -----------------------------------------------------------------

def analyze_primer(seq):
    """Tính toán Tm và %GC cho một trình tự primer."""
    
    # LÀM SẠCH VÀ CHUẨN HÓA TRÌNH TỰ: Chỉ giữ lại các ký tự A, T, C, G
    cleaned_seq = re.sub(r'[^ATCGatcg]', '', seq).upper()
    primer = Seq(cleaned_seq)
    
    # SỬA LỖI TM: Thêm strict=False để bỏ qua các base không xác định
    # Thay thế dnac=250.0 bằng c_seq=250.0 cho Biopython 1.85+
    tm = mt(primer, Na=50.0, c_seq=250.0, strict=False) 
    
    # Tính %GC
    gc_percent_fraction = gc_fraction(primer)
    gc_percent = gc_percent_fraction * 100 

    return {
        "sequence": str(primer),
        "length": len(primer),
        "tm": f"{tm:.2f} °C",
        "gc_percent": f"{gc_percent:.2f} %"
    }


@app.route('/api/primer/design', methods=['POST'])
def primer_design_api():
    # Lấy trình tự thô từ form
    full_sequence_raw = request.form['full_sequence']
    
    # LÀM SẠCH CHUỖI TỪ ĐẦU (chỉ giữ lại ATCG)
    full_sequence = re.sub(r'[^ATCGatcg]', '', full_sequence_raw).upper()

    if len(full_sequence) < 50:
        return jsonify({"error": "Trình tự đầu vào quá ngắn. Cần tối thiểu 50 bp để thiết kế primer."}), 400

    # Giả định primer Forward ở đầu và primer Reverse ở cuối.
    # Chiều dài primer lý tưởng: 20 bp
    primer_len = 20

    # 1. Primer Forward (1-20)
    primer_fwd_seq = full_sequence[:primer_len]
    fwd_analysis = analyze_primer(primer_fwd_seq)

    # 2. Primer Reverse (ngược và bổ sung, 20 bp cuối)
    # Lấy 20 bp cuối: full_sequence[-primer_len:]
    # Đảo ngược và bổ sung: Seq(chuỗi).reverse_complement()
    primer_rev_seq = str(Seq(full_sequence[-primer_len:]).reverse_complement())
    rev_analysis = analyze_primer(primer_rev_seq)
    
    # 3. Tính toán Product Size
    product_size = len(full_sequence)

    return jsonify({
        "success": True,
        "product_size": product_size,
        "forward": fwd_analysis,
        "reverse": rev_analysis
    })


@app.route('/jbrowse_view')
def jbrowse_view():
    # Lấy các tham số từ URL
    chrom = request.args.get('chrom')
    start = request.args.get('start')
    end = request.args.get('end')
    name = request.args.get('name')
    
    # Hàm render_template sẽ tìm và hiển thị tệp HTML
    # Bạn cần đảm bảo tệp jbrowse_view.html nằm trong thư mục 'templates'
    return render_template('jbrowse_view.html', 
                           chrom=chrom, start=start, end=end, name=name)

# -----------------------------------------------------------------
# 6. CHẠY ỨNG DỤNG
# -----------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)
