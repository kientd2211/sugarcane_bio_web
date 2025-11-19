import os
import subprocess
import json
import re 
from flask import Flask, render_template, jsonify, request, g
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN as mt
from Bio.SeqUtils import gc_fraction 
import psycopg2
import psycopg2.extras
from werkzeug.utils import secure_filename
from load_data import process_single_species_directory, connect_db

app = Flask(__name__)

# -----------------------------------------------------------------
# 1. CẤU HÌNH DATABASE VÀ DỮ LIỆU
# -----------------------------------------------------------------

# --- Cấu hình Database ---
# Lấy thông tin từ file load_data.py của bạn
DB_HOST = "localhost"
DB_PORT = "5432"  
DB_NAME = "gene_database"  
DB_USER = "kien"
DB_PASS = ""  # Mật khẩu rỗng
# -------------------------

def get_dynamic_genome_groups():
    """
    Lấy cấu trúc nhóm loài trực tiếp từ CSDL thay vì hard-code.
    """
    db = get_db()
    if db is None: return {}
    
    groups = {}
    try:
        with db.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            cursor.execute("SELECT name, common_name, group_name FROM species ORDER BY group_name")
            rows = cursor.fetchall()
            
            for row in rows:
                g_name = row['group_name'] if row['group_name'] else "Chưa phân loại"
                entry = {
                    "id": row['name'],
                    "name": row['common_name']
                }
                
                if g_name not in groups:
                    groups[g_name] = []
                groups[g_name].append(entry)
                
    except Exception as e:
        print(f"Lỗi lấy danh sách loài: {e}")
        
    return groups

# --- CẤU HÌNH UPLOAD ---
UPLOAD_FOLDER = 'data'
ALLOWED_EXTENSIONS = {'fasta', 'fna', 'fa', 'gbff', 'gb', 'gz'}
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# --- Cấu hình dữ liệu ---
GENOME_GROUPS = {
    "Mía (Saccharum & Sorghum)": [
        {"id": "S_spontaneum", "name": "Saccharum spontaneum (Mía hoang)"},
        {"id": "S_officinarum", "name": "Saccharum officinarum (Mía đường)"},
        {"id": "S_bicolor", "name": "Sorghum bicolor (Lúa miến)"}
    ],
    "Lúa (Oryza)": [
        {"id": "O_sativa", "name": "Oryza sativa Japonica (Lúa Nhật)"},
        {"id": "oryza_sativa_indica", "name": "Oryza sativa Indica (Lúa Việt/Ấn)"}
    ],
    "Cà phê (Coffee)": [
        {"id": "coffee_canephora_cultivar", "name": "Coffea canephora (Giống - Cultivar)"},
        {"id": "coffee_canephora_strain", "name": "Coffea canephora (Biến dị - Strain)"}
    ],
    "Nhân Sâm (Ginseng)": [
        {"id": "panax_ginseng", "name": "Panax ginseng (Nhân sâm)"}
    ],
    "Thanh Long (Dragon Fruit)": [
        {"id": "selenicereus_undatus", "name": "Selenicereus undatus (Thanh long)"}
    ],
    "Sầu Riêng (Durian)": [
        {"id": "durio_zibethinus", "name": "Durio zibethinus (Sầu riêng)"}
    ]
}

SPECIES_MAP = {}
for group, genomes in GENOME_GROUPS.items():
    for genome in genomes:
        SPECIES_MAP[genome['id']] = genome['name']

AVAILABLE_SPECIES = list(SPECIES_MAP.keys())

BLASTN_PATH = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastn.exe"

# -----------------------------------------------------------------
# 2. QUẢN LÝ KẾT NỐI CSDL VỚI FLASK
# -----------------------------------------------------------------

def get_db():
    """
    Mở một kết nối CSDL mới nếu chưa có cho context (g) này.
    """
    if 'db' not in g:
        try:
            g.db = psycopg2.connect(
                dbname=DB_NAME,
                user=DB_USER,
                password=DB_PASS,
                host=DB_HOST,
                port=DB_PORT
            )
        except psycopg2.OperationalError as e:
            print(f"LỖI KẾT NỐI CSDL: {e}")
            g.db = None
    return g.db

@app.teardown_appcontext
def close_db(e=None):
    """
    Đóng kết nối CSDL khi request kết thúc.
    """
    db = g.pop('db', None)
    if db is not None:
        db.close()

# -----------------------------------------------------------------
# 3. ROUTES CHÍNH (Đã cập nhật để dùng CSDL)
# -----------------------------------------------------------------

@app.route('/')
def home():
    # Lấy dữ liệu động từ DB
    dynamic_groups = get_dynamic_genome_groups()
    return render_template('index.html', genome_groups=dynamic_groups)

@app.route('/blast')
def blast():
    return render_template('blast.html', species_list=AVAILABLE_SPECIES)

@app.route('/primer')
def primer_design_page():
    return render_template('primer.html')

@app.route('/crispr_designer')
def crispr_designer():
    return render_template('crispr_tefor_designer.html')

@app.route('/gene_view')
def gene_view():
    # Lấy tham số từ URL
    species_folder = request.args.get('species') # vd: 'S_spontaneum'
    contig_id = request.args.get('contig_id')  # vd: 'JAKOGQ010001029.1'
    
    if not species_folder or not contig_id:
        return render_template('error.html', message="Thiếu thông tin loài hoặc Contig ID."), 400

    # Lấy kết nối CSDL
    db = get_db()
    if db is None:
        return render_template('error.html', message="Lỗi: Không thể kết nối đến cơ sở dữ liệu."), 500

    contig_data = None
    try:
        # Dùng cursor_factory để trả về kết quả dạng dict (giống JSON)
        with db.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            
            # 1. Truy vấn thông tin contig
            cursor.execute("SELECT * FROM contigs WHERE id = %s", (contig_id,))
            contig = cursor.fetchone()

            if not contig:
                return render_template('error.html', message=f"Không tìm thấy Contig ID: {contig_id}."), 404

            # 2. Truy vấn tất cả features của contig đó
            cursor.execute("SELECT * FROM features WHERE contig_id = %s ORDER BY start_pos", (contig_id,))
            features = cursor.fetchall()
            
            # 3. Gói dữ liệu lại để gửi cho template
            # Biến 'contig' (kiểu DictRow) đã chứa 'description', 'sequence_length', v.v.
            contig_data = dict(contig) 
            contig_data['features'] = [dict(f) for f in features] # Chuyển list features sang dict

    except (Exception, psycopg2.Error) as error:
        print(f"Lỗi truy vấn /gene_view: {error}")
        return render_template('error.html', message=f"Lỗi khi truy vấn CSDL: {error}"), 500

    # 4. Render template với dữ liệu từ CSDL
    return render_template('gene_view.html', 
                           species=species_folder, 
                           contig_id=contig_id, 
                           data=contig_data)


# --- ROUTE MỚI: TRANG UPLOAD ---
@app.route('/upload', methods=['GET', 'POST'])
def upload_data():
    if request.method == 'POST':
        # 1. Lấy thông tin form
        species_id = request.form['species_id'].strip() 
        common_name = request.form['common_name'].strip()
        group_name = request.form['group_name'].strip()
        
        # 2. Kiểm tra file
        if 'fasta_file' not in request.files:
            return "Chưa chọn file FASTA", 400
        file_fasta = request.files['fasta_file']
        file_gbff = request.files.get('gbff_file')

        if species_id and file_fasta:
            # 3. Tạo thư mục
            target_dir = os.path.join(app.config['UPLOAD_FOLDER'], species_id)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            
            # 4. Lưu file
            fasta_filename = "genome.fasta" # Đổi tên chuẩn để dễ xử lý
            file_fasta.save(os.path.join(target_dir, fasta_filename))
            
            if file_gbff:
                gbff_filename = "genome.gbff"
                file_gbff.save(os.path.join(target_dir, gbff_filename))

            # 5. CẬP NHẬT DATABASE (Metadata)
            db = get_db()
            try:
                with db.cursor() as cursor:
                    cursor.execute(
                        """
                        INSERT INTO species (name, common_name, group_name) 
                        VALUES (%s, %s, %s)
                        ON CONFLICT (name) DO UPDATE SET 
                        common_name = EXCLUDED.common_name,
                        group_name = EXCLUDED.group_name
                        """,
                        (species_id, common_name, group_name)
                    )
                db.commit()
            except Exception as e:
                return f"Lỗi CSDL: {e}", 500
            # 6. CHẠY XỬ LÝ DỮ LIỆU
            try:
                # Mở kết nối DB riêng cho tác vụ này
                task_conn = connect_db()
                if task_conn:
                    process_single_species_directory(task_conn, target_dir, species_id, common_name)
                    
                    task_conn.close()
                    print("Xử lý upload thành công!")
            except Exception as e:
                print(f"Lỗi xử lý nền: {e}")

            return render_template('success.html', message="Đã tải lên và lưu cấu hình. Vui lòng chạy script xử lý dữ liệu ở server để nạp chi tiết vào DB (hoặc chờ xử lý nền).")

    return render_template('upload.html')


# -----------------------------------------------------------------
# 4. API: TÌM KIẾM GEN/CONTIG (Đã cập nhật để dùng CSDL)
# -----------------------------------------------------------------
@app.route('/api/gene_search', methods=['GET'])
def gene_search():
    species_folder = request.args.get('species') # vd: 'S_spontaneum'
    search_term = request.args.get('contig_id')  # Tên gen, ID, hoặc mô tả

    if not species_folder or not search_term:
        return jsonify({"error": "Thiếu 'species' hoặc 'contig_id' (search_term)."}), 400

    if species_folder not in SPECIES_MAP:
        return jsonify({"error": f"Loài '{species_folder}' không hợp lệ."}), 404
    
    # Lấy tên loài đầy đủ từ map để truy vấn CSDL
    species_name = SPECIES_MAP[species_folder]

    db = get_db()
    if db is None:
        return jsonify({"error": "Lỗi kết nối CSDL."}), 500

    contig = None
    try:
        with db.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            
            # 1. Lấy species_id từ CSDL
            cursor.execute("SELECT id FROM species WHERE name = %s", (species_name,))
            species_row = cursor.fetchone()
            if not species_row:
                return jsonify({"error": f"Không tìm thấy loài {species_name} trong CSDL."}), 404
            species_id = species_row['id']

            # 2. Thử tìm bằng Contig ID chính xác (ưu tiên)
            cursor.execute(
                "SELECT * FROM contigs WHERE id = %s AND species_id = %s",
                (search_term, species_id)
            )
            contig = cursor.fetchone()

            # 3. Nếu không tìm thấy, thử tìm bằng mô tả (description)
            if not contig:
                # Dùng ILIKE để tìm kiếm không phân biệt hoa/thường
                search_pattern = f"%{search_term}%" 
                cursor.execute(
                    "SELECT * FROM contigs WHERE description ILIKE %s AND species_id = %s LIMIT 1",
                    (search_pattern, species_id)
                )
                contig = cursor.fetchone()

            # 4. Nếu vẫn không tìm thấy, thử tìm bằng tên gen hoặc locus_tag trong 'features'
            if not contig:
                # Truy vấn này có thể hơi chậm nếu bảng features rất lớn
                cursor.execute(
                    """
                    SELECT c.* FROM contigs c
                    JOIN features f ON c.id = f.contig_id
                    WHERE c.species_id = %s AND (f.gene_name = %s OR f.locus_tag = %s)
                    LIMIT 1
                    """,
                    (species_id, search_term, search_term)
                )
                contig = cursor.fetchone()

            # 5. Nếu không tìm thấy ở bất cứ đâu
            if not contig:
                return jsonify({"error": f"Không tìm thấy Contig/Gen ID '{search_term}' trong loài {species_folder}."}), 404
            
            # 6. Nếu tìm thấy contig, lấy tất cả features của nó
            cursor.execute(
                "SELECT * FROM features WHERE contig_id = %s ORDER BY start_pos",
                (contig['id'],)
            )
            features = cursor.fetchall()

            # 7. Trả về kết quả JSON
            return jsonify({
                "species": species_folder,
                "contig_id": contig['id'],
                "description": contig['description'],
                "sequence_length": contig['sequence_length'],
                # Chuyển list các DictRow features thành list các dict chuẩn
                "features": [dict(f) for f in features] 
            })

    except (Exception, psycopg2.Error) as error:
        print(f"Lỗi truy vấn /api/gene_search: {error}")
        return jsonify({"error": f"Lỗi khi truy vấn CSDL: {error}"}), 500


# -----------------------------------------------------------------
# 4.5. API: LẤY TRÌNH TỰ (Dùng cho gene_view)
# -----------------------------------------------------------------
@app.route('/api/get_sequence', methods=['GET'])
def get_sequence():
    """
    API này chỉ trả về trình tự thô (raw sequence) cho một contig_id.
    Tách riêng để trang gene_view tải metadata trước, sequence sau.
    """
    contig_id = request.args.get('contig_id')
    if not contig_id:
        return jsonify({"error": "Thiếu contig_id"}), 400

    db = get_db()
    if db is None:
        return jsonify({"error": "Lỗi kết nối CSDL."}), 500

    try:
        with db.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # CHỈ chọn cột sequence_data để tiết kiệm bộ nhớ
            cursor.execute("SELECT sequence_data FROM contigs WHERE id = %s", (contig_id,))
            row = cursor.fetchone()
            
            if row:
                # Trả về chuỗi sequence
                return jsonify({"contig_id": contig_id, "sequence": row['sequence_data']})
            else:
                return jsonify({"error": "Không tìm thấy contig"}), 404

    except (Exception, psycopg2.Error) as error:
        print(f"Lỗi truy vấn /api/get_sequence: {error}")
        return jsonify({"error": f"Lỗi khi truy vấn CSDL: {error}"}), 500
    
# -----------------------------------------------------------------
# 5. API: PHÂN TÍCH BLAST (Không thay đổi)
# (Phần này đã dùng file db riêng, không liên quan đến CSDL chính)
# -----------------------------------------------------------------
@app.route('/api/blast/search', methods=['POST'])
def blast_search():
    target_species = request.form.get('species') 
    sequence = request.form['sequence']
    evalue = request.form.get('evalue', '10')
    max_hits = request.form.get('max_hits', '10')

    if not os.path.exists(BLASTN_PATH):
        return jsonify({"error": "BLAST executable không tìm thấy. Vui lòng kiểm tra đường dẫn trong app.py"}), 500
        
    # Check này vẫn đúng vì AVAILABLE_SPECIES là list(SPECIES_MAP.keys())
    if target_species not in AVAILABLE_SPECIES:
        return jsonify({"error": "Loài mục tiêu không hợp lệ."}), 400
        
    # Đường dẫn này vẫn đúng vì nó dựa trên tên thư mục (target_species)
    blast_db = os.path.join("data", target_species, "genome_blast_db") 
    
    if not os.path.exists(blast_db + ".nhr"): # Kiểm tra 1 file CSDL BLAST
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
# 6. API: THIẾT KẾ PRIMER PCR (Không thay đổi)
# (Phần này không liên quan đến CSDL chính)
# -----------------------------------------------------------------

def analyze_primer(seq):
    """Tính toán Tm và %GC cho một trình tự primer."""
    cleaned_seq = re.sub(r'[^ATCGatcg]', '', seq).upper()
    primer = Seq(cleaned_seq)
    
    if len(primer) == 0:
        return {"sequence": "", "length": 0, "tm": "N/A", "gc_percent": "N/A"}

    tm = mt(primer, Na=50.0, c_seq=250.0, strict=False) 
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
    full_sequence_raw = request.form['full_sequence']
    full_sequence = re.sub(r'[^ATCGatcg]', '', full_sequence_raw).upper()

    if len(full_sequence) < 50:
        return jsonify({"error": "Trình tự đầu vào quá ngắn. Cần tối thiểu 50 bp để thiết kế primer."}), 400

    primer_len = 20
    primer_fwd_seq = full_sequence[:primer_len]
    fwd_analysis = analyze_primer(primer_fwd_seq)

    primer_rev_seq = str(Seq(full_sequence[-primer_len:]).reverse_complement())
    rev_analysis = analyze_primer(primer_rev_seq)
    
    product_size = len(full_sequence)

    return jsonify({
        "success": True,
        "product_size": product_size,
        "forward": fwd_analysis,
        "reverse": rev_analysis
    })

# -----------------------------------------------------------------
# 7. JBROWSE ROUTE (Không thay đổi)
# -----------------------------------------------------------------

@app.route('/jbrowse_view')
def jbrowse_view():
    chrom = request.args.get('chrom')
    start = request.args.get('start')
    end = request.args.get('end')
    name = request.args.get('name')
    
    return render_template('jbrowse_view.html', 
                           chrom=chrom, start=start, end=end, name=name)

# -----------------------------------------------------------------
# 8. CHẠY ỨNG DỤNG
# -----------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)