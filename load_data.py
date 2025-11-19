# load_data.py

import os
import sys
import re
import json
import psycopg2
from psycopg2.extras import execute_values
from Bio import SeqIO
import gzip
import glob

# --- CẤU HÌNH DATABASE ---
DB_HOST = "localhost"
DB_PORT = "5432"
DB_NAME = "gene_database"
DB_USER = "kien"
DB_PASS = "" 
# -------------------------

# --- CẤU HÌNH DỮ LIỆU ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIRECTORY = os.path.join(BASE_DIR, "data")

# Bản đồ tên thư mục sang tên loài (Dùng khi chạy script nạp hàng loạt)
SPECIES_MAP = {
    # 'S_spontaneum': 'Saccharum spontaneum',
    # 'S_officinarum': 'Saccharum officinarum',
    # 'S_bicolor': 'Sorghum bicolor',
    # 'O_sativa': 'Oryza sativa',
    # 'coffee_canephora_cultivar': 'Coffea canephora (Giống)',
    # 'coffee_canephora_strain': 'Coffea canephora (Biến dị)',
    # 'durio_zibethinus': 'Durio zibethinus (Sầu riêng)',
    # 'oryza_sativa_indica': 'Oryza sativa Indica (Lúa)',
    'panax_ginseng': 'Panax ginseng (Nhân sâm)', 
    'selenicereus_undatus': 'Selenicereus undatus (Thanh long)'
}


def connect_db():
    """Kết nối đến cơ sở dữ liệu PostgreSQL."""
    try:
        conn = psycopg2.connect(
            dbname=DB_NAME, user=DB_USER, password=DB_PASS, host=DB_HOST, port=DB_PORT
        )
        print("Kết nối CSDL thành công.")
        return conn
    except psycopg2.OperationalError as e:
        print(f"LỖI: Không thể kết nối tới CSDL '{DB_NAME}'.")
        print(f"Chi tiết lỗi: {e}")
        return None


def parse_location(location_string):
    """Phân tích chuỗi vị trí GenBank."""
    location_string = str(location_string)
    strand = "+"
    if location_string.startswith("complement"):
        strand = "-"
    numbers = [int(n) for n in re.findall(r"\d+", location_string)]
    if not numbers:
        return None, None, strand
    start = min(numbers)
    end = max(numbers)
    return start, end, strand


def find_files(directory):
    """Tìm file .gbff và .fna trong thư mục."""
    gbff_file = None
    fasta_file = None
    
    gb_patterns = ["*.gbff", "*.gb", "*.gbk", "*.gbff.gz", "*.gb.gz", "*.gbk.gz"]
    for ext in gb_patterns:
        found = list(glob.glob(os.path.join(directory, ext)))
        if found:
            gbff_file = found[0]
            break

    fa_patterns = ["*.fasta", ".fa", ".fna", ".fasta.gz", ".fna.gz", ".fa.gz"]
    for ext in fa_patterns:
        found = glob.glob(os.path.join(directory, ext))
        if found:
            fasta_file = found[0]
            break

    return gbff_file, fasta_file


# ==============================================================================
# HÀM XỬ LÝ MỘT LOÀI RIÊNG LẺ (HÀM MỚI BẠN CẦN)
# ==============================================================================
def process_single_species_directory(conn, species_dir, folder_name, common_name):
    print(f"\n--- Đang xử lý loài: {folder_name} ({common_name}) ---")
    
    if not os.path.isdir(species_dir):
        print(f"Bỏ qua: Không tìm thấy thư mục {species_dir}")
        return

    gbff_file_path, fasta_file_path = find_files(species_dir)

    if not fasta_file_path:
        print(f"CẢNH BÁO: Bỏ qua loài {folder_name} vì không tìm thấy file FASTA.")
        return

    try:
        with conn.cursor() as cursor:
            # 1. Đảm bảo loài tồn tại trong bảng 'species' và lấy ID
            cursor.execute(
                """
                INSERT INTO species (name, common_name) VALUES (%s, %s)
                ON CONFLICT (name) DO UPDATE SET common_name = EXCLUDED.common_name
                RETURNING id
                """,
                (folder_name, common_name),
            )
            species_id_row = cursor.fetchone()
            if species_id_row:
                species_id = species_id_row[0]
            else:
                cursor.execute("SELECT id FROM species WHERE name = %s", (folder_name,))
                species_id = cursor.fetchone()[0]
            
            conn.commit() # Commit metadata loài trước

            # 2. Đọc file FASTA
            sequences = {}
            print(f"  Đang đọc trình tự từ: {fasta_file_path}")
            handle_fa = gzip.open(fasta_file_path, "rt") if fasta_file_path.endswith(".gz") else open(fasta_file_path, "r")
            for record in SeqIO.parse(handle_fa, "fasta"):
                sequences[record.id] = str(record.seq)
            handle_fa.close()
            print(f"  Đã tải {len(sequences)} trình tự.")

            # 3. Xử lý GBFF (hoặc chỉ FASTA nếu thiếu GBFF)
            if not gbff_file_path:
                # Trường hợp chỉ có FASTA
                print(f"  Cảnh báo: Không có file GBFF. Chỉ nạp trình tự.")
                contigs_to_insert = []
                for contig_id, seq_data in sequences.items():
                    contigs_to_insert.append((contig_id, species_id, contig_id, len(seq_data), seq_data))
                
                if contigs_to_insert:
                    BATCH_SIZE_FASTA = 10 # An toàn cho FASTA only
                    for i in range(0, len(contigs_to_insert), BATCH_SIZE_FASTA):
                        batch = contigs_to_insert[i:i + BATCH_SIZE_FASTA]
                        psycopg2.extras.execute_values(
                            cursor,
                            """
                            INSERT INTO contigs (id, species_id, description, sequence_length, sequence_data) VALUES %s 
                            ON CONFLICT (id) DO UPDATE SET sequence_data = EXCLUDED.sequence_data
                            """,
                            batch, template="(%s, %s, %s, %s, %s)"
                        )
                conn.commit()
                return

            # Trường hợp có GBFF (Logic Streaming)
            print(f"  Đang đọc và chèn chú giải từ: {gbff_file_path}")
            
            BATCH_SIZE = 1 # GIỮ BATCH_SIZE = 1 ĐỂ TRÁNH SẬP SERVER
            
            contig_batch = []
            feature_batch = []
            records_processed = 0

            handle_gbff = gzip.open(gbff_file_path, "rt") if gbff_file_path.endswith(".gz") else open(gbff_file_path, "r")

            for record in SeqIO.parse(handle_gbff, "genbank"):
                records_processed += 1
                contig_id = record.id
                
                # Lấy trình tự từ map sequences
                sequence_data = sequences.get(contig_id)
                if not sequence_data:
                    base_id = contig_id.split(".")[0]
                    sequence_data = sequences.get(base_id)
                    if not sequence_data and len(record.seq) > 0:
                        sequence_data = str(record.seq)
                
                if not sequence_data:
                    print(f"  LỖI: Contig {contig_id} thiếu trình tự. Bỏ qua.")
                    continue

                # Thêm vào batch
                contig_batch.append((
                    contig_id, species_id, record.description, len(sequence_data), sequence_data
                ))

                # Xử lý Features
                for feature in record.features:
                    if feature.type not in ["gene", "CDS", "mRNA", "tRNA", "rRNA"]: continue
                    
                    loc_str = str(feature.location)
                    start, end, strand = parse_location(loc_str)
                    strand_int = 1 if strand == "+" else -1
                    
                    if start is None: continue
                    
                    # Sửa lỗi '?' ở Gene View: Ép kiểu string cho qualifiers
                    qualifiers = {}
                    for key, value in feature.qualifiers.items():
                        qualifiers[key] = " ".join(value) if isinstance(value, list) else str(value)

                    feature_batch.append((
                        contig_id, feature.type, loc_str, start, end, strand_int,
                        qualifiers.get("gene"), qualifiers.get("locus_tag"), 
                        qualifiers.get("product"), qualifiers.get("translation"),
                        json.dumps(qualifiers)
                    ))

                # Chèn Batch
                if len(contig_batch) >= BATCH_SIZE:
                    print(f"    ... Chèn lô {records_processed} ({len(contig_batch)} contigs)...")
                    psycopg2.extras.execute_values(cursor, 
                        "INSERT INTO contigs (id, species_id, description, sequence_length, sequence_data) VALUES %s ON CONFLICT (id) DO UPDATE SET sequence_data=EXCLUDED.sequence_data", 
                        contig_batch, template="(%s, %s, %s, %s, %s)")
                    
                    if feature_batch:
                        psycopg2.extras.execute_values(cursor, 
                            "INSERT INTO features (contig_id, type, location_string, start_pos, end_pos, strand, gene_name, locus_tag, product_name, translation, qualifiers) VALUES %s ON CONFLICT (id) DO NOTHING", 
                            feature_batch)
                    
                    conn.commit() # Commit ngay lập tức
                    contig_batch.clear()
                    feature_batch.clear()

            # Chèn nốt phần còn lại
            if contig_batch:
                print(f"    ... Chèn lô cuối ...")
                psycopg2.extras.execute_values(cursor, 
                    "INSERT INTO contigs (id, species_id, description, sequence_length, sequence_data) VALUES %s ON CONFLICT (id) DO UPDATE SET sequence_data=EXCLUDED.sequence_data", 
                    contig_batch, template="(%s, %s, %s, %s, %s)")
                if feature_batch:
                    psycopg2.extras.execute_values(cursor, 
                        "INSERT INTO features (contig_id, type, location_string, start_pos, end_pos, strand, gene_name, locus_tag, product_name, translation, qualifiers) VALUES %s ON CONFLICT (id) DO NOTHING", 
                        feature_batch)
                conn.commit()

            handle_gbff.close()
            print(f"  Hoàn tất {records_processed} bản ghi cho {folder_name}.\n")

    except Exception as error:
        print(f"LỖI KHI XỬ LÝ {folder_name}: {error}")
        conn.rollback()
        raise error # Ném lỗi ra để hàm gọi biết


# ==============================================================================
# HÀM MAIN CŨ (GIỜ GỌI HÀM process_single_species_directory)
# ==============================================================================
def process_genome_files(conn):
    """
    Hàm này chạy hàng loạt dựa trên SPECIES_MAP cấu hình ở đầu file.
    """
    print("Bắt đầu quá trình nạp dữ liệu hàng loạt...")
    
    for folder_name, species_name in SPECIES_MAP.items():
        species_dir = os.path.join(DATA_DIRECTORY, folder_name)
        try:
            # GỌI HÀM MỚI TÁCH RA Ở TRÊN
            process_single_species_directory(conn, species_dir, folder_name, species_name)
        except Exception as e:
            print(f" -> Bỏ qua loài {folder_name} do lỗi.")
            continue

    print("\nHoàn tất toàn bộ quá trình.")


# --- Thực thi hàm main khi file này được chạy trực tiếp ---
if __name__ == "__main__":
    db_conn = connect_db()
    if db_conn:
        process_genome_files(db_conn)
        db_conn.close()