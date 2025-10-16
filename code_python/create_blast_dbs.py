import subprocess
import os

# ---------------------------------------------------------------------
# BƯỚC 1: CẬP NHẬT ĐƯỜNG DẪN BLAST+ TRÊN MÁY BẠN
# ---------------------------------------------------------------------
# Ví dụ:
BLAST_BIN_PATH = r"C:\Program Files\NCBI\blast-2.17.0+\bin"

# Tên file thực thi
MAKEBLASTDB = os.path.join(BLAST_BIN_PATH, "makeblastdb.exe")

# Kiểm tra sự tồn tại của công cụ
if not os.path.exists(MAKEBLASTDB):
    print(f"LỖI CẤU HÌNH: Không tìm thấy makeblastdb.exe tại {MAKEBLASTDB}.")
    print("Vui lòng kiểm tra lại đường dẫn và cài đặt NCBI BLAST+.")
    exit()

# ---------------------------------------------------------------------
# BƯỚC 2: KHAI BÁO DỮ LIỆU ĐẦU VÀO VÀ TÊN DB ĐẦU RA
# ---------------------------------------------------------------------
species_list = ["S_spontaneum", "S_officinarum", "S_bicolor", "O_sativa"]
data_dir = "data"

print("--- BẮT ĐẦU TẠO BLAST DATABASE ---")

for species in species_list:
    # 1. Định vị file FASTA đầu vào
    fasta_file = os.path.join(data_dir, species, "genome.fasta")
    
    # 2. Đặt tên file DB đầu ra (sẽ được lưu trong thư mục của loài)
    db_out = os.path.join(data_dir, species, "genome_blast_db") 
    
    if not os.path.exists(fasta_file):
        print(f"CẢNH BÁO: Bỏ qua {species}. Không tìm thấy file {fasta_file}.")
        continue

    print(f"\nĐang tạo BLAST DB cho loài: {species}...")
    
    # -----------------------------------------------------------------
    # BƯỚC 3: CHẠY LỆNH makeblastdb BẰNG subprocess
    # -----------------------------------------------------------------
    try:
        command = [
            MAKEBLASTDB, 
            "-in", fasta_file,      # File FASTA đầu vào
            "-dbtype", "nucl",      # Loại DB: nucleotide (nucl)
            "-out", db_out,         # Tên cơ sở dữ liệu đầu ra
            "-title", f"{species} Genome DB" # Tiêu đề DB (tùy chọn)
        ]
        
        # Chạy lệnh
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        
        print(f"   -> Thành công. DB được lưu dưới tên: {db_out}")
        # print(result.stdout) # Có thể bỏ chú thích để xem output chi tiết của BLAST
        
    except subprocess.CalledProcessError as e:
        # Xử lý lỗi nếu makeblastdb thất bại
        print(f"   !!! LỖI khi tạo BLAST DB cho {species} !!!")
        print(f"   Error: {e.stderr}")
    except Exception as e:
        print(f"   Lỗi chung: {e}")

print("\n--- HOÀN THÀNH TẠO BLAST DATABASE CHO TẤT CẢ CÁC LOÀI ---")