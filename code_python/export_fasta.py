# export_fasta.py
import psycopg2

# Thông tin kết nối PostgreSQL
conn_params = {
    "dbname": "research_genome_db",
    "user": "postgres",
    "password": "ekgis@2907",
    "host": "124.158.4.142",
    "port": "5434"
}

# Kết nối
conn = psycopg2.connect(**conn_params)
c = conn.cursor()

# Xuất dữ liệu từ Chromosome
output_file = r"D:\biowork\sugarcane_bio_web\data_exported\sugarcane.fasta"
with open(output_file, 'w') as f:
    c.execute("SELECT chr_name, chr_length FROM Chromosome")
    chromosomes = c.fetchall()
    for chr_name, chr_length in chromosomes:
        # Giả định chuỗi đầy đủ không có trong CSDL, cần từ file genomic.fasta
        # Nếu có sequence trong Annotation, dùng cách khác (xem dưới)
        f.write(f">{chr_name}\n")
        # Thêm chuỗi giả (nếu không có, cần ghép từ genomic.fasta hoặc Annotation)
        f.write("N" * int(chr_length) + "\n")  # Placeholder, thay bằng chuỗi thực tế

conn.close()
print(f"File FASTA đã được xuất tại: {output_file}")