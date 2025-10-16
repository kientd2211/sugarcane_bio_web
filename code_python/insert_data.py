from Bio import SeqIO
import psycopg2

# Thông tin kết nối PostgreSQL
conn_params = {
    "dbname": "research_genome_db",
    "user": "postgres",
    "password": "ekgis@2907",
    "host": "124.158.4.142",
    "port": "5434"
}

conn = psycopg2.connect(**conn_params)
c = conn.cursor()

# Tạo bảng nếu chưa có
c.execute("""
    CREATE TABLE IF NOT EXISTS GenomeAssembly (
        assembly_id SERIAL PRIMARY KEY,
        variety_name VARCHAR(255),
        release_date TIMESTAMP,
        total_length DOUBLE PRECISION
    );
    CREATE TABLE IF NOT EXISTS Chromosome (
        chr_id SERIAL PRIMARY KEY,
        assembly_id INTEGER REFERENCES GenomeAssembly(assembly_id),
        chr_name VARCHAR(255),
        chr_length DOUBLE PRECISION
    );
    CREATE TABLE IF NOT EXISTS Gene (
        gene_id SERIAL PRIMARY KEY,
        chr_id INTEGER REFERENCES Chromosome(chr_id),
        gene_start_pos VARCHAR(50),
        gene_end_pos VARCHAR(50),
        name VARCHAR(255),
        description VARCHAR(255)
    );
    CREATE TABLE IF NOT EXISTS Annotation (
        annot_id SERIAL PRIMARY KEY,
        gene_id INTEGER REFERENCES Gene(gene_id),
        type VARCHAR(255),
        sequence VARCHAR(255)
    );
""")
conn.commit()

# Chèn vào GenomeAssembly và lấy assembly_id
c.execute("INSERT INTO GenomeAssembly (variety_name, release_date, total_length) VALUES (%s, %s, %s) RETURNING assembly_id", 
          ('Saccharum spontaneum', '2023-01-01', 0.0))
assembly_id = c.fetchone()[0]
conn.commit()

total_length = 0.0
for record in SeqIO.parse("data/genomic.gbff", "genbank"):
    chr_length = len(record.seq)
    total_length += chr_length
    c.execute("INSERT INTO Chromosome (assembly_id, chr_name, chr_length) VALUES (%s, %s, %s) RETURNING chr_id", 
              (assembly_id, record.id, chr_length))
    chr_id = c.fetchone()[0]
    conn.commit()  # Commit sau mỗi Chromosome
    
    for feature in record.features:
        if feature.type == 'gene' or feature.type == 'source':
            start = str(feature.location.start)
            end = str(feature.location.end)
            name = feature.qualifiers.get('gene', ['Unknown'])[0]
            description = feature.qualifiers.get('note', ['N/A'])[0]
            c.execute("INSERT INTO Gene (chr_id, gene_start_pos, gene_end_pos, name, description) VALUES (%s, %s, %s, %s, %s) RETURNING gene_id", 
                      (chr_id, start, end, name, description))
            gene_id = c.fetchone()[0]  # Lấy gene_id thực tế
            
            annot_id = 1
            annot_type = feature.type
            sequence = str(record.seq[feature.location.start:feature.location.end])[:255]  # Cắt ngắn nếu quá 255
            c.execute("INSERT INTO Annotation (gene_id, type, sequence) VALUES (%s, %s, %s) RETURNING annot_id", 
                      (gene_id, annot_type, sequence))
            annot_id = c.fetchone()[0]  # Lấy annot_id thực tế
            conn.commit()  # Commit sau mỗi Annotation

# Cập nhật total_length và commit cuối cùng
c.execute("UPDATE GenomeAssembly SET total_length = %s WHERE assembly_id = %s", (total_length, assembly_id))
conn.commit()

conn.close()
print("Dữ liệu đã được chèn thành công vào PostgreSQL.")