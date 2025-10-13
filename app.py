from flask import Flask, render_template, jsonify, request, send_from_directory
from Bio import SeqIO
import os
import subprocess

app = Flask(__name__)

# Load file GBFF khi khởi động ứng dụng
gbff_data = {}
gbff_file = "data/genomic.gbff"
if os.path.exists(gbff_file):
    for record in SeqIO.parse(gbff_file, "genbank"):
        gbff_data[record.id] = {
            "description": record.description,
            "sequence_length": len(record.seq),
            "features": [(f.type, str(f.location)) for f in record.features]
        }
else:
    print(f"Error: File {gbff_file} not found")

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/gene/<name>')
def gene_search(name):
    try:
        if not gbff_data:
            return jsonify({"error": f"File {gbff_file} not found or empty"})
        result = gbff_data.get(name)
        if result:
            return jsonify({
                "contig_id": name,
                "description": result["description"],
                "sequence_length": result["sequence_length"],
                "features": result["features"]
            })
        else:
            return jsonify({"error": f"No contig found with ID {name}"})
    except Exception as e:
        return jsonify({"error": f"Error: {str(e)}"})

@app.route('/blast')
def blast():
    return render_template('blast.html')

@app.route('/blast/search', methods=['POST'])
def blast_search():
    try:
        sequence = request.form['sequence']
        if not sequence or not sequence.strip():
            return jsonify({"error": "Please enter a sequence"})
        
        query_file = "data/query.fasta"
        with open(query_file, "w") as f:
            lines = sequence.strip().split("\n")
            if len(lines) > 1 and lines[0].startswith(">"):
                f.write("\n".join(lines) + "\n")
            else:
                f.write(">query\n" + sequence.strip() + "\n")
        
        if not os.path.exists(query_file):
            return jsonify({"error": "Failed to create query file"})
        
        blastn_path = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastn.exe"
        if not os.path.exists(blastn_path):
            return jsonify({"error": "BLAST executable not found. Check the path in app.py"})
        
        blast_db = "data/genome_blast"
        if not os.path.exists(blast_db + ".nhr"):
            return jsonify({"error": "BLAST database not found. Run makeblastdb first."})
        
        evalue = request.form.get('evalue', '10')
        max_hits = request.form.get('max_hits', '10')
        
        result = subprocess.run(
            [blastn_path, "-db", blast_db, "-query", query_file, "-outfmt", "6", "-evalue", evalue, "-max_target_seqs", max_hits],
            capture_output=True, text=True, check=True
        )
        
        print(f"BLAST output: {result.stdout}")
        return jsonify({"result": result.stdout.splitlines() if result.stdout else ["No hits found"]})
    except subprocess.CalledProcessError as e:
        print(f"BLAST error: {e.stderr}")
        return jsonify({"error": f"BLAST error: {e.stderr}"})
    except Exception as e:
        print(f"General error: {str(e)}")
        return jsonify({"error": f"Error: {str(e)}"})

@app.route('/jbrowse')
def jbrowse():
    return render_template('jbrowse.html')

@app.route('/config.json')
def config_json():
    return send_from_directory(app.static_folder, 'config.json')

if __name__ == '__main__':
    app.run(debug=True)