from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')

from Bio import SeqIO
@app.route('/gene/<name>')
def gene_search(name):
    # Giả lập: Tìm gene trong file GFF (cần file thực tế)
    return {"gene": name, "result": f"Found {name} in genome"}

if __name__ == '__main__':
    app.run(debug=True)