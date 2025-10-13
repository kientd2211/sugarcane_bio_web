function searchGene(event) {
    event.preventDefault();
    console.log("Search triggered");
    const contigId = document.getElementById('contig_id').value;
    console.log("Contig ID:", contigId);
    const resultDiv = document.getElementById('result');
    fetch(`/gene/${contigId}`)
        .then(response => response.json())
        .then(data => {
            console.log("Response:", data);
            if (data.error) {
                resultDiv.innerHTML = `<p style="color: red;">${data.error}</p>`;
            } else {
                resultDiv.innerHTML = `
                    <p><strong>Contig ID:</strong> ${data.contig_id}</p>
                    <p><strong>Description:</strong> ${data.description}</p>
                    <p><strong>Sequence Length:</strong> ${data.sequence_length} bp</p>
                    <p><strong>Features:</strong> ${data.features.map(f => `${f[0]}: ${f[1]}`).join(', ')}</p>
                `;
            }
        })
        .catch(error => {
            console.log("Fetch error:", error);
            resultDiv.innerHTML = `<p style="color: red;">Error: ${error}</p>`;
        });
}