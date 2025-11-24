document.addEventListener('DOMContentLoaded', () => {
    const input = document.getElementById('sequenceInput');
    const btn = document.getElementById('predictBtn');
    const loading = document.getElementById('loading');
    const resultDiv = document.getElementById('result');
    const errorDiv = document.getElementById('error');
    const errorMessage = document.getElementById('errorMessage');
    
    const nextNumberEl = document.getElementById('nextNumber');
    const patternNameEl = document.getElementById('patternName');
    const explanationEl = document.getElementById('explanation');
    const formulaEl = document.getElementById('formula');
    const candidatesListEl = document.getElementById('candidatesList');

    btn.addEventListener('click', handlePredict);
    input.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') handlePredict();
    });

    async function handlePredict() {
        const sequence = input.value.trim();
        if (!sequence) return;

        // Reset UI
        resultDiv.classList.add('hidden');
        errorDiv.classList.add('hidden');
        loading.classList.remove('hidden');
        candidatesListEl.innerHTML = ''; // Clear candidates

        try {
            const response = await fetch('/predict', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ sequence })
            });

            const data = await response.json();

            loading.classList.add('hidden');

            if (response.ok && !data.error) {
                nextNumberEl.textContent = data.formatted_number || data.next_number;
                patternNameEl.textContent = data.pattern;
                explanationEl.textContent = data.explanation || "No explanation provided.";
                formulaEl.textContent = data.formula;
                
                // Render Candidates
                if (data.candidates && data.candidates.length > 0) {
                    data.candidates.forEach(cand => {
                        const div = document.createElement('div');
                        div.className = 'candidate-item';
                        div.innerHTML = `
                            <div class="candidate-header">
                                <span class="candidate-next">${cand.next}</span>
                                <span class="candidate-conf">${Math.round(cand.confidence * 100)}%</span>
                            </div>
                            <div class="candidate-pattern">${cand.pattern}</div>
                            <div class="candidate-expl">${cand.explanation || ''}</div>
                        `;
                        candidatesListEl.appendChild(div);
                    });
                } else {
                    candidatesListEl.innerHTML = '<p>No other candidates found.</p>';
                }

                resultDiv.classList.remove('hidden');
            } else {
                showError(data.error || 'An unknown error occurred.');
            }
        } catch (err) {
            loading.classList.add('hidden');
            showError('Failed to connect to the server.');
            console.error(err);
        }
    }

    function showError(msg) {
        errorMessage.textContent = msg;
        errorDiv.classList.remove('hidden');
    }
});
