from flask import Flask, render_template, request, jsonify
import subprocess
import os
import json

app = Flask(__name__)

# Path to the CLI executable
CLI_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'seqfinder_cli.exe'))

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    if not data or 'sequence' not in data:
        return jsonify({'error': 'No sequence provided'}), 400

    sequence_str = data['sequence']
    # Split by comma or space
    try:
        # Normalize input: replace commas with spaces, then split
        # Allow characters like / and ^ for fractions and powers
        # We don't convert to float here, we let the CLI handle parsing
        parts = sequence_str.replace(',', ' ').split()
        
        # Basic validation: check if parts look like numbers or valid expressions
        import re
        valid_pattern = re.compile(r'^[0-9./\^+\-]+$')
        for p in parts:
            if not valid_pattern.match(p):
                 return jsonify({'error': 'Invalid characters in sequence'}), 400
                 
        numbers = parts # Pass raw strings
    except Exception:
        return jsonify({'error': 'Invalid sequence format'}), 400

    if not numbers:
        return jsonify({'error': 'Empty sequence'}), 400

    # Call the CLI
    try:
        result = subprocess.run([CLI_PATH] + numbers, capture_output=True, text=True, timeout=5)
        output = result.stdout.strip()
        
        if not output:
             return jsonify({'error': 'No output from solver'}), 500
             
        try:
            json_output = json.loads(output)
            return jsonify(json_output)
        except json.JSONDecodeError:
             return jsonify({'error': 'Invalid response from solver', 'raw': output}), 500

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)
