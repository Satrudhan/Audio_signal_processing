from flask import Flask, request, jsonify, render_template

app = Flask(__name__)

# Store received data
received_data = []

@app.route('/receive_data', methods=['POST'])
def receive_data():
    try:
        data = request.get_json()
        if data:
            received_data.append(data)  # Store data
        print("Received Data:", data)
        return jsonify({"message": "Data received"}), 200
    except Exception as e:
        print("Error:", str(e))
        return jsonify({"error": str(e)}), 500

@app.route('/')
def index():
    return render_template("index.html", data=received_data)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
