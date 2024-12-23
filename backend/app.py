from flask import Flask
from routes import main, health

# Initialize Flask app
app = Flask(__name__)

# Register Blueprints
app.register_blueprint(main.bp)
app.register_blueprint(health.bp)

if __name__ == '__main__':
    # Run the app in debug mode
    app.run(debug=True, host="0.0.0.0")
