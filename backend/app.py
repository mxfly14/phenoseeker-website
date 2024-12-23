from flask import Flask
from routes import main, health

app = Flask(__name__)

# Register Blueprints
app.register_blueprint(main.bp)
app.register_blueprint(health.bp)

if __name__ == '__main__':
    app.run(debug=True)
