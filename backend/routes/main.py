from flask import Blueprint, request, jsonify
from models.distance import find_closest_molecules

bp = Blueprint('main', __name__)

@bp.route('/find_closest', methods=['POST'])
def find_closest():
    try:
        data = request.json
        pos_control_id = data.get('positive_control_id')
        n = int(data.get('n', 5))  # Default to top 5 closest molecules
        
        if not pos_control_id or n <= 0:
            return jsonify({'error': 'Invalid input parameters'}), 400

        results = find_closest_molecules(pos_control_id, n)
        if 'error' in results:
            return jsonify(results), 404

        return jsonify(results)
    except Exception as e:
        return jsonify({'error': str(e)}), 500
