from flask import Blueprint

bp = Blueprint("health", __name__)


@bp.route("/health", methods=["GET"])
def health_check():
    return {"status": "healthy"}, 200
