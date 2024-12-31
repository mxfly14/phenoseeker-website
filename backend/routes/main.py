from flask import Blueprint, request, jsonify, render_template
from models.distance import find_closest_molecules, find_distance_to_dmso
import logging

# Create Blueprint
bp = Blueprint("main", __name__)

# Configure logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


@bp.route("/", methods=["GET"])
def home():
    """
    Render the welcome page with links.
    """
    return render_template("welcome.html")


@bp.route("/search", methods=["GET"])
def search_page():
    """
    Render the search page for finding closest molecules.
    """
    return render_template("search.html")


@bp.route("/about", methods=["GET"])
def about_page():
    """
    Render the about page with website information.
    """
    return render_template("about.html")


@bp.route("/find_closest", methods=["POST"])
def find_closest():
    """
    Find the closest molecules to a query molecule and calculate the distance to DMSO.

    Expects:
        - query: The molecule ID or InChI to search for.
        - query_type: "id" or "inchi".
        - n: Number of closest molecules to return.

    Returns:
        - results: List of closest molecules.
        - dmso_distance: Distance of the query molecule to DMSO.
    """
    try:
        # Log the incoming request
        logger.info("Received /find_closest request")
        data = request.json
        logger.debug(f"Request data: {data}")

        # Parse request data
        query = data.get("query")
        query_type = data.get("query_type", "id").lower()  # Default to "id"
        n = int(data.get("n", 5))  # Default to 5 closest molecules

        # Validate inputs
        if not query:
            logger.error("Query parameter is missing")
            return jsonify({"error": "Query parameter is required."}), 400
        if query_type not in ["id", "inchi"]:
            logger.error(f"Invalid query type: {query_type}")
            return jsonify({"error": "Query type must be 'id' or 'inchi'."}), 400
        if n <= 0:
            logger.error(f"Invalid number of molecules requested: {n}")
            return (
                jsonify(
                    {"error": "Number of closest molecules (n) must be greater than 0."}
                ),
                400,
            )

        # Log inputs
        logger.info(f"Query: {query}, Query Type: {query_type}, Number of results: {n}")

        # Find closest molecules
        results = find_closest_molecules(query, n, query_by=query_type)
        logger.debug(f"Closest molecules found: {results}")

        if not results:
            logger.warning("No results found for the given query")
            return jsonify({"error": "No results found for the given query."}), 404
        if "error" in results:
            logger.error(f"Error in finding molecules: {results['error']}")
            return jsonify({"error": results["error"]}), 404

        # Calculate the distance to DMSO
        dmso_distance = find_distance_to_dmso(query, query_type)
        logger.debug(f"Distance to DMSO: {dmso_distance}")

        # Return results and distance to DMSO
        response = {"results": results, "dmso_distance": dmso_distance}
        logger.info(f"Response: {response}")
        return jsonify(response)

    except KeyError as e:
        logger.error(f"Missing parameter: {str(e)}")
        return jsonify({"error": f"Missing parameter: {str(e)}"}), 400
    except ValueError as e:
        logger.error(f"Invalid parameter: {str(e)}")
        return jsonify({"error": f"Invalid parameter: {str(e)}"}), 400
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {str(e)}")
        return (
            jsonify({"error": "An unexpected error occurred. Please try again."}),
            500,
        )
