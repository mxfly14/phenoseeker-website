# PhenoSeeker Website

PhenoSeeker is a tool designed to find molecules phenotypically similar to your positive control. The app is currently in development and runs locally on port 5000.

## Features

- Search for molecules based on their phenotypes.
- Calculate the closest molecules to a given query.
- Display distance to DMSO for comparison.

## Prerequisites

To use this project, you'll need:

- Python 3.10.
- Poetry for dependency management.

## Installation

Follow these steps to set up the project:

1. Clone the repository:
   ```bash
   git clone https://github.com/mxfly14/phenoseeker.git
   cd phenoseeker
   ```

2. Install dependencies using Poetry:
   ```bash
   poetry install
   ```

3. Activate the environment:
   ```bash
   poetry shell
   ```

## Running the Application

1. Start the application:
   ```bash
   python backend/app.py
   ```

2. Open your browser and navigate to:
   ```
   http://localhost:5000
   ```

The app is in development mode and will run locally.

## Directory Structure

- `backend/`: Contains the Flask application.
  - `app.py`: The main entry point for running the app.
  - `routes/`: Directory containing route definitions.
  - `templates/`: HTML templates for rendering pages.

- `README.md`: This documentation file.

## Notes

- This application is under development. Contributions and suggestions are welcome!
- Contact the author via the GitHub repository for further details or issues.

## License

[MIT License](LICENSE)
