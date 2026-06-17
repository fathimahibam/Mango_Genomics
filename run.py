import os
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Only insert local precompiled libs under Windows with Python 3.14
if sys.platform == 'win32' and sys.version_info[:2] == (3, 14):
    libs_path = os.path.join(BASE_DIR, 'libs')
    if os.path.exists(libs_path) and libs_path not in sys.path:
        sys.path.insert(0, libs_path)

from app import app
from waitress import serve

if __name__ == '__main__':
    # Listen on PORT environment variable (for Render/Railway) or default 5000
    port = int(os.environ.get('PORT', 5000))
    print(f"Starting Waitress production server on http://0.0.0.0:{port} ...")
    serve(app, host='0.0.0.0', port=port, threads=8)
