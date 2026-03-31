from __future__ import annotations

import os
import subprocess
import sys
import time
import webbrowser


def main() -> None:
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    url = "http://127.0.0.1:8501"

    subprocess.Popen(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            "app.py",
            "--server.address",
            "127.0.0.1",
            "--server.port",
            "8501",
            "--server.headless",
            "true",
        ],
        cwd=repo_dir,
    )

    time.sleep(4)
    webbrowser.open(url)
    print(f"Streamlit app should be opening at {url}")
    print("Keep the launched Streamlit window/process running while you use the app.")


if __name__ == "__main__":
    main()
