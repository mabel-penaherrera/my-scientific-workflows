import json
import subprocess
from concurrent.futures import ProcessPoolExecutor

def run_simulation(entry):
    cmd = [
        "python",
        "run_simulations2.py",
        "--sequence", entry["sequence"],
        "--forcefield", entry["forcefield"],
        "--output", entry["output"],
        "--steps", "1000"
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"Completed: {entry['output']}")
    except subprocess.CalledProcessError as e:
        print(f"Failed: {entry['output']} - {e}")

if __name__ == "__main__":
    with open("batch_config.json", "r") as f:
        config = json.load(f)

    with ProcessPoolExecutor() as executor:
        executor.map(run_simulation, config)
