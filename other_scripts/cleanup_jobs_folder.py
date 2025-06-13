import subprocess
import os

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Executed: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing {command}: {e}")

if __name__ == "__main__":
    global_path = f'{os.getenv("RUN_PATH")}/jobs/'
    commands = [
        f"rm {global_path}stdall_*",
        f"rm {global_path}htcondor_submission_*.json",
        f"rm -r {global_path}tmp*",
        f"rm {global_path}logs/*"
    ]
    
    for cmd in commands:
        run_command(cmd)