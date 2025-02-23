import os
import sys
import subprocess

ENV_NAME = "DeepCV"

def run_command(command):
    """Runs a shell command and prints output in real-time."""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    for line in process.stdout:
        print(line, end="")  # Print command output in real-time
    
    process.wait()
    if process.returncode != 0:
        print(f"Error executing command: {command}")
        sys.exit(1)

def setup_conda_environment():
    """Creates a Conda environment and installs dependencies."""
    
    print("\n Checking if Conda environment exists...")
    env_list_cmd = "conda env list"
    env_list = subprocess.run(env_list_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout

    if ENV_NAME in env_list:
        print(f"\n Conda environment '{ENV_NAME}' already exists.")
    else:
        print("\n Creating Conda Environment...")
        run_command(f"conda env create -f environment.yml")

    print("\n Installing Dependencies from environment.yml...")
    run_command(f"conda env update --name {ENV_NAME} --file environment.yml --prune")

    print(f"\n Conda environment '{ENV_NAME}' is ready.")
    print(f"\n Before running the pipeline, activate Conda manually with:")
    print(f"   conda activate {ENV_NAME}\n")

def print_manual_run_instructions(input_vcf, disease_name):
    """Prints manual instructions for running DeepCV after Conda activation."""
    print("\nTo run the DeepCV pipeline, execute the following:")
    command = f"python DeepCV_Main.py {input_vcf}"
    if disease_name:
        command += f' --disease_name "{disease_name}"'
    
    print(f"\n   {command}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python setup_deepcv.py <input_vcf> [--disease_name <name>]")
        sys.exit(1)

    input_vcf = sys.argv[1]
    disease_name = sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == "--disease_name" else None

    setup_conda_environment()
    print_manual_run_instructions(input_vcf, disease_name)
