import os
import subprocess

# Initialize VOMS proxy
print("Initializing VOMS proxy...")
subprocess.run("voms-proxy-init --rfc --voms cms -valid 192:00", shell=True, check=True)
outbase_dir = "/eos/user/" + os.environ['USER'][:1] + "/" + os.environ['USER'] + "/HiggsDNA_v4PrelimProd/"
extra_dir = ""

samples = [
    # 2022preEE
    {
        "keyword": "Run2022C",
        "cmsdas": "/EGamma/Run2022C-22Sep2023-v1/NANOAOD",
        "year": "2022preEE",
        "nano": "12"
    },
    {
        "keyword": "Run2022D",
        "cmsdas": "/EGamma/Run2022D-22Sep2023-v1/NANOAOD",
        "year": "2022preEE",
        "nano": "12"
    },

    # 2022postEE
    {
        "keyword": "Run2022E",
        "cmsdas": "/EGamma/Run2022E-22Sep2023-v1/NANOAOD",
        "year": "2022postEE",
        "nano": "12"
    },
    {
        "keyword": "Run2022F",
        "cmsdas": "/EGamma/Run2022F-22Sep2023-v1/NANOAOD",
        "year": "2022postEE",
        "nano": "12"
    },
    {
        "keyword": "Run2022G",
        "cmsdas": "/EGamma/Run2022G-22Sep2023-v2/NANOAOD",
        "year": "2022postEE",
        "nano": "12"
    },

    # 2023preBPix
    {
        "keyword": "Run2023Cv1_EG0",
        "cmsdas": "/EGamma0/Run2023C-22Sep2023_v1-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv2_EG0",
        "cmsdas": "/EGamma0/Run2023C-22Sep2023_v2-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv3_EG0",
        "cmsdas": "/EGamma0/Run2023C-22Sep2023_v3-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv4_EG0",
        "cmsdas": "/EGamma0/Run2023C-22Sep2023_v4-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv1_EG1",
        "cmsdas": "/EGamma1/Run2023C-22Sep2023_v1-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv2_EG1",
        "cmsdas": "/EGamma1/Run2023C-22Sep2023_v2-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv3_EG1",
        "cmsdas": "/EGamma1/Run2023C-22Sep2023_v3-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Cv4_EG1",
        "cmsdas": "/EGamma1/Run2023C-22Sep2023_v4-v1/NANOAOD",
        "year": "2023preBPix",
        "nano": "12"
    },

    # 2023postBPix
    {
        "keyword": "Run2023Dv1_EG0",
        "cmsdas": "/EGamma0/Run2023D-22Sep2023_v1-v1/NANOAOD",
        "year": "2023postBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Dv2_EG0",
        "cmsdas": "/EGamma0/Run2023D-22Sep2023_v2-v1/NANOAOD",
        "year": "2023postBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Dv1_EG1",
        "cmsdas": "/EGamma1/Run2023D-22Sep2023_v1-v1/NANOAOD",
        "year": "2023postBPix",
        "nano": "12"
    },
    {
        "keyword": "Run2023Dv2_EG1",
        "cmsdas": "/EGamma1/Run2023D-22Sep2023_v2-v1/NANOAOD",
        "year": "2023postBPix",
        "nano": "12"
    },

    # 2024
    {
        "keyword": "Run2024C_EG0",
        "cmsdas": "/EGamma0/Run2024C-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024C_EG1",
        "cmsdas": "/EGamma1/Run2024C-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024D_EG0",
        "cmsdas": "/EGamma0/Run2024D-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024D_EG1",
        "cmsdas": "/EGamma1/Run2024D-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024E_EG0",
        "cmsdas": "/EGamma0/Run2024E-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024E_EG1",
        "cmsdas": "/EGamma1/Run2024E-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024F_EG0",
        "cmsdas": "/EGamma0/Run2024F-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024F_EG1",
        "cmsdas": "/EGamma1/Run2024F-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024G_EG0",
        "cmsdas": "/EGamma0/Run2024G-MINIv6NANOv15-v2/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024G_EG1",
        "cmsdas": "/EGamma1/Run2024G-MINIv6NANOv15-v2/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024H_EG0",
        "cmsdas": "/EGamma0/Run2024H-MINIv6NANOv15-v2/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024H_EG1",
        "cmsdas": "/EGamma1/Run2024H-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024Iv1_EG0",
        "cmsdas": "/EGamma0/Run2024I-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024Iv2_EG0",
        "cmsdas": "/EGamma0/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024Iv1_EG1",
        "cmsdas": "/EGamma1/Run2024I-MINIv6NANOv15-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },
    {
        "keyword": "Run2024Iv2_EG1",
        "cmsdas": "/EGamma1/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD",
        "year": "2024",
        "nano": "15"
    },

    # 2025
    {
        "keyword": "Run2025Cv1_EG0",
        "cmsdas": "/EGamma0/Run2025C-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv2_EG0",
        "cmsdas": "/EGamma0/Run2025C-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv1_EG1",
        "cmsdas": "/EGamma1/Run2025C-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv2_EG1",
        "cmsdas": "/EGamma1/Run2025C-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv1_EG2",
        "cmsdas": "/EGamma2/Run2025C-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv2_EG2",
        "cmsdas": "/EGamma2/Run2025C-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv1_EG3",
        "cmsdas": "/EGamma3/Run2025C-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Cv2_EG3",
        "cmsdas": "/EGamma3/Run2025C-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025D_EG0",
        "cmsdas": "/EGamma0/Run2025D-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025D_EG1",
        "cmsdas": "/EGamma1/Run2025D-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025D_EG2",
        "cmsdas": "/EGamma2/Run2025D-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025D_EG3",
        "cmsdas": "/EGamma3/Run2025D-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025E_EG0",
        "cmsdas": "/EGamma0/Run2025E-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025E_EG1",
        "cmsdas": "/EGamma1/Run2025E-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025E_EG2",
        "cmsdas": "/EGamma2/Run2025E-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025E_EG3",
        "cmsdas": "/EGamma3/Run2025E-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv1_EG0",
        "cmsdas": "/EGamma0/Run2025F-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv2_EG0",
        "cmsdas": "/EGamma0/Run2025F-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv1_EG1",
        "cmsdas": "/EGamma1/Run2025F-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv2_EG1",
        "cmsdas": "/EGamma1/Run2025F-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv1_EG2",
        "cmsdas": "/EGamma2/Run2025F-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv2_EG2",
        "cmsdas": "/EGamma2/Run2025F-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv1_EG3",
        "cmsdas": "/EGamma3/Run2025F-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025Fv2_EG3",
        "cmsdas": "/EGamma3/Run2025F-PromptReco-v2/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025G_EG0",
        "cmsdas": "/EGamma0/Run2025G-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025G_EG1",
        "cmsdas": "/EGamma1/Run2025G-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025G_EG2",
        "cmsdas": "/EGamma2/Run2025G-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    },
    {
        "keyword": "Run2025G_EG3",
        "cmsdas": "/EGamma3/Run2025G-PromptReco-v1/NANOAOD",
        "year": "2025",
        "nano": "15"
    }
]

# Iterate over each sample and execute the command
for sample in samples:
    parent_dir = outbase_dir + "/" + sample['year'] + "/" + extra_dir

    # Create the directory
    print(f"Creating directory: {parent_dir}")
    os.makedirs(parent_dir, exist_ok=True)
    os.chmod(parent_dir, 0o777)

    # Construct the command
    command = f"python submission/tools_HHbbgg/produce_one_data.py --keyword {sample['keyword']} --cmsdas {sample['cmsdas']} --parent-dir {parent_dir} --year {sample['year']} --nano {sample['nano']}"  # --memory 20GB"

    print(f"Executing: {command}")  # Print the command being executed
    subprocess.run(command, shell=True, check=True)

print("All jobs have been executed successfully.")
