import argparse
import glob
import os
import subprocess
from sys import argv
from tkinter import Pack

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("../build/RELEASE/")

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iter", type=int, default=64)
parser.add_argument("-o", "--out", default="")
parser.add_argument("--start", default="")

parsed_args, _ = parser.parse_known_args()
args = argv
args[0] = "GenerateFractureData.exe"
args.append("--start")
args.append("0")

files = glob.glob(f"{parsed_args.out}*_stress.csv")
error = 0
while len(files) < parsed_args.iter or error < parsed_args.iter:
    print(f"start: {len(files)}")
    args[-1] = f"{len(files)}"
    subprocess.run(args)
    files = glob.glob(f"{parsed_args.out}*_stress.csv")
    error += 1

