import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description='run the tests.')
parser.add_argument('path', metavar='in', type=str, help='path to file')
args = parser.parse_args()
py = sys.executable

print("\n\n")
print(py)
print("/*******-RUN TEST: " + args.path + "-*******/")
subprocess.call(py +" "+ args.path, shell=True)