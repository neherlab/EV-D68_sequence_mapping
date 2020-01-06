import sys

with open(sys.argv[1], "r") as f:
    reference = f.readline().split('\t')[0] #path to the best performing reference

print(reference)
