#!/usr/bin/env python3
# helper function to add gene names to the links file
import sys

next(sys.stdin)
print("chr1\tgene1\tpos1\tchr2\tgene2\tpos2")

for line in sys.stdin:
    chr1, pos1, chr2, pos2 = line.rstrip().split("\t")
    print(chr1, f"{chr1}:{pos1}", pos1, chr2, f"{chr2}:{pos2}", pos2, sep="\t")
