import sys
import json
import argparse
from collections import Counter

## 1 1.0 'hpv |vir C:0.00756124087213 A1:0.0226319332954 A3:0.922005282511 A2:0.00564503599358 A4:0.00548966803045 B2:0.0019679941996 D1:0.0019679941996 D2:0.00424672432544 D3:0.0129991195815 B1:0.0154850069916
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", dest="infile", type=str)
    parser.add_argument("-c", "--coinf", dest="coinf", action="store_true")
    parser.add_argument("-k", "--kollapse", dest="kollapse", action="store_true")

    return parser.parse_args()

def collapse(sketch):
    c = Counter()
    for i in sketch:
        c[i] += 1
    
    return " ".join( [ (str(i) + ":" + str(c[i])) for i in c  ] )

def j_to_v(j, isCoinf=False, kollapse=False):
    vec = "1" if isCoinf else "0"
    vec += " 1.0 "
    vec += ("`" + ("_".join(j["name"].split("|"))) )
    vec += " |f "
    if not kollapse:
        vec += " ".join( [ (str(i) + ":1") for i in j["sketches"] ] )
    else:
        vec += collapse(j["sketches"])

    return vec


if __name__ == "__main__":
   
    args = parse_args()
    x = []
    with open(args.infile, "r") as ifi:
        x = json.load(ifi)
    for i in x:
        print j_to_v(i, args.coinf, args.kollapse)
