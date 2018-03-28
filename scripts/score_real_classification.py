from __future__ import print_function
import sys
from collections import defaultdict

if __name__ == "__main__":

    lin_match_d = defaultdict(int)
    sublin_match_d = defaultdict(int)
    match_threshold = 0.01
    for line in sys.stdin:
        tokens = line.strip().split()
        lin_toks = tokens[3].strip().strip(";").split(";")
        sublin_toks = tokens[4].strip().strip(";").split(";")

        l_trip = False
        s_trip = False
        for i in lin_toks:
            t = i.split(":")
            if (float(t[1]) > match_threshold):
                if l_trip:
                    pass
                    #sys.stderr.write("Read matches to two or more lineages\n" + tokens[0])
                l_trip = True
                lin_match_d[t[0]] += 1

        for s in sublin_toks:
            t = s.split(":")
            if float(t[1]) > match_threshold:
                if s_trip:
                    pass
                    #sys.stderr.write("Read matches to two or more sublineages\n" + tokens[0])
                s_trip = True
                sublin_match_d[t[0]] += 1

    l_total = 0
    s_total = 0

    for i in lin_match_d:
        l_total += lin_match_d[i]
    for i in sublin_match_d:
        s_total += sublin_match_d[i]

    l_pct_d = defaultdict(float)
    s_pct_d = defaultdict(float)
    for i in lin_match_d:
        l_pct_d[i] = float(lin_match_d[i]) / float(l_total)

    for i in sublin_match_d:
        s_pct_d[i] = float(sublin_match_d[i]) / float(s_total)

    print(l_pct_d)
    print(s_pct_d)


