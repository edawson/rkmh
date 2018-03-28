from __future__ import print_function
import sys
from collections import defaultdict

if __name__ == "__main__":
    
    s_d = defaultdict(float)
    strains = sys.argv[1].replace("C", "C1").strip().split("_")[1:-2:2]
    amts = [float(i.split("-")[1]) for i in sys.argv[1].strip().split("_")[2:-2:2]]
    for i in range(0, len(strains)):
        s_d[strains[i]] = amts[i]
    for i in ["A1", "A2", "A3", "A4", "B1", "B2", "C1", "D1", "D2", "D3"]:
        if i not in s_d:
            s_d[i] = 0.0
    #print (s_d)

    amt_indices = sorted(range(len(amts)), key=lambda k: amts[k], reverse = True)

    lin_match_d = defaultdict(int)
    sublin_match_d = defaultdict(int)
    match_threshold = 0.01
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
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
    
    #print (s_d)
    #print(l_pct_d)
    #print(s_pct_d)
    s_pct_l = []
    s_strn_l = []
    for i in s_pct_d:
        s_pct_l.append(s_pct_d[i])
        s_strn_l.append(i)
    total_err = 0.0
    for i in s_d:
        err = abs(s_pct_d[i] - s_d[i])
        total_err += err

    s_p_indices = sorted(range(len(s_pct_l)), key=lambda k: s_pct_l[k], reverse = True)
    correct_primary = strains[amt_indices[0]] == s_strn_l[s_p_indices[0]]
    correct_secondary = strains[amt_indices[1]] == s_strn_l[s_p_indices[1]]
    diff = abs(amts[amt_indices[0]] - amts[amt_indices[1]])
    sdiff = 0.0
    if len(amts) >=3:
        sdiff = abs(amts[amt_indices[1]] - amts[amt_indices[2]])
    flipped = strains[amt_indices[0]] == s_strn_l[s_p_indices[1]] and strains[amt_indices[1]] == s_strn_l[s_p_indices[0]]
    primary_detected = s_pct_d[strains[amt_indices[0]]] > 0.05 and abs(s_pct_d[strains[amt_indices[0]]] - s_d[strains[amt_indices[0]]]) < 0.1
    secondary_detected =  abs(s_pct_d[strains[amt_indices[1]]] - s_d[strains[amt_indices[1]]]) < 0.1
    print (correct_primary, correct_secondary, flipped, \
            primary_detected, secondary_detected, len([i for i in s_d if s_d[i] > 0.005]), diff, sdiff, total_err, \
            s_d, s_pct_d, sys.argv[1])
            #":".join([str(amts[amt_indices[0]]), str(s_pct_l[s_p_indices[0]])]), \
            #":".join([str(amts[amt_indices[1]]), str(s_pct_l[s_p_indices[1]])]), \
            #":".join([str(strains[amt_indices[0]]), str(s_strn_l[s_p_indices[0]])]), \
            #":".join([str(strains[amt_indices[1]]), str(s_strn_l[s_p_indices[1]])]))
 
    

