from __future__ import print_function
import sys
from collections import defaultdict
import pprint

def dict_to_string(d):
    t = []
    for i in d:
        x = ":".join([str(i), str(d[i])])
        t.append(x)
    t = sorted(t, reverse = True, key = lambda x : float(x.split(":")[1]))
    return ";".join(t)

if __name__ == "__main__":

    lin_match_d = defaultdict(int)
    sublin_match_d = defaultdict(int)
    match_threshold = 0.005
    for line in sys.stdin:
        tokens = line.strip().split()
        read_len = int(tokens[2].strip().split("/")[1])
        hpv_match = int(tokens[2].strip().split("/")[0])
        
        ## We get some reads that look like ION Torrent barf - toss those:
        if read_len < 50 or hpv_match < 15:
            continue

        lin_toks = tokens[3].strip().strip(";").split(";")
        lin_kmer_counts = [int(i) for i in tokens[5].strip().strip(";").split(";")]
        sublin_toks = tokens[4].strip().strip(";").split(";")
        sublin_kmer_counts = [int(i) for i in tokens[6].strip().strip(";").split(";")]

        l_trip = False
        s_trip = False
        l_match = ""
        for i in range(0, len(lin_toks)):
            t = lin_toks[i].split(":")
            #if (float(t[1]) > match_threshold) and lin_kmer_counts[i] > 4:
            if lin_kmer_counts[i] > 5:
                if l_trip:
                    #l_match = ""
                    #sys.stderr.write("Read matches to two or more lineages\n" + tokens[0])
                    break
                else:
                    l_trip = True
                    l_match = t[0]
       
        s_match = ""
        for i in range(0, len(sublin_toks)):
            t = sublin_toks[i].split(":")
            if sublin_kmer_counts[i] > 2 and float(t[1]) > match_threshold:
            #if float(t[1]) > match_threshold:
                if s_trip:
                    #s_match = ""
                    #sys.stderr.write("Read matches to two or more sublineages" + tokens[0] + "\n")
                    break
                s_trip = True
                s_match = t[0]
        
        if l_match is not "" and s_match is not "" and l_match is not s_match[0]:
            old = ""
            if lin_kmer_counts[0] > 10 and sublin_kmer_counts[1] > 2 and lin_toks[0].split(":")[0] == sublin_toks[1].split(":")[0][0]:
                old = s_match 
                s_match = sublin_toks[1].split(":")[0]
                sys.stderr.write("Lin / Sublin mistmatch: " + l_match + " " + old)
                sys.stderr.write( " " + old + "->" + s_match + "\n")
            else:
                s_match = ""
        
        if l_match is not "":
            lin_match_d[l_match] += 1
 
        if s_match is not "":
            sublin_match_d[s_match] += 1

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

    low_read_lins = ""
    if l_total < 1000:
        low_read_lins = "WARN:low_lineage_counts:" + str(l_total)
    else:
        low_read_lins = "INFO:low_lineage_counts:" + str(l_total)
    low_read_sublins = ""
    if s_total < 1000:
        low_read_sublins = "WARN:low_sublineage_counts:" + str(s_total)
    else:
        low_read_sublins = "INFO:low_sublineage_counts:" + str(s_total)
    #pprint.pprint(l_pct_d)
    #pprint.pprint(s_pct_d)
    #pprint.pprint (sublin_match_d)
    print( dict_to_string(l_pct_d), dict_to_string(s_pct_d), dict_to_string(sublin_match_d), low_read_lins, low_read_sublins)

