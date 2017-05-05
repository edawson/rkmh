import sys

if __name__ == "__main__":

    str_d = {
        "coinf_correct" : 0,
        "coinf_incorrect" : 0,
        "pure_correct" : 0,
        "pure_incorrect" : 0
            }
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            tokens = line.strip().split(" ")
            if float(tokens[0]) > 0 and tokens[1] == "hpv":
                str_d["pure_correct"] += 1
            elif float(tokens[0]) < 0 and tokens[1] == "hpv":
                str_d["pure_incorrect"] += 1
            elif float(tokens[0]) < 0 and tokens[1] == "coinf":
                str_d["coinf_correct"] += 1
            elif float(tokens[0]) > 0 and tokens[1] == "coinf":
                str_d["coinf_incorrect"] += 1

    print "status", "correct"
    for k in str_d:
        for l in xrange(0, str_d[k]):
            print k.strip().split("_")[0], k.strip().split("_")[1]
