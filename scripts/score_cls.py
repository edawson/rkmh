import sys
from collections import defaultdict

if __name__ == "__main__":
    
    correct_lins = defaultdict(int)
    incorrect_lins = defaultdict(int)
    correct_sublins = defaultdict(int)
    incorrect_sublins = defaultdict(int)

    pcomp = defaultdict(int)
    ptotal = 0
    
    for line in sys.stdin:
        tokens = line.strip().split()

        if tokens[0] == "C":
            tokens[0] = "C1"

        p_lin_match = float(tokens[1].split(":")[1])
        p_sublin_match = float(tokens[2].split(":")[1])

        if tokens[0][0] == tokens[1][0] and p_lin_match > 0.0:
            correct_lins[tokens[0]] += 1
        else:
            incorrect_lins[tokens[0]] += 1

        if tokens[0] == tokens[2].split(":")[0] and p_sublin_match > 0.0:
            correct_sublins[tokens[0]] += 1
        else:
            incorrect_sublins[tokens[0]] += 1

        if p_sublin_match > 0.01:
            pcomp[tokens[2].split(":")[0]] += 1

    for i in pcomp:
        ptotal += pcomp[i]

    print ("Sublins: % correct sublineage: % correct lineage")
    total = 0
    for i in correct_sublins:
        #print (i, correct_sublins[i])
        x = ( float(correct_sublins[i]) / ( float(correct_sublins[i]) + float(incorrect_sublins[i]) ) )
        y = ( float(correct_lins[i]) / ( float(correct_lins[i]) + float(incorrect_lins[i]) ) )
        total += correct_sublins[i]
        total += incorrect_sublins[i]
        print (i, x, y)
    print()

    print ("Estimated sublineage composition")
    for i in pcomp:
       print (i, float(pcomp[i]) / float(ptotal))
