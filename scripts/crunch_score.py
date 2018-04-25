import sys
from collections import defaultdict

if __name__ == "__main__":

# 1020_B2_0-0.0494877072506_D3_0-0.0499643119336_C_0-0.900547980816_4060_.score
    s_d = {}

    err = defaultdict(float)
    est = defaultdict(float)

    strains = sys.argv[1].replace("C", "C1").strip().split("_")[1:-2:2]
    amts = [float(i.split("-")[1]) for i in sys.argv[1].strip().split("_")[2:-2:2]]
    for i in range(0, len(strains)):
        s_d[strains[i]] = amts[i] 
    for i in ["A1", "A2", "A3", "A4", "B1", "B2", "C1", "D1", "D2", "D3"]:
        if i not in s_d:
            s_d[i] = 0.0

    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            if line[0:2] in ["A1", "A2", "A3", "A4", "B1", "B2", "C1", "D1", "D2", "D3"]:
                tokens = line.strip().split()
                
                err[tokens[0]] = abs( float(tokens[1]) - s_d[tokens[0]] )
                est[tokens[0] ]= float(tokens[1])

    
    est_amt_l = []
    est_strain_l = []
    act_strain_l = []
    act_amt_l = []

    num_prim_correct = 0
    num_secondary_correct = 0
    
    total_err = 0.0
    for i in sorted(s_d):
        est_amt_l.append(est[i])
        est_strain_l.append(i)
        act_amt_l.append(s_d[i])
        act_strain_l.append(i)
        if est[i] > 0.005:
            print i, s_d[i], est[i], err[i]
        total_err += err[i]
    print total_err
    est_orders = [i[0] for i in sorted(enumerate(est_amt_l), key=lambda x:x[1], reverse = True)]
    act_orders = [i[0] for i in sorted(enumerate(act_amt_l), key=lambda x:x[1], reverse = True)]
    #print est_orders, act_orders
    s_est = [est_strain_l[i] for i in est_orders]
    s_act = [act_strain_l[i] for i in act_orders]
    if s_est[0] == s_act[0]:
        num_prim_correct += 1
    if s_est[1] == s_act[1]:
        num_secondary_correct += 1
    #print num_prim_correct, num_secondary_correct
    
