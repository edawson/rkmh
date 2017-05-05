import random
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="infile", required=True, help="A two column file in format LABEL\tPATH_TO_FILE with the strain label and the path to a file containing reads")
    parser.add_argument("-n", "--num-infs", dest="num", required=False, default=10, help="The number of simulated infections to generate", type=int)
    parser.add_argument("-c", "--coinfection", dest="coinfection", action="store_true", help="Generate coinfections (default is to generate pure samples).")
    parser.add_argument("-m", "--min-cov", dest="mincov", type=int, help="Minimum coverage for a strain", default=1000)
    return parser.parse_args()


def randproportions(numprops):

    ret_l = []

    if numprops == 1.0:
        return [1.0]
    init = random.uniform(0.04, 0.50)
    ret_l.append(init)
    for i in xrange(1, numprops - 1):
        rand_prop_in_interval = random.uniform(0.04, init)
        ret_l.append( rand_prop_in_interval )
        init = rand_prop_in_interval
    ret_l.append(1.0 - sum(ret_l))
    
    for i in ret_l:
        if i < 0:
            return []
    return ret_l

def help_negs(strain_d, mincov):
    ran_strain_list = []
    mixstr = ""
    while True:
        mixstr = ""
        ## Pick a random coverage > mincov
        ran_cov = random.randint(mincov, 10000)
        ## Pick a number of strains
        ran_num_strains = random.randint(2, len(strain_d))
        ## pick strains
        ran_strain_list = random.sample(strain_d.keys(), ran_num_strains)
        ## do some math - figure out the proportions for each strain
        prop_l = randproportions(ran_num_strains)
        if len(prop_l) == 0:
            continue

        ## create output line
        for i in xrange(0, ran_num_strains):
            line_str = "\t".join( [ ran_strain_list[i], str(prop_l[i]), str(ran_cov), strain_d[ ran_strain_list[i] ], "\n"] )
            mixstr = mixstr + line_str
        break

    return mixstr


def randmix(strain_d, coinfection, mincov):
    hf = help_negs
    if coinfection:
        mixstr = hf(strain_d, mincov)
    else:
    ## Pick a random strain and a random coverage > mincov
        ran_cov = random.randint(mincov, 10000)
        ran_strain = random.sample(strain_d.keys(), 1)[0]
        strain_fi = strain_d[ ran_strain ]
        mixstr = "\t".join( [ran_strain,  "1.0", str(ran_cov), strain_fi, "\n"] )


    return mixstr

if __name__ == "__main__":
    args = parse_args()
    
    strain_d = {}

    with open(args.infile) as ifi:
        for line in ifi:
            tokens = line.strip().split("\t")
            strain_d[ tokens[0] ] = tokens[1]

    for i in xrange(0, args.num):
        print randmix(strain_d, args.coinfection, args.mincov)
