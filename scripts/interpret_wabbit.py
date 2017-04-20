import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "-infile", dest="infile", type=str, required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()


    with open(args.infile, "r") as ifi:
        for line in ifi:
            tokens = line.split()
            call = float(tokens[0])
            if call < 0:
                print "STATUS:COINF. Sample", tokens[1], "is coinfected. Logit score was", tokens[0]
            elif call == 0:
                print "STATUS:UNKNOWN. Sample", tokens[1], "was inconclusive. Logit score was 0"
            else:
                print "STATUS:PURE. Sample", tokens[1], "is not coinfected. Logit score was", tokens[0]
