import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "-infile", dest="infile", type=str, required=True)
    parser.add_argument("-T", "--type", dest="type", type=str, required=False, help="They type of data to interpret")
    return parser.parse_args()

def interpret_binary(x):
    if x > 0:
        print "STATUS:COINF. Sample", tokens[1], "is coinfected. Logit score was", tokens[0]
    elif x == 0:
        print "STATUS:UNKNOWN. Sample", tokens[1], "was inconclusive. Logit score was 0"
    else:
        print "STATUS:PURE. Sample", tokens[1], "is not coinfected. Logit score was", tokens[0]   return False

def interpret_lineage(x):
    if x == 1:
        pass
    elif x == 2:
        pass
    elif x == 3:
        pass
    elif x == 4:
        pass

def interpret_sublineage(x):
    if x == 1:
        pass
    elif x == 2:
        pass
    elif x == 3:
        pass
    elif x == 4:
        pass
    elif x == 5:
        pass
    elif x == 6:
        pass
    elif x == 7:
        pass
    elif x == 8:
        pass
    elif x == 9:
        pass
    return False

if __name__ == "__main__":
    args = parse_args()


    with open(args.infile, "r") as ifi:
        for line in ifi:
            tokens = line.split()
            call = float(tokens[0])
            if args.type == "BINARY":
                interpret_binary(call)
            elif args.type == "SUB":
                interpret_sublineage(x)
            elif args.type == "LIN":
                interpret_lineage(x)
