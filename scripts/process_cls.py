import sys

if __name__ == "__main__":

    with open(sys.argv[1], "r") as ifi:

        for line in ifi:
            if "rand" in line:
                continue
            tokens = line.strip().split("\t")
            tokens = [i.strip().strip(";") for i in tokens]
            
            print(tokens[0].split("|")[2].split("_")[0], tokens[2].split(";")[0], tokens[3].split(";")[0])
