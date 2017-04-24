import sys
from collections import Counter
## 5 |strains A1:23146 C:377 B1:546 unclassified:211701 A3:133 A2:212 A4:2230 B2:1052 D2:551 D3:3685 D1:30293 |sketch sketchSize=1000 kmer=16

if __name__ ==  "__main__":
    for line in sys.stdin:
        x_d = Counter()
        tokens = line.split("|")
        features = tokens[1].split(" ")
        for i in features:
            if i.startswith("A"):
                x_d["A"] += int(i.strip().split(":")[1])
            elif i.startswith("B"):
                x_d["B"] += int(i.strip().split(":")[1])
            elif i.startswith("C"):
                x_d["C"] += int(i.strip().split(":")[1])
            elif i.startswith("D"):
                x_d["D"] += int(i.strip().split(":")[1])
            elif i.startswith("u"):
                x_d["U"] = int(i.strip().split(":")[1])
        total = sum([x_d[x] for x in x_d])

        #feat_l = [str(x + ":" + str( float(x_d[x]) / float(total) )) for x in x_d if x is not "U"]
        feat_l = [str(x + ":" + str( float(x_d[x]) / float(total) )) for x in x_d]
        #feat_l = [str(x + ":" + str((x_d[x]) )) for x in x_d]
        x_feat = "|vir " + " ".join(feat_l)
        #xtra_namespace = "|" + tokens[2]
        #print " ".join([ tokens[0].strip(), x_feat, xtra_namespace] ).strip()
        print " ".join([ tokens[0].strip(), x_feat] ).strip()
