import sys
import argparse
import os
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--strains", type=str, required=False, dest="strains", nargs="+")
    parser.add_argument("-i", "--input", type=str, required=True, dest="searchfiles", nargs="+")
    parser.add_argument("-K", "--kmer", type=int, required=False, default=16, dest="kmer")
    parser.add_argument("-N", "--sketchSize", type=int, required=False, default=1000, dest="sketchSize")
    parser.add_argument("-C", "--coinf", dest="coinf", required=True, type=int, nargs="+")
    parser.add_argument("-M", "--multiclass", dest="multiclass", required=False, action="store_true")
    parser.add_argument("-n", "--normalize", dest="normalize", required=False, action="store_true")
    parser.add_argument("-c", "--collapse", dest="collapse", required=False, action="store_true")
    parser.add_argument("-L", "--label", dest="label", required=False, type=str, default="hpv")
    parser.add_argument("-P", "--predict", dest="predict", required=False, type=str, help="A file 2-column TSV file with mappings from strains names to prediction labels (i.e. integer keys for the wabbit")
    return parser.parse_args()


def make_vw_classes(strainlist):
    vw_dict = {}
    for i in xrange(0, len(strainlist)):
        vw_dict[strainlist[i]] = i
    vw_dict["coinfected"] = len(strainlist)
    return vw_dict

def quantify_strains(strainlist, searchfile, collapse=False):
    
    str_d = Counter()
    if strainlist is not None:
        for i in strainlist:
            str_d[i] = 0

    with open(searchfile, "r") as sfi:
        for line in sfi:
            tokens = line.split("\t")
            #cls = tokens[1].strip().split(" ")[1]
            try:
                if not collapse:
                    cls = tokens[1].strip().split(" ")[1]
                else:
                    cls = tokens[1].strip().split(" ")[1][0]
            except IndexError:
                if not collapse:
                    cls = "unclassified"
                else:
                    cls = "U"
            str_d[cls] += 1;

    return str_d

def vw_line(str_d, isCoinfected, sketchSize, kmer, multiclass, class_d, label_str, normalize, arg_label):
    vw_s = []

    tot = 0
    for i in str_d:
        tot += str_d[i]
    
    if multiclass:
        vw_s.append( class_d[label_str] )
    elif isCoinfected:
        vw_s.append( "1")
    else:
        vw_s.append("-1")
    vw_s.append(" 1.0")
    vw_s.append( " '" + arg_label)
    vw_s.append(" |vir")

    strain_s = ""
    for i in str_d.keys():
        strain_s += " "
        strain_s += str(i)
        strain_s += ":"
        if normalize:
            strain_s += str( float(str_d[i]) / float(tot) )
        else:
            strain_s += str(str_d[i])

    vw_s.append(strain_s)
    #vw_s.append(" |sketch")
    #vw_s.append(" sketchSize=")
    #vw_s.append(sketchSize)
    #vw_s.append(" kmer=")
    #vw_s.append(kmer)

    return "".join([str(i) for i in vw_s])

if __name__ == "__main__":

    args = parse_args()

    class_d = {}
    if args.multiclass:
        class_d["A"] = "1"
        class_d["B"] = "2"
        class_d["C"] = "3"
        class_d["D"] = "4"
        class_d["coinfected"] = "5"
    if args.predict is not None:
        args.multiclass = True
        custom_key_d = {}
        with open(args.predict, "r") as ifi:
            for line in ifi:
                tokens = line.split("\t")
                custom_key_d[tokens[0]] = int(tokens[1])
                class_d = custom_key_d

    for i in xrange(0, len(args.searchfiles)):
        str_d = quantify_strains(args.strains, args.searchfiles[i], args.collapse)
        if args.coinf[i]:
            label_str = "coinfected"
        else:
            label_str = os.path.basename(args.searchfiles[i]).split("_")[0][0]    

        print vw_line(str_d, args.coinf[i], args.sketchSize, args.kmer, args.multiclass, class_d, label_str, args.normalize, args.label)
