import sys

def main():
    readname_file = sys.argv[1]

    zmw_dic = read_readname_file(readname_file)
    filter_subreads_bam(zmw_dic)

def read_readname_file(readname_file):
    zmw_dic = {}
    with open(readname_file, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            zmw = line.split("/")[0] + "/" + line.split("/")[1]
            zmw_dic[zmw] = 1
    return zmw_dic

def filter_subreads_bam(zmw_dic):
    for line in sys.stdin:
        line = line.rstrip("\n")
        if line.startswith("@"):
            print(line)
        else:
            ll = line.split("\t")
            zmw = ll[0].split("/")[0] + "/" + ll[0].split("/")[1]
            if zmw in zmw_dic:
                print(line)

if __name__ == "__main__":
    main()
