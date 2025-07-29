import sys

def main():
    deepconsensus_reads_ls = sys.argv[1:] # /rd/caiya/Tn5_duplex/C01.deepconsensus.cx3.fwd.fastq

    deepconsensus_dic = read_deepconsensus_reads(deepconsensus_reads_ls)
    read_input_ccs_bam(deepconsensus_dic)


def read_deepconsensus_reads(deepconsensus_reads_ls):
    deepconsensus_dic = {}
    for deepconsensus_reads in deepconsensus_reads_ls:
        c = 0
        with open(deepconsensus_reads, "r") as inf:
            for line in inf:
                c += 1
                line = line.rstrip("\n")
                if c % 4 == 1:
                    zmw = line.split("/")[1]
                    strand = line.split("/")[3]
                    deepconsensus_dic.setdefault(zmw, {})[strand] = ""
                elif c % 4 == 2:
                    deepconsensus_dic[zmw][strand] = line
                elif c % 4 == 0:
                    deepconsensus_dic[zmw][strand] = deepconsensus_dic[zmw][strand] + "\t" + line
    return deepconsensus_dic

def read_input_ccs_bam(deepconsensus_dic):
    for line in sys.stdin:
        line = line.rstrip("\n")
        if line.startswith("@"):
            print(line)
        else:
            ll = line.split("\t")
            zmw = ll[0].split("/")[1]
            strand = ll[0].split("/")[3]
            if zmw in deepconsensus_dic and strand in deepconsensus_dic[zmw]:
                rname = ll[0].split("/")[0] + "/" + ll[0].split("/")[1] + "/deepconsensus/" + ll[0].split("/")[3]
                seq_baseq = deepconsensus_dic[zmw][strand]
                out_line = rname + "\t" + "\t".join(ll[1:9]) + "\t" + seq_baseq + "\t" + "\t".join(ll[11:])
                print(out_line)

if __name__ == "__main__":
    main()