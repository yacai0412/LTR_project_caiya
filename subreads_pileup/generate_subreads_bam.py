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
                    deepconsensus_dic[zmw] = ""
                elif c % 4 == 2:
                    deepconsensus_dic[zmw] = line
                elif c % 4 == 0:
                    deepconsensus_dic[zmw] = deepconsensus_dic[zmw] + "\t" + line
    return deepconsensus_dic

def read_input_ccs_bam(deepconsensus_dic):
    out_zmw_dic = {}
    for line in sys.stdin:
        line = line.rstrip("\n")
        if line.startswith("@"):
            print(line)
        else:
            ll = line.split("\t")
            zmw = ll[0].split("/")[1]
            if zmw in deepconsensus_dic and zmw not in out_zmw_dic:
                seq_baseq = deepconsensus_dic[zmw]
                read_length = len(seq_baseq.split("\t")[0])
                rname = ll[0].split("/")[0] + "/" + ll[0].split("/")[1] + "/0_" + str(read_length)

                out_line = rname + "\t" + "\t".join(ll[1:9]) + "\t" + seq_baseq + "\t" + "\t".join(ll[11:])
                print(out_line)
                out_zmw_dic[zmw] = 1

if __name__ == "__main__":
    main()