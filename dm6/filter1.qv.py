import sys

def main():
    snv = sys.argv[1] # AA.deepconsensus.cx3.pbmm2.dm6.chr.overlap90.s.bam.ds.out
    ds_ss_tag = sys.argv[2] # ds / ss
    qv_filter = int(sys.argv[3]) # 20 30
    concordence = float(1 - 10 ** (qv_filter / (-10)))

    out = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".qv" + str(qv_filter) + ".out"
    final = open(out, "w")

    read_snv(snv, ds_ss_tag, concordence, final)
    final.close()

def filter_ds_qv(concordence, ll):
    plus_reads_ls = ll[4].split(",")
    minus_reads_ls = ll[5].split(",")

    plus_reads_dic = {}
    for read in plus_reads_ls:
        zmw = read.split("/")[1]
        rq = float(read.split(":")[2])
        if rq >= concordence:
            plus_reads_dic[zmw] = read
    
    minus_reads_dic = {}
    for read in minus_reads_ls:
        zmw = read.split("/")[1]
        rq = float(read.split(":")[2])
        if rq >= concordence:
            minus_reads_dic[zmw] = read

    out_dic = {}
    out_dic.setdefault("plus", [])
    out_dic.setdefault("minus", [])
    for zmw in plus_reads_dic:
        plus_read = plus_reads_dic[zmw]
        if zmw in minus_reads_dic:
            minus_read = minus_reads_dic[zmw]
            out_dic.setdefault("plus", []).append(plus_read)
            out_dic.setdefault("minus", []).append(minus_read)
    
    out_str = ""
    if len(out_dic["plus"]) > 0 and len(out_dic["plus"]) == len(out_dic["minus"]):
        out_str = ",".join(out_dic["plus"]) + "\t" + ",".join(out_dic["minus"])
    return out_str

def filter_ss_qv(concordence, ll):
    plus_reads_ls = ll[4].split(",")
    minus_reads_ls = ll[5].split(",")

    plus_ls = []
    for read in plus_reads_ls:
        if read == "NA":
            plus_ls = []
        else:
            rq = float(read.split(":")[2])
            if rq >= concordence:
                plus_ls.append(read)
    if len(plus_ls) > 0:
        plus_str = ",".join(plus_ls)
    else:
        plus_str = "NA"

    minus_ls = []
    for read in minus_reads_ls:
        if read == "NA":
            minus_ls = []
        else:
            rq = float(read.split(":")[2])
            if rq >= concordence:
                minus_ls.append(read)
    if len(minus_ls) > 0:
        minus_str = ",".join(minus_ls)
    else:
        minus_str = "NA"

    out_str = ""
    if plus_str != minus_str and (plus_str == "NA" or minus_str == "NA"):
        out_str = plus_str + "\t" + minus_str
    return out_str

def read_snv(snv, ds_ss_tag, concordence, final):
    with open(snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            
            if ds_ss_tag == "ds":
                out_str = filter_ds_qv(concordence, ll)
                if out_str != "":
                    out_line = "\t".join(ll[0:4]) + "\t" + out_str
                    final.write(out_line + "\n")
            
            elif ds_ss_tag == "ss":
                out_str = filter_ss_qv(concordence, ll)
                if out_str != "":
                    out_line = "\t".join(ll[0:4]) + "\t" + out_str
                    final.write(out_line + "\n")

if __name__ == "__main__":
    main()
