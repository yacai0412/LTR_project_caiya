import sys

def main():
    snv0 = sys.argv[1]

    ss_out = ".".join(snv0.split("/")[-1].split(".")[:-1]) + ".ss.out"
    ds_out = ".".join(snv0.split("/")[-1].split(".")[:-1]) + ".ds.out"
    error_out = ".".join(snv0.split("/")[-1].split(".")[:-1]) + ".error.out"

    ss_final = open(ss_out, "w")
    ds_final = open(ds_out, "w")
    error_final = open(error_out, "w")
    read_snv0(snv0, ss_final, ds_final, error_final)
    
    ss_final.close()
    ds_final.close()
    error_final.close()

def get_support_zmw_dic(support_reads):
    out_dic = {}
    if support_reads == "NA":
        out_dic["NA"] = "NA"
    else:
        ls0 = support_reads.split(",")
        for read in ls0:
            zmw = read.split("/")[1]
            strand0 = read.split(":")[0].split("/")[3]
            out_dic[zmw] = strand0
    return out_dic

def which_ds_ss_error(plus_reads_dic, minus_reads_dic):
    if "NA" in plus_reads_dic or "NA" in minus_reads_dic:
        out_tag = "ss"
    else:
        ds_zmw_dic = {}
        for zmw in plus_reads_dic:
            strand0 = plus_reads_dic[zmw]
            if zmw in minus_reads_dic and strand0 != minus_reads_dic[zmw]:
                ds_zmw_dic[zmw] = 1
        
        if len(ds_zmw_dic) > 0:
            out_tag = "ds"
        else:
            out_tag = "error"
    return out_tag

def get_ds_support_reads(plus_reads_dic, minus_reads_dic, plus_support_reads, minus_support_reads):
    ds_zmw_dic = {}
    for zmw in plus_reads_dic:
        strand0 = plus_reads_dic[zmw]
        if zmw in minus_reads_dic and strand0 != minus_reads_dic[zmw]:
            ds_zmw_dic[zmw] = 1
    
    out_reads0 = ""
    for i in [plus_support_reads, minus_support_reads]:
        ds_reads_ls = []
        ls0 = i.split(",")
        for read in ls0:
            zmw = read.split("/")[1]
            if zmw in ds_zmw_dic:
                ds_reads_ls.append(read)
        out_reads0 = out_reads0 + "\t" + ",".join(ds_reads_ls)
    return out_reads0

def read_snv0(snv0, ss_final, ds_final, error_final):
    with open(snv0, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            plus_reads_dic = get_support_zmw_dic(ll[4])
            minus_reads_dic = get_support_zmw_dic(ll[5])
            ds_ss_tag = which_ds_ss_error(plus_reads_dic, minus_reads_dic)
            if ds_ss_tag == "ss":
                ss_final.write(line + "\n")
            elif ds_ss_tag == "error":
                error_final.write(line + "\n")
            elif ds_ss_tag == "ds":
                ds_reads = get_ds_support_reads(plus_reads_dic, minus_reads_dic, ll[4], ll[5])
                out_line = "\t".join(ll[0:4]) + ds_reads
                ds_final.write(out_line + "\n")

if __name__ == "__main__":
    main()


