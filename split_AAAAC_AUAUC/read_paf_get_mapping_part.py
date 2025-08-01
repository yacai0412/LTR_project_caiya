import sys

def main():
    gene_te_paf = sys.argv[1]

    out = gene_te_paf + ".tab"

    paf_read_pos_dic = read_paf(gene_te_paf)
    output_parts_info(paf_read_pos_dic, out)



def read_paf(gene_te_paf):
    paf_read_pos_dic = {}
    with open(gene_te_paf, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0] + "\t" + ll[1]
            query_pos = ll[2] + "-" + ll[3]
            target_pos = ll[7] + "-" + ll[8]
            strand = ll[4]
            target_type = ll[5].split("|")[1]
            if target_type == "gene":
                targetname = ll[5].split("|")[0].split("_")[1]
            else:
                targetname = ll[5].split("|")[0]
            paf_read_pos_dic.setdefault(rname, {})[query_pos] = targetname, target_pos, strand, target_type # combine smae position mapping to same gene diff transcript
    return paf_read_pos_dic


def output_parts_info(paf_read_pos_dic, out):
    with open(out, "w") as final:
        for rname in paf_read_pos_dic:
            out_query_len = []
            out_query_pos = []
            out_target_name = []
            out_target_len = []
            out_target_pos = []
            out_strand = []
            out_type = []
            for query_pos in paf_read_pos_dic[rname]:
                query_len = str(int(query_pos.split("-")[1]) - int(query_pos.split("-")[0]))
                out_query_len.append(query_len)
                out_query_pos.append(query_pos)
                out_target_name.append(paf_read_pos_dic[rname][query_pos][0])
                out_target_pos.append(paf_read_pos_dic[rname][query_pos][1])
                target_len = str(int(paf_read_pos_dic[rname][query_pos][1].split("-")[1]) - int(paf_read_pos_dic[rname][query_pos][1].split("-")[0]))
                out_target_len.append(target_len)
                out_strand.append(paf_read_pos_dic[rname][query_pos][2])
                out_type.append(paf_read_pos_dic[rname][query_pos][3])
                
            out_line = rname + "\t" + ",".join(out_query_pos) + "\t" + ",".join(out_query_len) + "\t" + ",".join(out_target_name) + "\t" + ",".join(out_target_pos) + "\t" + ",".join(out_target_len) + "\t" + ",".join(out_strand) + "\t" + ",".join(out_type)
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
