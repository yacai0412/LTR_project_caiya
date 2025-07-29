import sys

def main():
    asm_paf = sys.argv[1] # 

    out = asm_paf.split("/")[-1] + ".asm.combine.tab"

    paf_read_pos_dic = read_paf(asm_paf)
    merged_data = merge_overlapping_asm_part(paf_read_pos_dic)
    output_parts_info(merged_data, out)


def export_full_cover_larger_range(query_pos0, query_pos1):
    query_start0 = int(query_pos0.split("-")[0])
    query_end0 = int(query_pos0.split("-")[1])
    query_start1 = int(query_pos1.split("-")[0])
    query_end1 = int(query_pos1.split("-")[1])
    if query_start0 <= query_start1 and query_end0 >= query_end1:
        return query_pos0
    elif query_start0 <= query_start1 and query_end0 >= query_end1:
        return query_pos1
    else:
        return 0


def read_paf(asm_paf):
    paf_read_pos_dic = {}
    with open(asm_paf, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0] + "\t" + ll[1]
            query_pos = ll[2] + "-" + ll[3]
            target_pos = ll[7] + "-" + ll[8]
            strand = ll[4]
            asm_chr = ll[5]

            if rname not in paf_read_pos_dic:
                paf_read_pos_dic.setdefault(rname, {})[query_pos] = asm_chr, target_pos, strand # combine smae position mapping to same gene diff transcript
            else:
                c = 0
                for query_pos0 in paf_read_pos_dic[rname]:
                    larger_query_pos = export_full_cover_larger_range(query_pos0, query_pos)
                    if larger_query_pos == query_pos0:
                        c += 1
                        break
                    elif larger_query_pos == query_pos:
                        paf_read_pos_dic[rname][query_pos0] = ""
                        paf_read_pos_dic.setdefault(rname, {})[query_pos] = asm_chr, target_pos, strand # combine smae position mapping to same gene diff transcript
                        c += 1
                        break
                    elif larger_query_pos == 0:
                        continue
                if c == 0:
                    paf_read_pos_dic.setdefault(rname, {})[query_pos] = asm_chr, target_pos, strand # combine smae position mapping to same gene diff transcript

    paf_read_pos_dic1 = {}
    for rname in paf_read_pos_dic:
        for query_pos in paf_read_pos_dic[rname]:
            if paf_read_pos_dic[rname][query_pos] != "":
                paf_read_pos_dic1.setdefault(rname, {})[query_pos] = paf_read_pos_dic[rname][query_pos]

    return paf_read_pos_dic1


def merge_overlapping_asm_part(data):
    merged_data = {}

    for read, positions in data.items():
        intervals = []

        for pos, details in positions.items():
            start, end = map(int, pos.split('-'))
            intervals.append((start, end, details))
        intervals.sort()

        merged = []
        current_start, current_end, current_details = intervals[0]

        for i in range(1, len(intervals)):
            start, end, details = intervals[i]

            # Check if there is an overlap and if the target names and strands are the same
            # if current_end - start > 50 and details[0] == current_details[0]: # same gene or TE (details[0]: gene name)
                # If there is an overlap and conditions match, merge the intervals
            if current_end - start > 50:
                left_range = (current_start, start)
                rigth_range = (current_end, end)

                left_length = left_range[1] - left_range[0]
                right_length = rigth_range[1] - rigth_range[0]

                if left_length < 50 and right_length < 50:
                    current_end = max(current_end, end)
                elif left_length < 50 and right_length >= 50:
                    merged.append((f"{current_start}-{current_end}", current_details))
                    current_start = current_end
                    current_end, current_details = end, details
                else: # (left_length >= 50 and right_length < 50) or (left_length >= 50 and right_length >= 50):
                    merged.append((f"{current_start}-{start}", current_details))
                    current_start, current_end, current_details = start, end, details

            else:
                # No sufficient overlap or different conditions, save the current interval
                merged.append((f"{current_start}-{current_end}", current_details))
                current_start, current_end, current_details = start, end, details


        # Add the last interval
        merged.append((f"{current_start}-{current_end}", current_details))

        # Store merged intervals into the new dictionary
        merged_data[read] = {pos: details for pos, details in merged}

    return merged_data




def normalize_target_name(name):
    if 'FBgn' in name:
        return name
    elif 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        return name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0]  # Extract the base name before '_'

def get_ltr_i_types(name):
    if "LTR" in name:
        return "LTR"
    elif "-I" in name or "_I" in name:
        return "I"
    elif "FBgn" in name:
        return "gene"
    else:
        return "unknown"

def output_parts_info(merged_data, out):
    with open(out, "w") as final:
        for rname in merged_data:
            out_query_len = []
            out_query_range = []
            out_target_name = []
            out_target_range = []
            out_target_len = []
            # out_target_pos = []
            out_strand = []
            # out_te_gene_type = []
            # out_name = []
            # out_LTR_I_gene = []
            
            for query_pos in merged_data[rname]: # asm_chr, target_pos, strand
                asm_chr, target_pos, strand = merged_data[rname][query_pos]
                # te_gene_type = asm_chr.split("|")[1]

                query_len = str(int(query_pos.split("-")[1]) - int(query_pos.split("-")[0]))
                out_query_len.append(query_len)
                out_query_range.append(query_pos)
                # out_target_name.append(merged_data[rname][query_pos][0])

                # target_chr_pos = asm_chr + ":" + target_pos
                # out_target_pos.append(target_chr_pos)

                out_target_name.append(asm_chr.split("|")[0])
                out_target_range.append(target_pos)

                target_len = str(int(target_pos.split("-")[1]) - int(target_pos.split("-")[0]))
                out_target_len.append(target_len)
                out_strand.append(strand)
            
                # out_te_gene_type.append(te_gene_type)
                # out_name.append(normalize_target_name(asm_chr.split("|")[0]))
                # out_LTR_I_gene.append(get_ltr_i_types(asm_chr.split("|")[0]))


            out_line = rname + "\t" + ",".join(out_query_range) + "\t" + ",".join(out_query_len) + "\t" + ",".join(out_target_name) + "\t" + ",".join(out_target_range) + "\t" + ",".join(out_target_len) + "\t" + ",".join(out_strand)
            # out_line = rname + "\t" + ",".join(out_query_range) + "\t" + ",".join(out_query_len) + "\t" + ",".join(out_target_name) + "\t" + ",".join(out_target_range) + "\t" + ",".join(out_target_len) + "\t" + ",".join(out_strand) + "\t" + ",".join(out_te_gene_type) + "\t" + ",".join(out_name) + "\t" + ",".join(out_LTR_I_gene)
            final.write(out_line + "\n")




if __name__ == "__main__":
    main()



