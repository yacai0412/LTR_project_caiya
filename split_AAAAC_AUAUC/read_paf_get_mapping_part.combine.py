import sys

def main():
    gene_te_paf = sys.argv[1]

    out = gene_te_paf.split("/")[-1] + ".combine.tab"

    paf_read_pos_dic = read_paf(gene_te_paf)
    merged_data = merge_overlapping_part(paf_read_pos_dic)
    output_parts_info(merged_data, out)





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


def merge_overlapping_part(data):
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
            if current_end - start > 50 and details[0] == current_details[0]:
                # If there is an overlap and conditions match, merge the intervals
                current_end = max(current_end, end)  # Extend the end if necessary
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
            out_query_pos = []
            out_target_name = []
            out_target_len = []
            out_target_pos = []
            out_strand = []
            out_type = []
            for query_pos in merged_data[rname]:
                query_len = str(int(query_pos.split("-")[1]) - int(query_pos.split("-")[0]))
                out_query_len.append(query_len)
                out_query_pos.append(query_pos)
                out_target_name.append(merged_data[rname][query_pos][0])
                out_target_pos.append(merged_data[rname][query_pos][1])
                target_len = str(int(merged_data[rname][query_pos][1].split("-")[1]) - int(merged_data[rname][query_pos][1].split("-")[0]))
                out_target_len.append(target_len)
                out_strand.append(merged_data[rname][query_pos][2])
                out_type.append(merged_data[rname][query_pos][3])
            
            
            normalized_target_names = [normalize_target_name(t) for t in out_target_name]
            ltr_i_gene_types = [get_ltr_i_types(t) for t in out_target_name]

            out_line = rname + "\t" + ",".join(out_query_pos) + "\t" + ",".join(out_query_len) + "\t" + ",".join(out_target_name) + "\t" + ",".join(out_target_pos) + "\t" + ",".join(out_target_len) + "\t" + ",".join(out_strand) + "\t" + ",".join(out_type) + "\t" + ",".join(normalized_target_names) + "\t" + ",".join(ltr_i_gene_types)
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
