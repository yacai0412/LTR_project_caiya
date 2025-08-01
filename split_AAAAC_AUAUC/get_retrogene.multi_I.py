import sys
import pandas as pd


def main():
    intab = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.tab

    out_retrogene = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".retrogene.multi_I.tab"


    data = pd.read_csv(intab, delimiter='\t', header=None, 
                    names=['read_name', 'read_length', 'query_span', 'query_part_length', 
                            'target_name', 'target_span', 'target_length', 'strand', 'type1', 'name', 'type2'])

    retrogene_ls = classify_retrogene(data)
    pd.DataFrame(retrogene_ls).to_csv(out_retrogene, index=False, header=False, sep='\t')



def normalize_LTR_name(name):
    if 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        name0 = name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0].split("-")[0].split("_")[0]
        return name0  # Extract the base name before '_'


def get_ltr_i_types(name):
    if "LTR" in name:
        return "LTR"
    elif "-I" in name or "_I" in name:
        return "I"
    elif "FBgn" in name:
        return "gene"
    else:
        return "unknown"


def classify_retrogene(data):
    retrogene_ls = []
    for index, row in data.iterrows():
        targets = row['target_name'].split(',')
        strands = row['strand'].split(',')
        type1s = row['type1'].split(',')
        type2s = row['type2'].split(',')

        rname = row['read_name'].split("/")[0] + "/" + row['read_name'].split("/")[1] + "/deepconsensus/" + row['read_name'].split("/")[3]
        row['read_name'] = rname

        normalized_LTR_names = [normalize_LTR_name(t) for t in targets if not t.startswith("FBgn")]
        unique_LTR_names = set(normalized_LTR_names)
        
        # Check if all strands are the same and only one unique LTR name
        # if "TE" in type1s and "gene" in type1s and all(s == strands[0] for s in strands) and len(unique_LTR_names) == 1 and "unknown" not in type2s:
        if "TE" in type1s and "gene" in type1s and len(unique_LTR_names) == 1 and "unknown" not in type2s:
            ltr_positions = [i for i, t in enumerate(type2s) if t == "LTR"]
            gene_positions = [i for i, t in enumerate(type2s) if t == "gene"]

            if len(gene_positions) != 0 and all((pos1 != 0 and pos1 != len(targets) -1) for pos1 in gene_positions):
                if len(ltr_positions) != 0:
                    if ltr_positions and all((pos == 0 or pos == len(targets) -1) for pos in ltr_positions):
                        retrogene_ls.append(row)
                else:
                    retrogene_ls.append(row)

    return retrogene_ls

if __name__ == "__main__":
    main()



