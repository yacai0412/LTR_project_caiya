import sys
import pandas as pd

def main():
    intab = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.tab

    out_ecc1LTR = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".ecc1LTR.multi_I.tab"
    out_ecc2LTR = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".ecc2LTR.multi_I.tab"


    data = pd.read_csv(intab, delimiter='\t', header=None, 
                    names=['read_name', 'read_length', 'query_span', 'query_part_length', 
                            'target_name', 'target_span', 'target_length', 'strand', 'type'])

    ecc1LTR_data, ecc2LTR_data = classify_ecc(data)
    pd.DataFrame(ecc1LTR_data).to_csv(out_ecc1LTR, index=False, header=False, sep='\t')
    pd.DataFrame(ecc2LTR_data).to_csv(out_ecc2LTR, index=False, header=False, sep='\t')


def normalize_name(name):
    if 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        return name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0]  # Extract the base name before '_'


def get_ltr_i_types(name):
    if "LTR" in name:
        return "LTR"
    else:
        return "I"

def classify_ecc(data):
    ecc1LTR_data = []
    ecc2LTR_data = []

    for index, row in data.iterrows():
        targets = row['target_name'].split(',')
        strands = row['strand'].split(',')
        types = row['type'].split(',')

        normalized_names = [normalize_name(t) for t in targets]
        unique_names = set(normalized_names)
        
        ltr_i_types = [get_ltr_i_types(t) for t in targets]
        combined_ltr_i_type = ",".join(ltr_i_types)
        row['out_ltr_i_type'] = combined_ltr_i_type

        # Check if all strands are the same and only one unique LTR name
        if all(s == types[0] for s in types) and all(s == strands[0] for s in strands) and len(unique_names) == 1:
            ltr_positions = [i for i, t in enumerate(targets) if 'LTR' in t]
            if ltr_positions and all(0 < pos < len(targets) - 1 for pos in ltr_positions):
                ltr_count = sum('LTR' in t for t in targets)
                i_count = sum('_I' in t or '-I' in t for t in targets)

                if ltr_count == 1 and i_count >= 2:
                    ecc1LTR_data.append(row)
                elif ltr_count == 2 and i_count >= 2 and ltr_positions[1] - ltr_positions[0] ==  1:
                    ecc2LTR_data.append(row)
    return ecc1LTR_data, ecc2LTR_data


if __name__ == "__main__":
    main()
