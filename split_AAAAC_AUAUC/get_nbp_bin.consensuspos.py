import sys
import pandas as pd


def main():
    in_refpos1 = sys.argv[1] # AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.combine.BEL_LTR.LTR_gene.breakpoint
    bin_size = int(sys.argv[2]) # 5

    # out = ".".join(in_refpos1.split("/")[-1].split(".")[:-1]) + "." + str(bin_size) + "bp.bin.refpos1"
    out = in_refpos1 + "." + str(bin_size) + "bp.bin.refpos"
    get_binsize_refpos1(in_refpos1, bin_size, out)


def get_binsize_refpos1(in_refpos1, bin_size, out):
    df = pd.read_csv(in_refpos1, sep='\t', header=None)
    df.columns = ['LTR', 'Strand', 'Pos', 'Count']
    binned_results = []

    for (ltr, strnad), group in df.groupby(['LTR', 'Strand']):
        binned_group = bin_data(group, bin_size)
        binned_results.append(binned_group)

    final_binned_df = pd.concat(binned_results, ignore_index=True)
    final_binned_df.to_csv(out, sep='\t', index=False, header=False)


def bin_data(df, bin_size):
    binned_data = []

    # Get the minimum and maximum positions
    min_position = df['Pos'].min()
    max_position = df['Pos'].max()

    # Iterate over the positions in non-overlapping bins of size 'bin_size'
    for start in range(min_position, max_position + 1, bin_size):
        end = start + bin_size
        # print(str(start) + "\t" + str(end))
        
        # Select the rows within this bin
        bin_df = df[(df['Pos'] >= start) & (df['Pos'] < end)]

        if not bin_df.empty:
            # Sum Count
            sum_count = bin_df['Count'].sum()
            
            # Calculate the median position
            median_position = bin_df['Pos'].median()
            
            ltr = bin_df['LTR'].iloc[0]
            strand = bin_df['Strand'].iloc[0]

            # Append the aggregated data to the list
            binned_data.append([ltr, strand, median_position, sum_count])

    # Convert the binned data back into a DataFrame
    binned_df = pd.DataFrame(binned_data, columns=['LTR', 'Strand', 'Pos', 'sum_count'])
    
    return binned_df


if __name__ == "__main__":
    main()
