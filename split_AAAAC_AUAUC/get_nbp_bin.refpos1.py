import sys
import pandas as pd


def main():
    in_refpos1 = sys.argv[1]
    bin_size = int(sys.argv[2])

    out = ".".join(in_refpos1.split("/")[-1].split(".")[:-1]) + "." + str(bin_size) + "bp.bin.refpos1"
    get_binsize_refpos1(in_refpos1, bin_size, out)


def get_binsize_refpos1(in_refpos1, bin_size, out):
    df = pd.read_csv(in_refpos1, sep='\t', header=None)
    df.columns = ['Label', 'Type', 'LTR', 'Position', 'mutation_count', 'depth', 'mutation_rate']
    binned_results = []

    for (label, type_, ltr), group in df.groupby(['Label', 'Type', 'LTR']):
        binned_group = bin_data(group, bin_size)
        binned_results.append(binned_group)

    final_binned_df = pd.concat(binned_results, ignore_index=True)
    final_binned_df.to_csv(out, sep='\t', index=False, header=False)


def bin_data(df, bin_size):
    binned_data = []

    # Get the minimum and maximum positions
    min_position = df['Position'].min()
    max_position = df['Position'].max()

    # Iterate over the positions in non-overlapping bins of size 'bin_size'
    for start in range(min_position, max_position, bin_size):
        end = start + bin_size
        
        # Select the rows within this bin
        bin_df = df[(df['Position'] >= start) & (df['Position'] < end)]

        if not bin_df.empty:
            # Sum Col5 and Col6 in this bin
            sum_col5 = bin_df['mutation_count'].sum()
            sum_col6 = bin_df['depth'].sum()

            if sum_col6 == 0:
                mean_mutation_rate = 0
            else:
                mean_mutation_rate = sum_col5 / sum_col6
            
            # Calculate the median position
            median_position = bin_df['Position'].median()
            
            # Use the first Label, Type, and LTR (assuming they are constant)
            label = bin_df['Label'].iloc[0]
            type_ = bin_df['Type'].iloc[0]
            ltr = bin_df['LTR'].iloc[0]
            
            # Append the aggregated data to the list
            binned_data.append([label, type_, ltr, median_position, sum_col5, sum_col6, mean_mutation_rate])

    # Convert the binned data back into a DataFrame
    binned_df = pd.DataFrame(binned_data, columns=['Label', 'Type', 'LTR', 'Position', 'sum_mutation_count', 'sum_depth', 'mean_mutation_rate'])
    
    return binned_df


if __name__ == "__main__":
    main()
