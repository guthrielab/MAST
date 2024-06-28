import pandas as pd
import sys

input_file = sys.argv[1]
output_base_name = sys.argv[2]
output_dir = sys.argv[3]

df = pd.read_csv(input_file, sep='\t')
df_mutations = pd.read_csv('/Users/maximfedorov/Downloads/Lab/Summer_Project/all_resistant_variants.csv')

df_mutations['position_end'] = df_mutations.apply(
    lambda row: row['position'] + len(row['reference_nucleotide']),
    axis=1)

def find_matching_range(row, ranges):
    matching_row = ranges[
        (ranges['position'] == row['POS']) & 
        (ranges['reference_nucleotide'] == row['REF']) &
        (ranges['alternative_nucleotide'] == row['ALT'])
    ]
    return matching_row if not matching_row.empty else None

merged_rows = []
for _, pos_row in df.iterrows():
    match = find_matching_range(pos_row, df_mutations)
    if match is not None:
        for _, range_row in match.iterrows():
            merged_row = pos_row.to_dict()
            merged_row.update(range_row.to_dict())
            merged_rows.append(merged_row)

print(merged_rows)
df_merged = pd.DataFrame(merged_rows)
df_merged.to_csv(f'{output_dir}/{output_base_name}_report.csv', index=False)
