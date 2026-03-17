#!/usr/bin/env python3
import pandas as pd
import sys
import os
from docx import Document
from jinja2 import Template
import pysam

input_file        = sys.argv[1]
output_base_name  = sys.argv[2]
patient_dir       = sys.argv[3]
fasta_file        = sys.argv[4]
mutations_csv     = sys.argv[5]
lineage_csv       = sys.argv[6]
template_docx     = sys.argv[7]
patient_info_csv  = sys.argv[8]
bam_file          = sys.argv[9]   # staged directly by Nextflow alongside its .bai

resistances = {}

df           = pd.read_csv(input_file, sep='\t')
df_mutations = pd.read_csv(mutations_csv)
df_lineage   = pd.read_csv(lineage_csv)

df.columns           = df.columns.str.upper()
df_mutations.columns = df_mutations.columns.str.upper()
df_lineage.columns   = df_lineage.columns.str.upper()

print(f"Loaded {len(df)} variants, {len(df_mutations)} mutation references, {len(df_lineage)} lineage references")

# ── Exact match (fast vectorized merge) ──────────────────────────────────────
exact_merged = pd.merge(
    df,
    df_mutations,
    left_on  = ['POS', 'REF', 'ALT'],
    right_on = ['POSITION', 'REFERENCE_NUCLEOTIDE', 'ALTERNATIVE_NUCLEOTIDE'],
    how      = 'inner'
)
for _, row in exact_merged.iterrows():
    resistances[row['VARIANT']] = row['DRUG']
    print(f"Exact match: {row['VARIANT']} at position {row['POS']} ({row['REF']}>{row['ALT']})")

print(f"Exact matches found: {len(exact_merged)}")

# ── Fuzzy indel/position-shift matching (vectorized) ─────────────────────────
# Strategy: for each TSV variant, find CSV rows within ±1 position, then
# apply the indel sub-checks. We use a merge on a position range rather than
# a nested Python loop, reducing the search space dramatically.

# Add position windows to the mutations reference table
df_mutations['POS_LOW']  = df_mutations['POSITION'] - 1
df_mutations['POS_HIGH'] = df_mutations['POSITION'] + 1

# Cross-join candidates within ±1 position using a merge on an integer key,
# then filter down — much faster than iterating every pair.
df['_key'] = 1
df_mutations['_key'] = 1
candidates = pd.merge(df, df_mutations, on='_key').drop(columns='_key')
df.drop(columns='_key', inplace=True)
df_mutations.drop(columns='_key', inplace=True)

# Keep only rows where positions are within ±1
candidates = candidates[
    (candidates['POS'] >= candidates['POS_LOW']) &
    (candidates['POS'] <= candidates['POS_HIGH'])
]

print(f"Candidate pairs after position filter: {len(candidates)}")

for _, c in candidates.iterrows():
    variant = c['VARIANT']
    if variant in resistances:
        continue  # already matched by exact merge

    tsv_pos = c['POS'];   tsv_ref = c['REF'];   tsv_alt = c['ALT']
    csv_pos = c['POSITION']; csv_ref = c['REFERENCE_NUCLEOTIDE']; csv_alt = c['ALTERNATIVE_NUCLEOTIDE']

    matched = False

    # Exact (already caught above, but included for completeness)
    if tsv_pos == csv_pos and tsv_ref == csv_ref and tsv_alt == csv_alt:
        matched = True

    # Deletion check
    elif len(str(tsv_ref)) > len(str(tsv_alt)) and len(str(csv_ref)) > len(str(csv_alt)):
        tsv_deleted = str(tsv_ref)[len(str(tsv_alt)):]
        csv_deleted = str(csv_ref)[len(str(csv_alt)):]
        if tsv_deleted in csv_deleted or csv_deleted in tsv_deleted:
            matched = True

    # Insertion check
    elif len(str(tsv_alt)) > len(str(tsv_ref)) and len(str(csv_alt)) > len(str(csv_ref)):
        tsv_inserted = str(tsv_alt)[len(str(tsv_ref)):]
        csv_inserted = str(csv_alt)[len(str(csv_ref)):]
        if tsv_inserted in csv_inserted or csv_inserted in tsv_inserted:
            matched = True

    # Position-shifted SNP check
    elif str(tsv_ref) in str(csv_ref) or str(csv_ref) in str(tsv_ref):
        if str(tsv_alt) in str(csv_alt) or str(csv_alt) in str(tsv_alt):
            matched = True

    if matched:
        resistances[variant] = c['DRUG']
        print(f"Fuzzy match: {variant} at position {tsv_pos} ({tsv_ref}>{tsv_alt}) "
              f"with CSV position {csv_pos} ({csv_ref}>{csv_alt})")

print(f"Total resistance variants found: {len(resistances)}")

# ── Lineage from mutation list ────────────────────────────────────────────────
exact_lin = pd.merge(
    df,
    df_lineage,
    left_on  = ['POS', 'REF', 'ALT'],
    right_on = ['POSITION', 'REFERENCE_NUCLEOTIDE', 'ALTERNATIVE_NUCLEOTIDE'],
    how      = 'inner'
)

lineage_detected = False
if not exact_lin.empty:
    resistances['Lineage'] = exact_lin['LIN'].iloc[0]
    lineage_detected = True
    print(f"Lineage detected from mutation list: {exact_lin['LIN'].iloc[0]}")

# ── Lineage from BAM reference position ──────────────────────────────────────
# The BAI is staged by Nextflow alongside the BAM — pysam finds it automatically.
if not lineage_detected:
    lineage_reference_positions = {
        420008: ('4.9', 'A')
    }

    if os.path.exists(bam_file):
        print(f"Checking BAM for lineage reference positions: {bam_file}")
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            for pos, (lineage, expected_base) in lineage_reference_positions.items():
                try:
                    coverage = bam.count_coverage(
                        contig=bam.references[0],
                        start=pos - 1,
                        stop=pos
                    )
                    base_counts = {
                        'A': coverage[0][0],
                        'C': coverage[1][0],
                        'G': coverage[2][0],
                        'T': coverage[3][0]
                    }
                    total_bases = sum(base_counts.values())
                    print(f"Position {pos}: coverage={total_bases}, bases={base_counts}")

                    if total_bases >= 10:
                        consensus = max(base_counts.items(), key=lambda x: x[1])[0]
                        consensus_fraction = base_counts[consensus] / total_bases
                        if consensus == expected_base and consensus_fraction >= 0.9:
                            resistances['Lineage'] = lineage
                            print(f"Detected lineage {lineage} from reference position {pos}")
                            lineage_detected = True
                            break
                except Exception as e:
                    print(f"Error checking position {pos}: {e}")
                    continue
            bam.close()
        except Exception as e:
            print(f"Error processing BAM file {bam_file}: {e}")
    else:
        print(f"BAM file not found at expected path: {bam_file}")

    if not lineage_detected:
        print("No lineage detected from mutations or reference positions")

# ── Build report ──────────────────────────────────────────────────────────────
os.makedirs(patient_dir, exist_ok=True)

doc = Document(template_docx)

df_pat = pd.read_csv(patient_info_csv)
row    = df_pat[df_pat['Barcode'] == output_base_name]
if row.empty:
    print(f"No record found with Barcode: {output_base_name}", file=sys.stderr)
    sys.exit(1)
patient_info = row.to_dict(orient='records')[0]

status_fields = {
    'Ethambutol':'Susceptible','Ethambutol_g':'None',
    'Ethionamide':'Susceptible','Ethionamide_g':'None',
    'Pyrazinamide':'Susceptible','Pyrazinamide_g':'None',
    'Isoniazid':'Susceptible','Isoniazid_g':'None',
    'Rifampicin':'Susceptible','Rifampicin_g':'None',
    'Streptomycin':'Susceptible','Streptomycin_g':'None',
    'Ciprofloxacin':'Susceptible','Ciprofloxacin_g':'None',
    'Ofloxacin':'Susceptible','Ofloxacin_g':'None',
    'Moxifloxacin':'Susceptible','Moxifloxacin_g':'None',
    'Amikacin':'Susceptible','Amikacin_g':'None',
    'Kanamycin':'Susceptible','Kanamycin_g':'None',
    'Capreomycin':'Susceptible','Capreomycin_g':'None',
    'Bedaquiline':'Susceptible','Bedaquiline_g':'None',
    'Linezolid':'Susceptible','Linezolid_g':'None'
}

drugs = [
    'Ethambutol', 'Ethionamide', 'Pyrazinamide', 'Isoniazid', 'Rifampicin',
    'Streptomycin', 'Ciprofloxacin', 'Ofloxacin', 'Moxifloxacin', 'Amikacin',
    'Kanamycin', 'Capreomycin', 'Bedaquiline', 'Linezolid'
]

tsv_data = {'Sample': output_base_name, 'Lineage': 'Unknown'}
for drug in drugs:
    tsv_data[f'{drug}_Status']   = 'Susceptible'
    tsv_data[f'{drug}_Mutation'] = 'None'

for mut, drug in resistances.items():
    if drug in status_fields:
        status_fields[drug]          = 'Resistant'
        status_fields[f'{drug}_g']   = mut
        tsv_data[f'{drug}_Status']   = 'Resistant'
        tsv_data[f'{drug}_Mutation'] = mut
    elif mut == 'Lineage':
        status_fields['Lineage'] = drug
        tsv_data['Lineage']      = drug

context = {**patient_info, **status_fields}

for para in doc.paragraphs:
    if '{{' in para.text:
        para.text = Template(para.text).render(context)

for table in doc.tables:
    for row in table.rows:
        for cell in row.cells:
            if '{{' in cell.text:
                cell.text = Template(cell.text).render(context)

# Save DOCX report
out_path = os.path.join(patient_dir, f'{output_base_name}_report.docx')
doc.save(out_path)
print(f"Saved DOCX file to: {out_path}")

# Save TSV report
tsv_df   = pd.DataFrame([tsv_data])
tsv_path = os.path.join(patient_dir, f'{output_base_name}_results.tsv')
tsv_df.to_csv(tsv_path, sep='\t', index=False)

# Print summary
print(f"\nResistance Profile Summary for {output_base_name}:")
print(f"Lineage: {tsv_data['Lineage']}")
print("\nDrug Resistance Status:")
for drug in drugs:
    status   = tsv_data[f'{drug}_Status']
    mutation = tsv_data[f'{drug}_Mutation']
    if status == 'Resistant':
        print(f"  {drug}: {status} (Mutation: {mutation})")
    else:
        print(f"  {drug}: {status}")
