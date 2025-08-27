#!/usr/bin/env python3
import pandas as pd
import sys
import os
from Bio import SeqIO        
from docx import Document
from jinja2 import Template
import pysam

input_file        = sys.argv[1]
output_base_name  = sys.argv[2]
patient_dir       = sys.argv[3]
fasta_file        = sys.argv[4]   
resistances       = {}

df           = pd.read_csv(input_file, sep='\t')
df_mutations = pd.read_csv('../../../MAST/Data/all_resistant_variants.csv')
df_lineage   = pd.read_csv('../../../MAST/Data/Lineage.csv')

df.columns           = df.columns.str.upper()
df_mutations.columns = df_mutations.columns.str.upper()
df_lineage.columns   = df_lineage.columns.str.upper()

merged = pd.merge(
    df,
    df_mutations,
    left_on  = ['POS','REF','ALT'],
    right_on = ['POSITION','REFERENCE_NUCLEOTIDE','ALTERNATIVE_NUCLEOTIDE'],
    how      = 'inner'
)
for _, row in merged.iterrows():
    resistances[row['VARIANT']] = row['DRUG']

merged_lin = pd.merge(
    df,
    df_lineage,
    left_on  = ['POS','REF','ALT'],
    right_on = ['POSITION','REFERENCE_NUCLEOTIDE','ALTERNATIVE_NUCLEOTIDE'],
    how      = 'inner'
)

# Check for lineage-specific reference position for samples without lineage
lineage_detected = False
if not merged_lin.empty:
    resistances['Lineage'] = merged_lin['LIN'].iloc[0]
    lineage_detected = True

# If no lineage detected from mutation list, check ref position
if not lineage_detected:
    lineage_reference_positions = {
        420008: ('4.9', 'A')
    }
    # Search for BAM file in parent dirs
    bam_file = None
    bam_patterns = [
        f"*{output_base_name}*.bam",
        f"aligned_sorted_{output_base_name}.bam",
        f"*sorted*.bam"
    ]
    search_dirs = ['.', patient_dir]
    current_dir = os.getcwd()
    for i in range(1,4):
        parent_dir = os.path.abspath(os.path.join(current_dir, *['..'] *i))
        if os.path.exists(parent_dir):
            search_dirs.append(parent_dir)
    work_dirs = []
    for search_dir in search_dirs[:]:
        if os.path.exists(search_dir):
            for root, dirs, files in os.walk(search_dir):
                if 'work' in dirs:
                    work_dirs.append(os.path.join(root, 'work'))
                if any(f.endswith('.bam') for f in files):
                    work_dirs.append(root)
    search_dirs.extend(work_dirs)
    search_dirs = list(set([d for d in search_dirs if os.path.exists(d)]))

    print(f"Searching for BAM files in directories: {search_dirs}")

    for search_dir in search_dirs:
        for pattern in bam_patterns:
            try:
                for root, dirs, files in os.walk(search_dir):
                    for file in files:
                        if (file.endswith('.bam') and
                            output_base_name in file and
                            'sorted' in file.lower()):
                            potential_bam = os.path.join(root, file)
                            try:
                                with pysam.AlignmentFile(potential_bam, "rb") as test_bam:
                                    if test_bam.nreferences > 0:
                                        bam_file = potential_bam
                                        print(f"Found BAM file: {bam_file}")
                                        break
                            except:
                                continue
                    if bam_file:
                        break
                if bam_file:
                    break
            except (PermissionError, OsError) as e:
                print(f"Cannot access {search_dir}: {e}")
                continue
        if bam_file:
            break
    if bam_file and os.path.exists(bam_file):
        print(f"Using BAM file for lineage reference check: {bam_file}")

        try:
            # Check if BAM file has index and/or create one
            bai_file = bam_file + '.bai'
            if not os.path.exists(bai_file):
                print(f"Indexing BAM file: {bam_file}")
                pysam.index(bam_file)
            # Open BAM file
            bam = pysam.AlignmentFile(bam_file, "rb")
            for pos, (lineage, expected_base) in lineage_reference_positions.items():
                try:
                    # Check coverage at ref position
                    coverage = bam.count_coverage(contig=bam.references[0], start=pos-1, stop=pos)
                    base_counts = {
                        'A': coverage[0][0],
                        'C': coverage[1][0],
                        'G': coverage[2][0],
                        'T': coverage[3][0]
                    }

                    total_bases = sum(base_counts.values())

                    print(f"Position {pos}: coverage={total_bases}, bases={base_counts}")

                    if total_bases >= 10: # Minimum coverage threshold
                        consensus = max(base_counts.items(), key=lambda x: x[1])[0]
                        consensus_fraction = base_counts[consensus] / total_bases

                        # Check if the consensus matches the expected reference base
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

    if not lineage_detected:
        print("No lineage detected from mutations or reference positions")

if not os.path.isabs(patient_dir):
    cwd = os.getcwd()
    marker = os.sep + 'work' + os.sep
    if marker in cwd:
        pipeline_root = cwd.split(marker)[0]
        patient_dir = os.path.join(pipeline_root, patient_dir)
    else:
        patient_dir = os.path.abspath(patient_dir)

template_path     = '../../../MAST/Data/Report_Template.docx'
patient_info_path = '../../../MAST/Data/patient_info.csv'
os.makedirs(patient_dir, exist_ok=True)

doc = Document(template_path)

df_pat = pd.read_csv(patient_info_path)
row    = df_pat[df_pat['Barcode'] == output_base_name]
if row.empty:
    print(f"No record found with Barcode: {output_base_name}", file=sys.stderr)
    sys.exit(1)
patient_info = row.to_dict(orient='records')[0]

status_fields = {
    'Ethambutol':'Susceptible','Ethambutol_g':'None',
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

# Create TSV output data structure
tsv_data = {
    'Sample': output_base_name,
    'Lineage': 'Unknown'
}

# Initialize all drugs as susceptible in TSV
drugs = ['Ethambutol', 'Pyrazinamide', 'Isoniazid', 'Rifampicin', 'Streptomycin',
        'Ciprofloxacin', 'Ofloxacin', 'Moxifloxacin', 'Amikacin', 'Kanamycin',
        'Capreomycin', 'Bedaquiline', 'Linezolid']

for drug in drugs:
    tsv_data[f'{drug}_Status'] = 'Susceptible'
    tsv_data[f'{drug}_Mutation'] = 'None'

for mut, drug in resistances.items():
    if drug in status_fields:
        status_fields[drug]        = 'Resistant'
        status_fields[f'{drug}_g'] = mut
        tsv_data[f'{drug}_Status'] = 'Resistant'
        tsv_data[f'{drug}_Mutation'] = mut
    elif mut == 'Lineage':
        status_fields['Lineage'] = drug
        tsv_data['Lineage'] = drug

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
tsv_df = pd.DataFrame([tsv_data])
tsv_path = os.path.join(patient_dir, f'{output_base_name}_results.tsv')
tsv_df.to_csv(tsv_path, sep='\t', index=False)

# Print summary to console
print(f"\nResistance Profile Summary for {output_base_name}:")
print(f"Lineage: {tsv_data['Lineage']}")
print("\nDrug Resistance Status:")
for drug in drugs:
    status = tsv_data[f'{drug}_Status']
    mutation = tsv_data[f'{drug}_Mutation']
    if status == 'Resistant':
        print(f"  {drug}: {status} (Mutation: {mutation})")
    else:
        print(f"  {drug}: {status}")





