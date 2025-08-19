#!/usr/bin/env python3
import pandas as pd
import sys
import os
from Bio import SeqIO        
from docx import Document
from jinja2 import Template

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
if not merged_lin.empty:
    resistances['Lineage'] = merged_lin['LIN'].iloc[0]


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
    'Capreomycin':'Susceptible','Capreomycin_g':'None'
}


for mut, drug in resistances.items():
    if drug in status_fields:
        status_fields[drug]        = 'Resistant'
        status_fields[f'{drug}_g'] = mut
    elif mut == 'Lineage':
        status_fields['Lineage'] = drug

context = {**patient_info, **status_fields}

for para in doc.paragraphs:
    if '{{' in para.text:
        para.text = Template(para.text).render(context)

for table in doc.tables:
    for row in table.rows:
        for cell in row.cells:
            if '{{' in cell.text:
                cell.text = Template(cell.text).render(context)

out_path = os.path.join(patient_dir, f'{output_base_name}_report.docx')
doc.save(out_path)
print(f"Saved DOCX file to: {out_path}")







