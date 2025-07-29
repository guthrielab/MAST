import pandas as pd
import sys
import os
from Bio import SeqIO
from docx import Document
from jinja2 import Template

# Declaring inputs
input_file        = sys.argv[1]
output_base_name  = sys.argv[2]
patient_dir       = sys.argv[3]          
fasta_file        = sys.argv[4]
resistances       = {}

# Read tables (relative to cwd)
df           = pd.read_csv(input_file, sep='\t')
df_mutations = pd.read_csv('../../../Data/all_resistant_variants.csv')
df_lineage   = pd.read_csv('../../../Data/Lineage.csv')

# Normalize all DataFrames columns to uppercase
df.rename(columns=str.upper, inplace=True)
df_mutations.rename(columns=str.upper, inplace=True)
df_lineage.rename(columns=str.upper, inplace=True)

def find_matching_range(row, ranges):
    match = ranges[
        (ranges['POSITION'] == row['POS']) & 
        (ranges['REFERENCE_NUCLEOTIDE'] == row['REF']) &
        (ranges['ALTERNATIVE_NUCLEOTIDE'] == row['ALT'])
    ]
    return match if not match.empty else None

for _, pos_row in df.iterrows():
    match = find_matching_range(pos_row, df_mutations)
    if match is not None:
        for _, range_row in match.iterrows():
            variant = range_row['VARIANT']
            drug    = range_row['DRUG']
            resistances[variant] = drug 

for _, pos_row in df.iterrows():
    match = find_matching_range(pos_row, df_lineage)
    if match is not None:
        for _, range_row in match.iterrows():
            resistances['Lineage'] = range_row['LIN']

# Load sequence from FASTA
sequences = [str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")]
sequence  = sequences[0]

def get_complement(seq):
    comp_map = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(comp_map.get(n, n) for n in seq)

def translate_gene_sequence(sequence):
    start_codons = {'TAC','GTG','AAC','TTG','ATG'}
    if sequence is True:
        return 'indel'
    # find start
    start = next((i for i in range(0, len(sequence)-2, 3)
                  if sequence[i:i+3].upper() in start_codons), None)
    if start is None:
        return 'indel'
    # codon map
    aa_to_codons = {
        'A':['GCT','GCC','GCA','GCG'], 'C':['TGT','TGC'], 'D':['GAT','GAC'],
        'E':['GAA','GAG'], 'F':['TTT','TTC'], 'G':['GGT','GGC','GGA','GGG'],
        'H':['CAT','CAC'], 'I':['ATT','ATC','ATA'], 'K':['AAA','AAG'],
        'L':['TTA','TTG','CTT','CTC','CTA','CTG'], 'M':['ATG'], 'N':['AAT','AAC'],
        'P':['CCT','CCC','CCA','CCG'], 'Q':['CAA','CAG'],
        'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
        'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
        'T':['ACT','ACC','ACA','ACG'], 'V':['GTT','GTC','GTA','GTG'],
        'W':['TGG'], 'Y':['TAT','TAC'], '*':['TAA','TAG','TGA']
    }
    codon_to_aa = {codon: aa for aa, codons in aa_to_codons.items() for codon in codons}
    aas = []
    for i in range(start, len(sequence)-2, 3):
        aas.append(codon_to_aa.get(sequence[i:i+3].upper(), '?'))
    return ''.join(aas)

def get_mutated_sequence(df, positions, seq):
    df_gene = df[(df['POS'] >= positions[0]) & (df['POS'] <= positions[1])]
    for _, row in df_gene.iterrows():
        if 'del' in row['INFO'] or 'ins' in row['INFO']:
            return True
        pos = int(row['POS']) - 1
        seq = seq[:pos] + row['ALT'] + seq[pos+len(row['ALT']):]
    return seq[positions[0]:positions[1]]

# Apply mutations & complements
katG = get_mutated_sequence(df, [2153888,2156111], sequence)
if katG is not True: katG = get_complement(katG[::-1])

ethA = get_mutated_sequence(df, [4326003,4327473], sequence)
if ethA is not True: ethA = get_complement(ethA[::-1])

gid = get_mutated_sequence(df, [4407527,4408202], sequence)
if gid is not True: gid = get_complement(gid[::-1])

pncA = get_mutated_sequence(df, [2288680,2289241], sequence)
if pncA is not True: pncA = get_complement(pncA[::-1])

rpoB = get_mutated_sequence(df, [759806,763325], sequence)

gene_dict = {'katG': katG, 'ethA': ethA, 'gid': gid, 'pncA': pncA, 'rpoB': rpoB}
genes_map = {
    'katG':'Isoniazid','ethA':'Ethionamide','gid':'Streptomycin',
    'pncA':'Pyrazinamide','rpoB':'Rifampicin'
}

for key, seq_val in gene_dict.items():
    aa_seq = translate_gene_sequence(seq_val)
    if key == 'rpoB':
        ref_aa = translate_gene_sequence(sequence[759806:763325])
        if aa_seq != ref_aa or aa_seq == 'indel':
            resistances[key] = genes_map[key]
    else:
        if aa_seq.count('*') > 1 or aa_seq == 'indel':
            resistances[key] = genes_map[key]

# Load and fill DOCX template
template_path      = '../../../Data/Report_Template.docx'
patient_info_path  = '../../../Data/patient_info.csv'
output_dir         = patient_dir

os.makedirs(output_dir, exist_ok=True)
doc = Document(template_path)

# Patient lookup
df_pat = pd.read_csv(patient_info_path)
row   = df_pat[df_pat['Barcode'] == output_base_name]
if row.empty:
    print(f"No record found with Barcode: {output_base_name}")
    sys.exit(1)
patient_info = row.to_dict(orient='records')[0]

# Prepare context
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
        status_fields[drug]    = 'Resistant'
        status_fields[drug + '_g'] = mut
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


out_path = os.path.join(output_dir, f'{output_base_name}_report.docx')
doc.save(out_path)
print(f"Saved DOCX file to: {out_path}")




