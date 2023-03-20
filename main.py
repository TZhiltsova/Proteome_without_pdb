from pyfaidx import Fasta
import csv
from collections import Counter

# reading of the file
# with amyloids motifs
waltz_seq = []       # list for short motifs of amyloid
waltz_dict = {}      # dictionary for sequences from WaltzDB
with open('WALTZ_DB_amyloid_seq') as waltz_db:
    for line in waltz_db:
        waltz_seq = [line.strip() for line in waltz_db]
    for number, amiloid_motif in enumerate(waltz_seq):
        if amiloid_motif.strip().isalpha() == False:    # deleting empty lines
            break
        else:
            waltz_dict[number] = amiloid_motif


proteome_seq = {}  # dict for all sequences which are shorter than 1000 residues
with Fasta('proteome_without_pdb.fasta') as genes:
    for number_of_peptide in range((len(genes.keys()))):
        proteome_seq[genes[number_of_peptide][:].name] = genes[number_of_peptide][:]

# finding peptides with amyloid motifs
new_csv = []
amyloid_motifs_list_for_header = {}  # remembering amyloid motifs for each peptide, which have amyloid-forming region
amount_of_amyloids_motifs = {}
list_of_amyl_motifs_for_counter = []
filter_seq = {}  # dict for seq that have amyloid-forming regions
for key, seq in proteome_seq.items():
    counter_for_motifs = 0
    amyloid_motifs_list = []   # updated list for all amyloid motifs in one peptide
    for val in waltz_dict.values():
        if val in str(seq):
            counter_for_motifs += 1
            list_of_amyl_motifs_for_counter.append(val)
            filter_seq[key] = seq
            pos = str(seq).find(val)+1
            val = str(pos) + val
            amyloid_motifs_list.append(val)
            name = key.split('|')
            amount_of_amyloids_motifs[key] = counter_for_motifs
            if name[1] not in new_csv:
                new_csv.append(name[1])
        amyloid_motifs_list_for_header[key] = amyloid_motifs_list
counter_of_motifs = Counter(list_of_amyl_motifs_for_counter)

table = []
with open('uniprot-proteome.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    count = 0
    for rows in file:
        if count == 0:
            table.append(rows)
            count += 1
        for code in new_csv:
            if code in rows:
                table.append(rows)

table_short = []
t = 1
for line in table:
    entry = line.index('Entry')
    entry_name = line.index('Entry Name')
    reviewed = line.index('Reviewed')
    protein_names = line.index('Protein names')
    gene_names = line.index('Gene Names')
    organism = line.index('Organism')
    length = line.index('Length')
    break

for elem in table:
    table_short_words = []
    table_short_words.append(elem[entry])
    table_short_words.append(elem[entry_name])
    table_short_words.append(elem[reviewed])
    table_short_words.append(elem[protein_names])
    table_short_words.append(elem[gene_names])
    table_short_words.append(elem[organism])
    table_short_words.append(elem[length])
    table_short.append(table_short_words)

with open('table_prot_amyl_NO_PDB.tsv', 'w') as tab:
    for elem in table_short:
        for words in elem:
            tab.write(words + '\t')
        tab.write('\n')

with open('amount_of_motifs_No_PDB.tsv', 'w') as amount:
    for motif, count in counter_of_motifs.items():
        amount.write(motif + '\t')
        amount.write(str(count) + '\n')

with open('statistic_on_amount_of_motifs_in_every_peptide_NO_PDB.tsv', 'w') as statistic:
    statistic.write('Peptide' + '\t')
    statistic.write('Amount of motifs' + '\n')
    for peptide, count in amount_of_amyloids_motifs.items():
        statistic.write(peptide + '\t')
        statistic.write(str(count) + '\n')