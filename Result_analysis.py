from pyfaidx import Fasta

proteome_seq = {}
with Fasta('proteome_without_pdb.fasta') as genes:
    n = 0
    for record in genes:
        line = '>' + record.long_name
        proteome_seq[line] = genes[n][:]
        n += 1

with open('proteome_without_pdb.result') as result:
    peptide_pos = {}
    n = 0
    dict_file = {}
    for line in result:
        position_amyl = []
        if line.find('#', 0, 1) == 0:
            continue
        if line.find('>', 0, 1) == 0:
            key = line
            n += 1
            start = line.find('amyloid_seq=')+12
            amyl_list = line[start:].split(sep='|')
            for elem in amyl_list:
                i = 0
                alphabet = True
                while alphabet == True:
                    alphabet = elem[:i+1].isdigit()
                    amyloid_num = elem[:i]
                    i += 1
                position_amyl.append(amyloid_num)
            peptide_pos[key] = position_amyl
            lines = []
        else:
            lines.append(line.replace('\n', '').split(sep='\t'))
            dict_file[key] = lines

high_index = {}
for key, val in peptide_pos.items():
    n = 0
    for position in val:
        for key_table, val_table in dict_file.items():
            if key in key_table:
                for pos_num in val_table:
                    if pos_num[0] == position:
                        for i in range(0, 6):
                            if float(val_table[val_table.index(pos_num)+i][2]) >= 0.9 or float(val_table[val_table.index(pos_num)+i][3]) >= 0.9:
                                n += 1
                        if 5 <= n <= 6:
                            high_index[key_table] = val_table

proteome_fasta = {}
k = 0
for key, val in proteome_seq.items():
    for key_high, val_high in high_index.items():
        if key in key_high:
            k += 1
            proteome_fasta[key] = val


print(k)
with open('filter_result_0.9.txt', 'w') as result:
    for key, val in high_index.items():
        result.write(key + '\n')
        for lines in val:
            for i in range(len(lines)):
                result.write(lines[i] + '\t')
            result.write('\n')

with open('proteome_frame_filter_0.9.fasta', 'w') as frame:
    for key, val in proteome_fasta.items():
        frame.write(key + '\n')
        frame.write(str(val)+ '\n')
