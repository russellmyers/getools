
def read_fasta(file_name):
    txt = read_txt_file(file_name)
    lines = txt.split('\n')
    sequences = []
    sequence = ''
    fasta_header = ''
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '>':
            if len(sequence) > 0:
                sequences.append((fasta_header,sequence))
            fasta_header = line
            sequence = ''
        else:
            sequence += line

    if len(sequence) > 0:
       sequences.append((fasta_header, sequence))

    return sequences

def list_to_file(file_name, ar):
    with open(file_name, 'w') as f:
        f.write('\n'.join([str(x) for x in ar]))

def file_to_list(file_name):
    with open(file_name) as f:
        txt = f.read()
        return txt.split('\n')

def string_to_file(file_name, s):
    with open(file_name, 'w') as f:
        f.write(s)

def read_txt_file(file_name):
    with open(file_name) as f:
        txt = f.read()
    return txt