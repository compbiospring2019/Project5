# Utils for COMP 5970 Project 5

import os

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


def read_fasta(file_path, dir=None):
    """
    Reads a fasta sequence
    :return: A fasta sequence as a string, in upper case
    """
    if file_path is None:
        return None
    sequence = ''
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        for line in f:
            sequence += line.strip()

    return sequence.upper()


def read_directory_contents(path, file_extension):
    """
    Lists all files ending with the file_extension in the directory path
    :return: list of file names
    """
    if not os.path.isdir(path):
        # This is not a valid directory!
        raise Exception('Not a valid directory!')

    # Return a list of files with the file extension
    ls_dir = os.listdir(path)
    return [file_name for file_name in ls_dir if file_name.endswith(file_extension)]


def read_pssm(file_path, dir=None):
    """
    Reads in the right half of a PSSM file
    :return: list of dictionaries, where each dict is a row of the matrix.
        Keys are amino acids, values are the value in the row of the matrix
    """
    if dir:
        file_path = os.path.join(dir, file_path)
    pssm = []
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        if title in ['', '\n']:
            title = f.readline()

        # Get the list of amino acids on the top
        headers = f.readline().strip().split()[20:]

        # Now, read each line of the matrix into a dictionary
        for line in f:
            if line in ['', '\n']:
                break
            line_list = line.strip().split()
            row = {'this-acid': line_list[1]}
            line_list = line_list[20:-2]
            for acid_num in range(len(headers)):
                print('working on {}'.format(line_list[acid_num + 2]))
                row[headers[acid_num]] = int(line_list[acid_num + 2])
            pssm.append(row)

    # Returns a list of dictionaries, where each dict is a row of the matrix
    return pssm


def read_tmalign(file_path, dir=None):
    """
    Reads in a tmalign file
    :return: Dictionary with sequence names and TM scores
    """
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        # Skip the header
        line = f.readline()

        while not line.startswith('TM-score'):
            line = f.readline()

        tm_score_1 = get_tm_score_from_tmalign(line)
        line = f.readline()
        tm_score_2 = get_tm_score_from_tmalign(line)

    return (tm_score_1 + tm_score_2) / 2.0


def get_tm_score_from_tmalign(line):
    line_parts = line.split()
    return float(line_parts[1])
