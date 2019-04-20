# Utils for Project 3
import os
from random import sample

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


def read_sequence(file_path, dir=None):
    # Read a biological sequence or RSA sequence from a file
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
    if not os.path.isdir(path):
        # This is not a valid directory!
        raise Exception('Not a valid directory!')

    # Return a list of files with the file extension
    ls_dir = os.listdir(path)
    return [file_name for file_name in ls_dir if file_name.endswith(file_extension)]


def read_pssm(file_path, dir=None):
    if dir:
        file_path = os.path.join(dir, file_path)
    pssm = []
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        if title in ['', '\n']:
            title = f.readline()

        # Get the list of amino acids on the top
        headers = f.readline().strip().split()[:20]

        # Now, read each line of the matrix into a dictionary
        for line in f:
            if line in ['', '\n']:
                break
            line_list = line.strip().split()[:22]
            row = {'this-acid': line_list[1]}
            for acid_num in range(len(headers)):
                row[headers[acid_num]] = int(line_list[acid_num + 2])
            pssm.append(row)

    # Returns a list of dictionaries, where each dict is a row of the matrix
    return pssm


def read_dist(file_path, dir=None):
    # Prior - probability of class label
    # Dists - list of means and standard deviations
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        prior = float(f.readline())
        dists = [[float(x) for x in line] for line in f]
    return prior, dists


def split_files(pssm_list, ss_list):
    # Divide the .pssm files into training and testing sets
    test_correct_pssm_files(pssm_list, ss_list)
    pssm_train = sample(pssm_list, int(0.75 * len(pssm_list)))
    pssm_test = [pssm_name for pssm_name in pssm_list if pssm_name not in pssm_train]
    return pssm_train, pssm_test


def test_correct_pssm_files(pssm_list, ss_list):
    # Make sure each .pssm has a corresponding .ss
    for pssm_name in pssm_list:
        if pssm_name.replace('.pssm', '.ss') not in ss_list:
            raise Exception('PSSM files don\'t match up with .ss files: {}'.format(pssm_name))


def compute_percentages(ss_prediction):
    label_c_len = float(len([x for x in ss_prediction if x == 'C']))
    label_e_len = float(len([x for x in ss_prediction if x == 'E']))
    percents = {
        'C': label_c_len / len(ss_prediction),
        'E': label_e_len / len(ss_prediction)
    }
    percents['H'] = 1.0 - percents['C'] - percents['E']
    return percents
