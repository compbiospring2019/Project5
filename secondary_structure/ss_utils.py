# Utils for Project 3
import os
from random import sample

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


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


def compute_percentages(ss_prediction):
    label_c_len = float(len([x for x in ss_prediction if x == 'C']))
    label_e_len = float(len([x for x in ss_prediction if x == 'E']))
    percents = {
        'C': label_c_len / len(ss_prediction),
        'E': label_e_len / len(ss_prediction)
    }
    percents['H'] = 1.0 - percents['C'] - percents['E']
    return percents
