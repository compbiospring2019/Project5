# Utils for Project 2
import os
import json

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


# Read a biological sequence or RSA sequence from a file:
def read_sequence(file_path, dir=None):
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


def read_json():
    # Reads in the JSON object from model.json
    json_file = os.path.join(parent_directory, 'sa_model.json')
    with open(json_file, 'r') as file:
        contents = json.load(file)

    return contents


def compute_percentages(rsa_prediction):
    label_b_len = float(len([x for x in rsa_prediction if x == 'B']))
    percents = {
        'B': label_b_len / len(rsa_prediction)
    }
    percents['E'] = 1.0 - percents['B']
    return percents
