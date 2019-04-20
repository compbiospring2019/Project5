import os
import sys
import solvent_accessibility.sa_utils as utils
from solvent_accessibility.decision_tree import DecisionTree


err_msg = '''
Please enter two directory names (absolute paths)
containing sequences for Decision Tree training data
(with double quotes around them if they have spaces).
The directory with FASTA files should come first, 
followed by the path to the .sa files.'''


def parse_args():
    if len(sys.argv) < 2:
        print(err_msg)
        sys.exit()

    if len(sys.argv) == 2:
        return sys.argv[1], None
    return sys.argv[1], sys.argv[2]


def classify(test_file, test_file_dir=None):
    if test_file_dir:
        test_file = os.path.join(test_file_dir, test_file)
    test_sequence = utils.read_sequence(test_file)

    # Read in the model
    json_model = utils.read_json()

    # Classify the test sequence
    predicted_rsa = DecisionTree.classify_sequence(test_sequence, json_model)

    return utils.compute_percentages(predicted_rsa)


if __name__ == '__main__':
    classify()
