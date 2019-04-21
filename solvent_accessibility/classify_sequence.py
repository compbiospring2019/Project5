import os
import solvent_accessibility.sa_utils as utils
from solvent_accessibility.decision_tree import DecisionTree


# Read in the model
json_model = utils.read_json()


def classify(test_file, test_file_dir=None):
    if test_file_dir:
        test_file = os.path.join(test_file_dir, test_file)
    test_sequence = utils.read_sequence(test_file)

    # Classify the test sequence
    predicted_rsa = DecisionTree.classify_sequence(test_sequence, json_model)

    return utils.compute_percentages(predicted_rsa)


if __name__ == '__main__':
    classify()
