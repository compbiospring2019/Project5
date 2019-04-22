from feature_matrix import build_feature_matrix
import sys
import utils


MODEL = None


def main():
    pssm_dir, pssm_files, fasta_dir = parse_args()

    feature_matrix = build_feature_matrix(pssm_files, pssm_dir, fasta_dir)

    tm_score = predict_tm_score(feature_matrix[0])

    print('\nCalculated TM-Score for {} and {} is {}'.format(pssm_files[0], pssm_files[1], tm_score))


def predict_tm_score(feature_row):
    global MODEL
    if MODEL is None:
        MODEL = utils.read_model()

    # Start with the intercept
    tm_score = MODEL.get('intercept', 0.0)

    # Get rid of the tm-score, if it's in the dictionary
    actual_tm_score = feature_row.pop('tm-score', None)

    # Compute the sum of weight vector values and feature values
    for seq_num in feature_row.keys():
        for feat_type in feature_row[seq_num].keys():
            for feat_name in feature_row[seq_num][feat_type].keys():
                tm_score += feature_row[seq_num][feat_type][feat_name] * \
                    MODEL[seq_num][feat_type][feat_name]

    return tm_score


err_msg = '''
Please enter the directory containing the two PSSM files
followed by the file names for the two PSSM files within
that directory (Note: just the file names, not directories).
The fourth argument should be the directory containing
the .fasta files corresponding to the two PSSM files.'''


def parse_args():
    if len(sys.argv) < 5:
        print(err_msg)
        sys.exit()

    # Return the directory and two PSSM file names
    return sys.argv[1], [sys.argv[2], sys.argv[3]], sys.argv[4]


if __name__ == '__main__':
    main()
