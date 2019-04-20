from feature_matrix import build_feature_matrix
from random import sample
import sys
import utils


def main():
    # Parse args and split data into training and testing
    pssm_list, tm_align_list, fasta_list, pssm_dir, tm_align_dir, fasta_dir = parse_args()
    pssm_train = sample(pssm_list, int(0.75 * len(pssm_list)))
    pssm_test = [pssm for pssm in pssm_list if pssm not in pssm_train]
    # train(pssm_train, pssm_dir, tm_align_dir)

    # TODO: For testing purposes, remove later
    feature_matrix = build_feature_matrix(pssm_train, pssm_dir, fasta_dir, tm_align_dir)


# TODO: Edit error message
err_msg = '''
Please enter two directory names (absolute paths)
containing sequences for linear regression training data
(with double quotes around them if they have spaces).
The directory with PSSM files should come first, 
followed by the path to the tmalign files.'''


def parse_args():
    if len(sys.argv) < 4:
        print(err_msg)
        sys.exit()

    try:
        # Get the lists of pssm and rr file names
        pssm = utils.read_directory_contents(sys.argv[1], '.pssm')
        tm_align = utils.read_directory_contents(sys.argv[2], '_tmalign')
        fasta = utils.read_directory_contents(sys.argv[3], '.fasta')
    except:
        # Given paths are not valid directories
        print(err_msg)
        sys.exit()

    # Return list of pssm & tm_align files, and their parent directories
    return pssm, tm_align, fasta, sys.argv[1], sys.argv[2], sys.argv[3]


if __name__ == '__main__':
    main()