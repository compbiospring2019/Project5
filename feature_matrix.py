from solvent_accessibility.classify_sequence import classify as classify_sa
from secondary_structure.classify_ss import classify as classify_ss
import utils


def build_feature_matrix(pssm_files, pssm_dir, fasta_dir, tm_align_dir=None):
    """
    Builds a feature matrix from the PSSM file data
    :return: A list of dictionaries
    """
    matrix = []
    # Store data that has already been read in
    # sequence_info: {sequence_name: {SS:{E:val, ...}, SA:{E:val}, PSSM:{amino_acid: avg PSSM val}}
    sequence_info = {}

    for pssm_file in pssm_files:
        sequence_name = pssm_file.replace('.pssm', '')

        # If a sequence hasn't been read in yet, read in sequence, and add to sequence_info
        if sequence_name not in sequence_info:
            info = calc_sequence_info(pssm_file, pssm_dir, fasta_dir)
            sequence_info[sequence_name] = info

        # Loop through pairs of sequences
        for other_pssm in pssm_files:
            if other_pssm == pssm_file:
                break
            other_seq_name = other_pssm.replace('.pssm', '')
            feature = {}
            # TODO: Add in features from sequence_info
            if tm_align_dir:
                # Read in TM score, add to feature
                file_name = '{}_{}_tmalign'.format(sequence_name, other_seq_name)
                feature['tm-score'] = utils.read_tmalign(file_name, tm_align_dir)

            matrix.append(feature)
    return matrix


def calc_sequence_info(pssm_file, pssm_dir, fasta_dir):
    """
    Calculates information about a sequence
    :return: Dictionary - {SS:{E:val, ...}, SA:{E:val}, PSSM:{amino_acid: avg PSSM val}
    """
    info = {}

    # Calc PSSM averages
    info['pssm'] = calc_pssm_averages(pssm_file, pssm_dir)

    info['ss'] = classify_ss(utils.read_pssm(pssm_file, pssm_dir))

    info['sa'] = classify_sa(pssm_file.replace('.pssm', '.fasta'), fasta_dir)

    return info


def calc_pssm_averages(pssm_file, pssm_dir):
    pssm = utils.read_pssm(pssm_file, pssm_dir)
    pssm_averages = {}
    for row in pssm:
        for key in row.keys():
            # For each amino acid (key), calculate the sum of the entire column
            if key != 'this-acid':
                if key not in pssm_averages:
                    pssm_averages[key] = 0.0
                pssm_averages[key] += row[key]

    # Scale down numbers
    for key in pssm_averages.keys():
        pssm_averages[key] /= (100 * len(pssm))

    return pssm_averages
