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
        print(pssm_file)
        sequence_name = pssm_file.replace('.pssm', '')
        # Check if the PSSM has been read in yet
        if sequence_name not in sequence_info:
            # Read in sequence, add to sequence_info
            info = calc_sequence_info()
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


def calc_sequence_info():
    """
    Calculates information about a sequence
    :return: Dictionary - {SS:{E:val, ...}, SA:{E:val}, PSSM:{amino_acid: avg PSSM val}
    """
    # TODO: Calc PSSM averages

    # TODO: Calc SS percentages

    # TODO: Calc SA percentages

    return None
