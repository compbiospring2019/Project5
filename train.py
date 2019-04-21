from feature_matrix import build_feature_matrix
from math import exp, log
from random import normalvariate, sample
import utils

# Global variables
SAMPLE_SIZE = 250
STEP_SIZE = 0.001


def gradient_descent(pssm_train, pssm_dir, fasta_dir, tm_align_dir):
    # Build the feature matrix
    feature_matrix = build_feature_matrix(pssm_train, pssm_dir, fasta_dir, tm_align_dir)

    w_vector = new_w_vector(feature_matrix[0])
    gradient_vector = None
    count = 0

    print('Training the model...')
    while not reached_top(w_vector, gradient_vector):
        count += 1
        print('{}st loop!'.format(count))
        gradient_vector = calc_gradient(w_vector, feature_matrix)
        w_vector = update_w(w_vector, gradient_vector)

    # Save the model to the file
    utils.write_model(w_vector)


def calc_gradient(w_vector, matrix):
    """
    Calculates the gradient based on (SAMPLE_SIZE) training examples
    :return: gradient vector
    """
    gradient_vector = new_w_vector(matrix[0])
    training_data = sample(matrix, SAMPLE_SIZE)

    for training_example in training_data:
        for seq_num in w_vector.keys():
            if seq_num == 'intercept':
                w_vector['intercept'] += training_example['tm-score'] - calc_sum(w_vector, training_example)
                continue

            for feat_type in w_vector[seq_num].keys():
                for feat_name in w_vector[seq_num][feat_type].keys():
                    w_vector[seq_num][feat_type][feat_name] += \
                        (training_example['tm-score'] - calc_sum(w_vector, training_example)) * \
                        training_example[seq_num][feat_type][feat_name]

    return gradient_vector


def update_w(w_vector, gradient_vector):
    """
    Updates each w value in the w_vector
    :return: w_vector
    """
    for seq_num in w_vector.keys():
        if seq_num == 'intercept':
            w_vector['intercept'] += 2 * STEP_SIZE * gradient_vector['intercept']
            continue

        for feat_type in w_vector[seq_num].keys():
            for feat_name in w_vector[seq_num][feat_type].keys():
                w_vector[seq_num][feat_type][feat_name] += \
                    2 * STEP_SIZE * gradient_vector[seq_num][feat_type][feat_name]

    return w_vector


def reached_top(w_vector, gradient_vector):
    """
    Check if we've reached the top of the mountain
    :return: boolean
    """
    if not gradient_vector:
        return False
    # TODO: Fix this.
    if normalvariate(0, 5) < 15:
        return False
    print('Reached the top!')
    return True


def calc_max_conditional_likelihood(w_vector, feature_matrix):
    sum_mcl = 0.0

    for feature in feature_matrix:
        feature_sum = calc_sum(w_vector, feature)
        sum_mcl += feature['class'] * feature_sum - log(1 + exp(feature_sum))

    return sum_mcl


def new_w_vector(feature):
    """
    Make a new w vector.
    :return: Nested dictionary, mirroring a feature in the feature matrix
    """
    w_vector = {}

    # Mirroring feature organization, randomly generate starting w values
    for seq_num in feature.keys():
        if seq_num == 'tm-score':
            continue
        w_vector[seq_num] = {}

        for feat_type in feature[seq_num].keys():
            w_vector[seq_num][feat_type] = {}

            for feat_name in feature[seq_num][feat_type].keys():
                w_vector[seq_num][feat_type][feat_name] = normalvariate(0, 4)

    w_vector['intercept'] = normalvariate(0, 4)

    return w_vector


def new_zero_vector(feature):
    """
    Make a new vector initialized with all zeros.
    :return: Nested dictionary, mirroring a feature in the feature matrix
    """
    vector = {}

    # Mirroring feature organization, randomly generate starting w values
    for seq_num in feature.keys():
        if seq_num == 'tm-score':
            continue
        vector[seq_num] = {}

        for feat_type in feature[seq_num].keys():
            vector[seq_num][feat_type] = {}

            for feat_name in feature[seq_num][feat_type].keys():
                vector[seq_num][feat_type][feat_name] = 0.0

    vector['intercept'] = 0.0

    return vector


def calc_sum(w_vector, feature_vector):
    """
    Calculates the sum w0 + SUM(wi * xi)
    """
    sum_of_w = w_vector['intercept']

    for seq_num in w_vector.keys():
        if seq_num == 'intercept':
            continue

        for feat_type in w_vector[seq_num].keys():
            for feat_name in w_vector[seq_num][feat_type].keys():
                sum_of_w += w_vector[seq_num][feat_type][feat_name] * feature_vector[seq_num][feat_type][feat_name]

    return sum_of_w
