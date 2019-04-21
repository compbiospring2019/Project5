from feature_matrix import build_feature_matrix
from math import exp, log
from random import normalvariate, sample
import utils

# Global variables
SAMPLE_SIZE = 10
STEP_SIZE = 0.001


def gradient_descent(pssm_train, pssm_dir, fasta_dir, tm_align_dir):
    # Build the feature matrix
    feature_matrix = build_feature_matrix(pssm_train, pssm_dir, fasta_dir, tm_align_dir)

    w_vector = new_w_vector(feature_matrix[0])
    gradient_vector = None

    print('Training the model...')
    while not reached_top(w_vector, gradient_vector):
        gradient_vector = calc_gradient(w_vector, feature_matrix)
        w_vector = update_w(w_vector, gradient_vector)

    # Save the model to the file
    utils.write_model(w_vector)


def calc_gradient(w_vector, matrix):
    """
    Calculates the gradient based on (SAMPLE_SIZE) training examples
    :return: gradient vector
    """
    gradient_vector = [0] * len(w_vector)
    training_data = sample(matrix, SAMPLE_SIZE)

    for training_example in training_data:
        # Calculate P(Y=1|X,w)
        p_hat = 1.0 / (1 + exp(calc_sum(w_vector, training_example)))

        # Deal with w0
        gradient_vector[0] += training_example['class'] - p_hat

        # For each feature, calculate the gradient
        for i in range(1, len(w_vector)):
            gradient_vector[i] += training_example[i-1] * (training_example['class'] - p_hat)

    return gradient_vector


def update_w(w_vector, gradient_vector):
    """
    Updates each w value in the w_vector
    :return: w_vector
    """
    for index in range(len(w_vector)):
        w_vector[index] += STEP_SIZE * gradient_vector[index]
    return w_vector


def reached_top(w_vector, gradient_vector):
    """
    Check if we've reached the top of the mountain
    :return: boolean
    """
    if not gradient_vector:
        return False
    for i in range(len(gradient_vector)):
        if gradient_vector[i] > 0.005:
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

    return w_vector


def calc_sum(w_vector, feature_vector):
    """
    Calculates the sum w0 + SUM(wi * xi)
    """
    sum_of_w = w_vector[0]
    for i in range(1, len(w_vector)):
        sum_of_w += w_vector[i] * feature_vector[i - 1]

    return sum_of_w
