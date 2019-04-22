from feature_matrix import build_feature_matrix
from random import normalvariate, sample
import utils

# Global variables
SAMPLE_SIZE = 750
STEP_SIZE = 0.001
SMALLEST_SQUARED_ERROR = None
BEST_MODEL = None
ITERATIONS_SINCE_SMALLEST = 0


def gradient_descent(pssm_train, pssm_dir, fasta_dir, tm_align_dir):
    # Build the feature matrix
    print('Building feature matrix...')
    feature_matrix = build_feature_matrix(pssm_train, pssm_dir, fasta_dir, tm_align_dir)

    # Batch
    # global SAMPLE_SIZE
    # SAMPLE_SIZE = len(feature_matrix)

    w_vector = new_w_vector(feature_matrix[0])
    count = 0

    print('Training the model...')
    while not reached_top(w_vector, feature_matrix):
        count += 1
        gradient_vector = calc_gradient(w_vector, feature_matrix)
        w_vector = update_w(w_vector, gradient_vector)

    print('Gradient Descent completed in {} iterations'.format(count))
    print('with final step size {} and sample size {}'.format(STEP_SIZE, SAMPLE_SIZE))

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
    new_w = new_zero_vector(w_vector)
    for seq_num in w_vector.keys():
        if seq_num == 'intercept':
            new_w['intercept'] = w_vector['intercept'] + 2 * STEP_SIZE * gradient_vector['intercept']
            continue

        for feat_type in w_vector[seq_num].keys():
            for feat_name in w_vector[seq_num][feat_type].keys():
                new_w[seq_num][feat_type][feat_name] = \
                    w_vector[seq_num][feat_type][feat_name] + 2 * STEP_SIZE * gradient_vector[seq_num][feat_type][feat_name]

    return new_w


def reached_top(w_vector, feature_matrix):
    """
    Check if we've reached the top of the mountain
    :return: boolean
    """
    global SMALLEST_SQUARED_ERROR, BEST_MODEL, ITERATIONS_SINCE_SMALLEST, STEP_SIZE
    squared_error = calc_squared_error(w_vector, feature_matrix)

    if SMALLEST_SQUARED_ERROR is None or squared_error < SMALLEST_SQUARED_ERROR:
        # New lowest squared error found
        SMALLEST_SQUARED_ERROR = squared_error
        BEST_MODEL = feature_matrix
        ITERATIONS_SINCE_SMALLEST = 0
        return False

    ITERATIONS_SINCE_SMALLEST += 1

    if ITERATIONS_SINCE_SMALLEST == 10:
        STEP_SIZE /= 10

    if ITERATIONS_SINCE_SMALLEST > 50:
        # We haven't found a better mountaintop in a while, so we've probably reached it
        print('Reached the top!')
        return True
    return False


def calc_squared_error(w_vector, feature_matrix):
    sum_error = 0.0

    for feature in feature_matrix:
        difference = feature['tm-score'] - calc_sum(w_vector, feature)
        sum_error += difference ** 2

    return sum_error / len(feature_matrix)


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
                w_vector[seq_num][feat_type][feat_name] = normalvariate(0, 0.5)

    w_vector['intercept'] = normalvariate(0, 0.5)

    return w_vector


def new_zero_vector(feature):
    """
    Make a new vector initialized with all zeros.
    :return: Nested dictionary, mirroring a feature in the feature matrix
    """
    vector = {}

    # Mirroring feature organization, randomly generate starting w values
    for seq_num in feature.keys():
        if seq_num in ['tm-score', 'intercept']:
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
