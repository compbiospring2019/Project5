from classify import predict_tm_score
from feature_matrix import build_feature_matrix


def test_model(pssm_test, pssm_dir, fasta_dir, tm_align_dir):
    # Build the feature matrix
    feature_matrix = build_feature_matrix(pssm_test, pssm_dir, fasta_dir, tm_align_dir)

    # Calculate squared errors
    squared_errors = [calc_squared_error(feature) for feature in feature_matrix]

    average_error = sum(squared_errors) / len(squared_errors)

    print('\nAverage of squared errors on the test set: {}'.format(average_error))
    print('Maximum squared error: {}'.format(max(squared_errors)))
    print('Minimum squared error: {}'.format(min(squared_errors)))


def calc_squared_error(feature):
    actual_tm_score = feature.pop('tm-score', 0.0)
    predicted_tm_score = predict_tm_score(feature)

    return (actual_tm_score - predicted_tm_score) ** 2
