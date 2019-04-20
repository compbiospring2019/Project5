from solvent_accessibility.amino_acids import get_amino_acid, rsa_labels, options
from solvent_accessibility.node import Node
import math
from random import sample
import solvent_accessibility.sa_utils as utils


class DecisionTree(object):
    root = Node()
    feature_matrix = []
    fasta_train = None
    fasta_test = None
    fasta_dir = ''
    sa_dir = ''

    def __init__(self, fasta_list, sa_list, fasta_dir, sa_dir):
        # Split files into training and testing data
        self.split_files(fasta_list, sa_list)

        # Save the directory paths
        self.fasta_dir = fasta_dir
        self.sa_dir = sa_dir

    def split_files(self, fasta_list, sa_list):
        # Split files into testing and training data
        self.test_correct_fasta_files(fasta_list, sa_list)
        self.fasta_train = sample(fasta_list, int(0.75 * len(fasta_list)))
        self.fasta_test = [f_name for f_name in fasta_list if f_name not in self.fasta_train]

    def test_correct_fasta_files(self, fasta_list, sa_list):
        # Make sure a .sa file exists for each .fasta file
        for fasta_name in fasta_list:
            if fasta_name.replace('.fasta', '.sa') not in sa_list:
                raise Exception('FASTA files don\'t match up with .sa files: {}'.format(fasta_name))

    def build_feature_matrix(self):
        # For each fasta file in the training data, read the sequences and add them to the feature matrix
        for fasta_name in self.fasta_train:
            fasta = utils.read_sequence(fasta_name, self.fasta_dir)
            sa = utils.read_sequence(fasta_name.replace('.fasta', '.sa'), self.sa_dir)
            for index in range(len(fasta)):
                # Create the AA object
                acid = get_amino_acid(fasta[index].upper())

                # Add the RSA label
                rsa_binary = rsa_labels[sa[index]]
                acid['rsa-label'] = rsa_binary

                # Add the acid to the matrix
                self.feature_matrix.append(acid)

    def build_tree(self, current_node=None):
        # Recursive 'build the tree' funct starting with root node
        if current_node is None:
            # Start the tree with the root node
            self.root.attributes_left = [key for key in options.keys() if key not in ['name', 'rsa-label']]
            self.root.molecules = self.feature_matrix
            current_node = self.root

        # Choose best attribute
        best_attribute = self.compute_best_attribute(current_node)
        if best_attribute is None:
            # Data is perfectly classified, or there are no more attributes
            return

        # Create children nodes based on that attribute
        children = current_node.make_children(best_attribute)

        # Recursive calls to those children nodes
        for child in children:
            if len(child.molecules):
                self.build_tree(child)

    def compute_best_attribute(self, current_node):
        # Based on current node, select the best attribute from the leftover attributes
        total_outcomes = [mol['rsa-label'] for mol in current_node.molecules]
        total_entropy = self.compute_entropy(total_outcomes)

        if total_entropy == 0:
            # Data is perfectly classified
            return None

        attribute_gain = {}
        max_gain_attr = None
        for attr in current_node.attributes_left:
            attribute_gain[attr] = total_entropy
            for attr_val in range(options[attr]):
                # Calculate the portion of the weighted average
                # So, (Probability of attr=attr_val) * (Entropy(attr=attr_val))
                relevant_molecules = [mol for mol in current_node.molecules if mol[attr] == attr_val]
                probability = float(len(relevant_molecules)) / len(current_node.molecules)
                relevant_rsas = [mol['rsa-label'] for mol in relevant_molecules]
                entropy = self.compute_entropy(relevant_rsas)
                attribute_gain[attr] -= probability * entropy
            # Keep track of the attribute with the maximum gain
            if not max_gain_attr:
                max_gain_attr = attr
            elif attribute_gain[max_gain_attr] < attribute_gain[attr]:
                max_gain_attr = attr

        return max_gain_attr

    def compute_entropy(self, outcome_list):
        # Outcome list: the list of outcomes (y/n) for a given data set
        # or attribute

        # Calc number of each outcome
        outcome_count = {}
        for outcome in outcome_list:
            if outcome not in outcome_count:
                outcome_count[outcome] = 0
            outcome_count[outcome] += 1

        # Calculate the entropy (sum of p*log(p))
        sum = 0
        for outcome in outcome_count.keys():
            probability = float(outcome_count[outcome]) / len(outcome_list)
            sum -= probability * math.log(probability, 2)

        return sum

    def evaluate_model(self):
        # Keeping track of metrics (true condition, predicted condition)
        metrics = {
            (0, 0): 0,  # True negative
            (0, 1): 0,  # False positive
            (1, 0): 0,  # False negative
            (1, 1): 0   # True positive
        }
        # For each fasta file in the testing data, walk the tree for each amino acid
        for fasta_name in self.fasta_test:
            fasta = utils.read_sequence(fasta_name, self.fasta_dir)
            sa = utils.read_sequence(fasta_name.replace('.fasta', '.sa'), self.sa_dir)
            for index in range(len(fasta)):
                # Test this amino acid against our decision tree
                amino_acid = fasta[index]
                expected_result = rsa_labels[sa[index]]
                calculated_result = self.walk_tree(get_amino_acid(amino_acid))
                # print('Acid {}, expected {}, calculated {}'.format(amino_acid, expected_result, calculated_result))
                metrics[expected_result, calculated_result] += 1

        self.calculate_eval_metrics(metrics)

    def walk_tree(self, feature_vector):
        # print('Walking the tree...')
        # print('Feature vector: {}'.format(feature_vector))
        current_node = self.root
        while current_node.children:
            attr_value = feature_vector[current_node.attribute]
            # print('{} and {}'.format(current_node.attribute, attr_value))
            current_node = current_node.children[attr_value]

        # print('Current node\'s rsa-label = {}\n'.format(current_node.molecules[0]['rsa-label']))

        return current_node.calc_rsa_value()

    def calculate_eval_metrics(self, metrics):
        metrics = self.clarify_metrics(metrics)
        print('True positive count:  {} ({:2.4}%)'.format(metrics['tp'], (float(metrics['tp']) * 100 /metrics['sum'])))
        print('True negative count:  {} ({:2.4}%)'.format(metrics['tn'], (float(metrics['tn']) * 100 /metrics['sum'])))
        print('False positive count: {} ({:2.4}%)'.format(metrics['fp'], (float(metrics['fp']) * 100 /metrics['sum'])))
        print('False negative count: {} ({:2.4}%)'.format(metrics['fn'], (float(metrics['fn']) * 100 /metrics['sum'])))

        precision = float(metrics['tp']) / (metrics['tp'] + metrics['fp'])
        recall = float(metrics['tp']) / (metrics['tp'] + metrics['fn'])
        accuracy = float(metrics['tp'] + metrics['tn']) / (metrics['tp'] + metrics['tn'] + metrics['fp'] + metrics['fn'])
        f1 = float(2 * precision * recall) / (precision + recall)
        mcc = float((metrics['tp'] * metrics['tn']) - (metrics['fp'] * metrics['fn']))
        denominator = math.sqrt((metrics['tp'] + metrics['fp']) * (metrics['tp'] + metrics['fn']) *
                         (metrics['tn'] + metrics['fp']) * (metrics['tn'] + metrics['fn']))
        if denominator == 0:
            # Note: this probably won't ever happen, but just in case...
            mcc = 'Divide by zero'
        else:
            mcc /= denominator

        print('Precision: {}'.format(precision))
        print('Recall:    {}'.format(recall))
        print('Accuracy:  {}'.format(accuracy))
        print('F1:        {}'.format(f1))
        print('MCC:       {}'.format(mcc))

    def clarify_metrics(self, metrics):
        # Switch the metrics dictionary into more human-readable keys
        metrics['tn'] = metrics.pop((0, 0), 0)
        metrics['fp'] = metrics.pop((0, 1), 0)
        metrics['fn'] = metrics.pop((1, 0), 0)
        metrics['tp'] = metrics.pop((1, 1), 0)

        sum = 0
        for key in metrics.keys():
            sum += metrics[key]

        metrics['sum'] = sum

        return metrics

    def to_json(self, current_node=None):
        # Convert the DT object into a JSON object (python dictionary)
        if current_node is None:
            # Start with the root node
            current_node = self.root

        if current_node.attribute is '' and not current_node.children:
            # This may be a leaf node
            if current_node.rsa_value is not None:
                # This is a leaf node
                return {'rsa-label': current_node.rsa_value}
            # This was not a leaf node
            return None

        json_node = {
            'attribute': current_node.attribute,
            'children': []
        }

        for child in current_node.children:
            # Create the json objects for the children nodes
            json_node['children'].append(self.to_json(child))

        return json_node

    @classmethod
    def classify_sequence(cls, sequence, json_object):
        # Classify the amino acids in the sequence via the json tree
        rsa_binary_list = []
        for molecule in sequence:
            rsa_value = cls.walk_json(get_amino_acid(molecule), json_object)
            rsa_binary_list.append(rsa_value)

        # Convert the binary labels back into E/B labels
        rsa_list = [rsa_labels[val] for val in rsa_binary_list]
        return ''.join(rsa_list)

    @classmethod
    def walk_json(cls, feature_vector, json_tree):
        # Walk through the tree and find the RSA label for the given feature vector
        current_node = json_tree

        while 'attribute' in current_node:
            # For each non-leaf node, grab the appropriate child node
            attr_value = feature_vector[current_node['attribute']]
            current_node = current_node['children'][attr_value]

        return current_node['rsa-label']
