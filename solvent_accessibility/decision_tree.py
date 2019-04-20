from solvent_accessibility.amino_acids import get_amino_acid, rsa_labels


class DecisionTree(object):
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
