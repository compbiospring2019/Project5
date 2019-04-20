from solvent_accessibility.amino_acids import options


class Node(object):
    parent = None           # Pointer to parent node
    children = []           # Pointers to children nodes (index is attribute value)
    attributes_left = []    # List of attribute (strings) left
    molecules = []          # List of amino acids (dicts) that fit into this subtree
    attribute = ''          # The attribute choice that this node represents
    rsa_value = None          # The rsa-value this node represents, if it is a leaf node

    def __init__(self, parent_node=None):
        self.parent = parent_node

    def get_outcomes(self, attribute, attr_value):
        if attribute not in self.attributes_left:
            raise Exception('Attribute {} not in {}'.format(attribute, self.attributes_left))

        # Go through molecules and get the rsa-label for that
        # value of the attribute
        outcomes = []
        for molecule in self.molecules:
            if molecule[attribute] == attr_value:
                outcomes.append(molecule['rsa-label'])
        return outcomes

    def make_children(self, attribute):
        # Create & return the children nodes of this node based on the current attribute
        if attribute not in self.attributes_left:
            raise Exception('Attribute {} not in {}'.format(attribute, self.attributes_left))

        # Set this node's attribute
        self.attribute = attribute

        # Create the children
        children = []
        # For each value of the attribute, set parent attributes left, and molecules
        # Note: index in array of children is the attribute's value
        for attr_value in range(options[attribute]):
            child = Node(self)
            child.attributes_left = [attr for attr in self.attributes_left if attr != attribute]
            child.molecules = [mol for mol in self.molecules if mol[attribute] == attr_value]
            children.append(child)

        self.children = children

        return children

    def calc_rsa_value(self):
        # This node is a leaf node. Calculate the RSA value for this node.
        if self.rsa_value is not None:
            return self.rsa_value

        rsa_values = [mol['rsa-label'] for mol in self.molecules]
        values_list = []
        max_val = None
        for value in range(options['rsa-label']):
            values_list.append(rsa_values.count(value))
            if max_val is None or values_list[value] > values_list[max_val]:
                max_val = value

        self.rsa_value = max_val

        return max_val
