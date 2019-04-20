# For each option, 0 means NO and 1 means YES
amino_acids = {
    'A': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 1,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'C': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'D': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 1,
        'positive': 0,
        'negative': 1,
        'small': 1,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'E': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 1,
        'positive': 0,
        'negative': 1,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'F': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 1,
        'proline': 0
    },
    'G': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 1,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'H': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 1,
        'positive': 1,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 1,
        'proline': 0
    },
    'I': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 1,
        'aromatic': 0,
        'proline': 0
    },
    'K': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 1,
        'positive': 1,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'L': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 1,
        'aromatic': 0,
        'proline': 0
    },
    'M': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'N': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'P': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 1
    },
    'Q': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'R': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 1,
        'positive': 1,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'S': {
        'hydrophobic': 0,
        'polar': 1,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 1,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'T': {
        'hydrophobic': 1,
        'polar': 1,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 0,
        'proline': 0
    },
    'V': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 1,
        'tiny': 0,
        'aliphatic': 1,
        'aromatic': 0,
        'proline': 0
    },
    'W': {
        'hydrophobic': 1,
        'polar': 0,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 1,
        'proline': 0
    },
    'Y': {
        'hydrophobic': 1,
        'polar': 1,
        'charged': 0,
        'positive': 0,
        'negative': 0,
        'small': 0,
        'tiny': 0,
        'aliphatic': 0,
        'aromatic': 1,
        'proline': 0
    }
}

options = {
    'hydrophobic': 2,
    'polar': 2,
    'charged': 2,
    'positive': 2,
    'negative': 2,
    'small': 2,
    'tiny': 2,
    'aliphatic': 2,
    'aromatic': 2,
    'proline': 2,
    'rsa-label': 2
}


def get_amino_acid(name):
    acid_dictionary = {'name': name}
    acid_dictionary.update(amino_acids[name])
    return acid_dictionary


rsa_labels = {
    'B': 0,
    'E': 1,
    0: 'B',
    1: 'E'
}
