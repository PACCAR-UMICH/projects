import random


def exercise_01(dna_sequence: str):
    """Calculating GC Content

    This function calculates the GC content of a given DNA sequence.
    The GC content is the percentage of nucleotides in the DNA sequence that are either guanine (G) or cytosine (C).

    Parameters:
    dna_sequence (str): The DNA sequence to calculate the GC content for.

    Returns:
    float: The GC content of the DNA sequence as a percentage.
    """
    gc_content = dna_sequence.count('G') + dna_sequence.count('C')
    total_length = len(dna_sequence)
    gc_percentage = (gc_content / total_length) * 100
    return gc_percentage


def exercise_02(dna_sequence: str, motif: str):
    """Finding Motif in DNA

    This function finds the positions of a motif within a DNA sequence.

    Parameters:
    dna_sequence (str): The DNA sequence to search within.
    motif (str): The motif to find within the DNA sequence.

    Returns:
    list: A list of positions where the motif occurs within the DNA sequence.
    """
    positions = []
    motif_length = len(motif)
    for i in range(len(dna_sequence) - motif_length + 1):
        if dna_sequence[i:i+motif_length] == motif:
            positions.append(i)
    return positions

def exercise_03(rna_sequence):
    codon_table = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y',
        'UAC': 'Y', 'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'UAA': 'Stop',
        'UAG': 'Stop', 'UGA': 'Stop', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L',
        'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R',
        'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I',
        'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
        'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V',
        'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A',
        'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein = ""
    print(rna_sequence)
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if codon_table[codon] == 'Stop':
            break
        protein += codon_table[codon]
    return protein

def exercise_04(dna_sequence1: str, dna_sequence2: str):
    """Calculating Hamming Distance between DNA sequences.
    This function calculates the Hamming distance between two DNA sequences.

    The Hamming distance is the number of positions at which the corresponding symbols are different.

    Parameters:
    dna_sequence1 (str): The first DNA sequence.
    dna_sequence2 (str): The second DNA sequence.

    Returns:
    int: The Hamming distance between the two DNA sequences.

    Note:
    The two DNA sequences must be of the same length.
    """
    distance = 0
    for el1, el2 in zip(dna_sequence1, dna_sequence2):
        if el1 != el2:
            distance += 1
    return distance

def exercise_05(dna_sequence: str, mutation_probability: float):
    """Simulates random mutations in a DNA sequence.

    This function simulates random mutations in a DNA sequence with a given probability.

    Parameters:
    dna_sequence (str): The DNA sequence.
    mutation_probability (float): The probability of a mutation at each position.

    Returns:
    str: The mutated DNA sequence.
    """
    nucleotides = ['A', 'T', 'C', 'G']
    mutated_sequence = ""
    for nucleotide in dna_sequence:
        if random.random() < mutation_probability:
            mutated_sequence += random.choice([nt for nt in nucleotides if nt != nucleotide])
        else:
            mutated_sequence += nucleotide
    return mutated_sequence


if __name__ == "__main__":
    # Testing exercise_01
    dna_sequence = "AGCTAGCTAGCTAGCTAGCT"
    gc_percentage = exercise_01(dna_sequence)
    print(f"The GC content is: {gc_percentage}%")
    # Expected output
    print("Expected output: 50%")
    print()

    # Testing exercise_02
    dna_sequence = "AGCTAGCTAGCTAGCTAGCT"
    motif = "AGCT"
    positions = exercise_02(dna_sequence, motif)
    print(f"The motif '{motif}' occurs at positions: {positions}")
    # Expected output
    print("Expected output: [0, 4, 8, 12, 16]")
    print()


    # Testing exercise_03
    rna_sequence = "AUGUUUUCUUACAAGUGGAUACCCACAGUGGACCAUGAAGGAUAA"
    protein_sequence = exercise_03(rna_sequence)
    print(f"The protein sequence is: {protein_sequence}")
    # Expected output
    print("Expected output: MFSYKWIPTVDHEG")
    print()

    # Testing exercise_04
    dna_sequence1 = "AGCTAGCTAGCTAGCTAGCT"
    dna_sequence2 = "AGCTAGCTAGGTAGCTAGCT"
    hamming_distance = exercise_04(dna_sequence1, dna_sequence2)
    print(f"The Hamming distance is: {hamming_distance}")
    # Expected output
    print("Expected output: 1")
    print()

    # Testing exercise_05
    dna_sequence = "ATCGAATCG"
    mutation_probability = 0.1
    mutated_sequence = exercise_05(dna_sequence, mutation_probability)
    print(f"The mutated DNA sequence is: {mutated_sequence}")
    # Example output (since it is random)
    print("Example output: ATCGAATGG")
    print()
