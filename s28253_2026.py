# Album number: s28253
# Date: 06/05/2026
# Description: Random DNA sequence generator in FASTA format
# Author: Meryem Orhan

import random


def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length."""
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence


def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary of sequence statistics.
    Keys: 'A', 'C', 'G', 'T' (float values, %),
          'GC' (float value, %)."""
    length = len(sequence)
    stats = {
        'A': round(sequence.count('A') / length * 100, 2),
        'C': round(sequence.count('C') / length * 100, 2),
        'G': round(sequence.count('G') / length * 100, 2),
        'T': round(sequence.count('T') / length * 100, 2),
    }
    stats['GC'] = round(stats['G'] + stats['C'], 2)
    return stats


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Gets an integer from the user in a range.
    In case of an error, repeats the question."""
    while True:
        user_input = input(prompt)
        try:
            value = int(user_input)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def insert_name(sequence: str, name: str) -> str:
    """Inserts a name at a random position in the sequence.
    Name written in lowercase letters."""
    name_lower = name.lower()
    position = random.randint(0, len(sequence))
    return sequence[:position] + name_lower + sequence[position:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a formatted FASTA record as a string."""
    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"
    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]
    return header + "\n" + "\n".join(lines) + "\n"


def find_motif(sequence: str, motif: str) -> list:
    """Searches for all occurrences of a motif in the sequence.
    Returns a list of positions (1-based, biological convention)."""
    positions = []
    motif = motif.upper()
    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            positions.append(i + 1)
    return positions


def get_complement(sequence: str) -> str:
    """Returns the complementary DNA strand (A<->T, C<->G)."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_map[nuc] for nuc in sequence)


def get_reverse_complement(sequence: str) -> str:
    """Returns the reverse complementary DNA strand."""
    return get_complement(sequence)[::-1]


def transcribe(sequence: str) -> str:
    """Transcribes a DNA sequence to mRNA (T -> U replacement)."""
    return sequence.replace('T', 'U')


def main():
    """Main function that orchestrates the program flow."""
    length = validate_positive_int("Enter sequence length: ")

    seq_id = input("Enter sequence ID: ")
    while ' ' in seq_id:
        print("Error: ID cannot contain whitespace.")
        seq_id = input("Enter sequence ID: ")

    description = input("Enter a description of the sequence: ")
    name = input("Enter your name: ")

    sequence = generate_sequence(length)
    sequence_with_name = insert_name(sequence, name)
    stats = calculate_stats(sequence)

    fasta_content = format_fasta(seq_id, description, sequence_with_name)
    filename = f"{seq_id}.fasta"

    with open(filename, 'w') as f:
        f.write(fasta_content)

    print(f"\nSequence saved to file: {filename}")
    print(f"Sequence statistics (n={length}):")
    print(f"  A: {stats['A']:.2f}%")
    print(f"  C: {stats['C']:.2f}%")
    print(f"  G: {stats['G']:.2f}%")
    print(f"  T: {stats['T']:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")

    motif = input("\nEnter a motif to search (or press Enter to skip): ")
    if motif:
        positions = find_motif(sequence, motif)
        if positions:
            print(f"Motif '{motif.upper()}' found at {len(positions)} position(s): {positions}")
        else:
            print(f"Motif '{motif.upper()}' not found in the sequence.")

    comp = get_complement(sequence)
    rev_comp = get_reverse_complement(sequence)

    with open(filename, 'a') as f:
        f.write(format_fasta(f"{seq_id}_complement", "Complementary strand", comp))
        f.write(format_fasta(f"{seq_id}_revcomp", "Reverse complementary strand", rev_comp))

    print(f"\nComplement and reverse complement added to {filename}")

    mrna = transcribe(sequence)
    with open(filename, 'a') as f:
        f.write(format_fasta(f"{seq_id}_mRNA", "mRNA transcript", mrna))

    print(f"mRNA transcript added to {filename}")


if __name__ == "__main__":
    main()