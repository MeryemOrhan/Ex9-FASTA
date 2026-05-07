# Album number: s28253
# Date: 06/05/2026
# Description: Random DNA sequence generator in FASTA format
# Author: Meryem Orhan

import random
import csv


def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length."""
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence


def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary of sequence statistics.
    Keys: 'A', 'C', 'G', 'T' (float values, %),
          'GC' (float value, %), 'gc_ratio_A' (float value, %)."""
    length = len(sequence)
    stats = {
        'A': round(sequence.count('A') / length * 100, 2),
        'C': round(sequence.count('C') / length * 100, 2),
        'G': round(sequence.count('G') / length * 100, 2),
        'T': round(sequence.count('T') / length * 100, 2),
    }
    stats['GC'] = round(stats['G'] + stats['C'], 2)
    stats['gc_ratio_A'] = stats['GC']
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


def translate(sequence: str) -> str:
    """Translates a DNA sequence into an amino acid sequence using the standard codon table.
    Stops translation at the first stop codon (*).
    Reference: https://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html"""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T','ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3].upper()
        aa = codon_table.get(codon, '?')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def sliding_window_gc(sequence: str, window_size: int) -> list:
    """Calculates GC content in a sliding window along the sequence.
    Returns a list of (start_position, gc_content) tuples.
    Positions are 1-based (biological convention)."""
    results = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        gc = round((window.count('G') + window.count('C')) / window_size * 100, 2)
        results.append((i + 1, gc))
    return results


def save_sliding_window_csv(results: list, filename: str) -> None:
    """Saves sliding window GC content results to a CSV file."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['start_position', 'gc_content'])
        writer.writerows(results)


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

    protein = translate(sequence)

    filename = f"{seq_id}.fasta"

    with open(filename, 'w') as f:
        f.write(format_fasta(seq_id, description, sequence_with_name))
        f.write(format_fasta(f"{seq_id}_complement", "Complementary strand", get_complement(sequence)))
        f.write(format_fasta(f"{seq_id}_revcomp", "Reverse complementary strand", get_reverse_complement(sequence)))
        f.write(format_fasta(f"{seq_id}_mRNA", "mRNA transcript", transcribe(sequence)))
        f.write(format_fasta(f"{seq_id}_protein", "Translated protein sequence", protein))
        f.write("# EOF_1\n")

    print(f"\nSequence saved to file: {filename}")
    print(f"Sequence statistics (n={length}):")
    print(f"  A: {stats['A']:.2f}%")
    print(f"  C: {stats['C']:.2f}%")
    print(f"  G: {stats['G']:.2f}%")
    print(f"  T: {stats['T']:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")
    print(f"\nProtein sequence (length={len(protein)}): {protein[:60]}{'...' if len(protein) > 60 else ''}")

    motif = input("\nEnter a motif to search (or press Enter to skip): ")
    if motif:
        positions = find_motif(sequence, motif)
        if positions:
            print(f"Motif '{motif.upper()}' found at {len(positions)} position(s): {positions}")
        else:
            print(f"Motif '{motif.upper()}' not found in the sequence.")

    window_size = validate_positive_int("\nEnter sliding window size: ", min_val=1, max_val=length)
    sw_results = sliding_window_gc(sequence, window_size)
    csv_filename = f"{seq_id}_gc_sliding.csv"
    save_sliding_window_csv(sw_results, csv_filename)
    print(f"Sliding window GC analysis saved to {csv_filename}")


if __name__ == "__main__":
    main()