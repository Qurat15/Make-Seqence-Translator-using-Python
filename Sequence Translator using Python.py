def transcribe(dna_seq):
    # Replace 'T' with 'U' to transcribe DNA to RNA
    rna_seq = dna_seq.replace('T', 'U')
    return rna_seq

def translate_sequence(rna_seq):
    genetic_code = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGA': '*',
        'UGU': 'C', 'UGC': 'C', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    protein_seq = ""
    i = 0
    while i < len(rna_seq) - 2:
        codon = rna_seq[i:i + 3]
        amino_acid = genetic_code.get(codon, "X")  # "X" for unknown
        if amino_acid == '*':
            break
        protein_seq += amino_acid
        i += 3

    return protein_seq

def generate_reading_frames(seq):
    reading_frames = []

    # Iterate through all three possible reading frames
    for i in range(3):
        frame = seq[i:]
        reading_frames.append(frame)
        rna_seq = transcribe(frame)
        protein_seq = translate_sequence(rna_seq)
        print(f"Reading Frame {i + 1} (Starting from {frame[:3]}):")
        print("Transcribed RNA Sequence:", rna_seq)
        print("Translated Protein Sequence:", protein_seq)
        print()

    return reading_frames

input_sequence = input("Enter a DNA or RNA sequence (uppercase letters only): ")
print("Input Sequence:", input_sequence)
print()

rna_sequence = transcribe(input_sequence)
protein_sequence = translate_sequence(rna_sequence)
print("Transcribed RNA Sequence:", rna_sequence)
print("Translated Protein Sequence:", protein_sequence)
print()

reading_frames = generate_reading_frames(input_sequence)