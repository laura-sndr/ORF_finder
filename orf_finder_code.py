# BAM WS23/24 ORF Finder 3

# help with biopython and ideas: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/00%20-%20Tutorial%20-%20Index.html

# Import necessary modules for DNA sequence handling, file operations, and regular expressions
from Bio import SeqIO                                                               # Biopython module for sequence input/output
from Bio.Seq import Seq                                                             # Biopython module for sequence manipulation
import argparse                                                                     # Module for handling command-line arguments
import csv                                                                          # Module for handling CSV files
import re                                                                           # Module for working with regular expressions

# Function to ensure the DNA sequence length is a multiple of 3
def pad_sequence(sequence):
    while len(sequence) % 3 != 0:  # While loop ensures the sequence length is a multiple of 3
        sequence += 'N'  # Add 'N' nucleotide to the sequence to pad it
    return sequence  # Return the padded sequence

# Function to find positions of regulatory elements within the DNA sequence
def find_regulatory_elements(sequence_str, regulatory_element):
    regulatory_positions = []  # Initialize an empty list to store regulatory element positions
    regulatory_pattern = re.compile(fr'({regulatory_element})', re.IGNORECASE)  # Compile a regex pattern to find the regulatory element (more efficient for repeated use)
    #regulatory_pattern: regular expression pattern is compiled using re.compile, pattern is designed to match the specified regulatory_element (ignoring case)
    #f-string: creates a regular expression pattern by inserting the value of regulatory_element into the pattern
    #re.IGNORECASE: specifies a flag that makes the regular expression pattern case-insensitive
    regulatory_matches = regulatory_pattern.finditer(sequence_str)  # Find all matches of the regulatory element in the sequence

    # Iterate over the matches and store their positions
    for match in regulatory_matches:
        regulatory_positions.append(match.start() + 1)  # Store the start position of each regulatory element match, #match.start() + 1: start() method of a match object returns the starting position of the match. Adding 1 is done to convert from 0-based indexing to 1-based indexing

    return regulatory_positions  # Return the list of regulatory element positions

# Function to find start and stop codons within the DNA sequence
def find_start_and_stop_codons(sequence_str, start_codons, stop_codons):
    start_positions = []  # Initialize an empty list to store start codon positions
    stop_positions = []  # Initialize an empty list to store stop codon positions

    # Iterate over the DNA sequence in all 3 reading frames
    for frame in range(3):
        seq_frame = sequence_str[frame:]  # Extract the sequence in the current reading frame
        for i in range(0, len(seq_frame), 3):  # Iterate over the sequence in steps of 3 (codon length)
            codon = seq_frame[i:i + 3]  # Extract a codon (sequence of 3 nucleotides)
            if str(codon) in start_codons:  # Check if the codon is a start codon
                start_pos = i + frame  # Calculate the start position of the codon
                start_positions.append(start_pos)  # Store the start codon position
            elif str(codon) in stop_codons:  # Check if the codon is a stop codon
                stop_pos = i + frame  # Calculate the stop position of the codon
                stop_positions.append(stop_pos)  # Store the stop codon position

    return start_positions, stop_positions  # Return lists of start and stop codon positions

# Function to translate a DNA sequence into an amino acid sequence based on the organism type
def translate_sequence(sequence, organism_type):
    if organism_type.lower() == "prokaryotic":  # Check the type of organism
        return Seq(sequence).translate(table=11)  # Translate the sequence using the prokaryotic genetic code (table 11)
    else:
        return Seq(sequence).translate()  # Translate the sequence using the default genetic code

#function processes an open reading frame (ORF) by extracting its nucleotide sequence, padding it to ensure it's a multiple of 3, generating a string representing the frame, and finally translating it into an amino acid sequence
def process_open_reading_frame(start_pos, stop_pos, sequence_str, organism_type, regulatory_positions):
    orf_sequence = sequence_str[start_pos:stop_pos]  # Extract the ORF sequence
    orf_sequence = pad_sequence(orf_sequence)  # Pad the sequence to ensure it's a multiple of 3

    orf_frame = f"{start_pos + 1}*{stop_pos}*+"  # Generate a string representing the ORF frame, e.g: generated orf_frame string might look like this: "101*200*+", where 101 is the start position, 200 is the stop position and the + gives additional information about the ORF
    orf_sequence_aa = translate_sequence(orf_sequence, organism_type)  # Translate the ORF sequence to amino acids

    return orf_frame, str(orf_sequence), str(orf_sequence_aa)  # Return the ORF frame, nucleotide sequence, and translated amino acid sequence

# Function to find open reading frames (ORFs) in the DNA sequence
# with help from ChatGPT
def find_open_reading_frames(sequence, organism_type, regulatory_positions):
    start_codons = ["ATG"]  # List of start codons
    stop_codons = ["TAA", "TGA", "TAG"]  # List of stop codons
    sequence_str = str(sequence)  # Convert the DNA sequence object to a string
    start_positions, stop_positions = find_start_and_stop_codons(sequence_str, start_codons, stop_codons)  # Find start and stop codon positions

    orfs = []  # Initialize a list to store identified ORFs

    # Iterate over potential start positions to identify ORFs
    for start_pos in start_positions:
        closest_regulatory_position = "Not found"  # Initialize variables to track regulatory element distance
        closest_regulatory_distance = "Not found"

        # Iterate over regulatory positions to find the closest regulatory element
        for regulatory_pos in regulatory_positions:
            if 0 < start_pos - regulatory_pos <= 500:  # Check if the distance between start codon and regulatory element is within a threshold
                distance = abs(start_pos - regulatory_pos)  # Calculate the distance
                if closest_regulatory_distance == "Not found" or (0 < distance <= 500 and distance < closest_regulatory_distance):
                    closest_regulatory_distance = distance  # Update the closest distance
                    closest_regulatory_position = regulatory_pos + 1  # Update the closest regulatory position

        # Iterate over potential stop positions to identify ORFs
        for stop_pos in stop_positions:
            if stop_pos > start_pos and stop_pos - start_pos >= 50:  # Check if the stop codon is after the start codon and the ORF length is above a threshold
                # Process the ORF and get its frame, nucleotide sequence, and translated amino acid sequence
                orf_frame, orf_sequence, orf_sequence_aa = process_open_reading_frame(start_pos, stop_pos, sequence_str, organism_type, regulatory_positions)

                nested_orf = "No"  # Assume the ORF is not nested initially

                # Check if the current ORF is nested within any existing ORF
                for existing_orf in orfs:
                    if existing_orf[3] < start_pos < existing_orf[4] or existing_orf[3] < stop_pos < existing_orf[4] or start_pos < existing_orf[3] and stop_pos > existing_orf[4]:
                        nested_orf = "Yes"  # Update the nested status if the ORF is nested within another ORF
                        break  # No need to check further

                # Append the details of the identified ORF to the list of ORFs
                orfs.append((orf_frame, orf_sequence, orf_sequence_aa, start_pos, stop_pos, closest_regulatory_distance, closest_regulatory_position, nested_orf))
                break  # Move to the next potential ORF

    return orfs  # Return the list of identified ORFs

#help with argparse: https://realpython.com/command-line-interfaces-python-argparse/

# Function to handle the main execution flow
def main():
    parser = argparse.ArgumentParser(description='Find open reading frames and regulatory elements in a DNA sequence.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output CSV file')
    parser.add_argument('organism_type', choices=['eukaryotic', 'prokaryotic', 'viral'], help='Type of organism')

    args = parser.parse_args()  # Parse the command line arguments

    sequences = list(SeqIO.parse(args.input_file, "fasta"))  # Read sequences from the input FASTA file

    if len(sequences) != 1:  # Ensure there is exactly one DNA sequence in the file
        print("Please provide a file with only one DNA sequence.")
        return  # Exit if more or fewer than one sequence is found

    sequence = str(sequences[0].seq.upper())  # Extract and convert the DNA sequence to uppercase

    # Define regulatory elements for different organism types
    #help and ideas: Cookbook from https://biopython-tutorial.readthedocs.io/en/latest/notebooks/19%20-%20Cookbook%20-%20Cool%20things%20to%20do%20with%20it.html
    regulatory_element_prokaryotic = "TA[AT][AT][AT]T"
    regulatory_element_eukaryotic = "TATA[AT]A[AT]"

    # Find regulatory element positions based on the organism type
    regulatory_positions = find_regulatory_elements(sequence, regulatory_element_prokaryotic if args.organism_type.lower() == "prokaryotic" else regulatory_element_eukaryotic)

    #help with csv: https://docs.python.org/3/library/csv.html
    with open(args.output_file, 'w', newline='') as csvfile:  # Open the output CSV file for writing
        writer = csv.writer(csvfile)  # Create a CSV writer

        # Write appropriate headers based on the organism type
        if args.organism_type.lower() == "prokaryotic":
            writer.writerow(['Index', 'ORF Nucleotide Sequence', 'ORF Translated Sequence', 'Start Codon Position', 'Stop Codon Position', 'Pribnow-Box Distance from ORF', 'Pribnow-Box Position in Sequence', 'Nested ORF'])
        elif args.organism_type.lower() == "eukaryotic":
            writer.writerow(['Index', 'ORF Nucleotide Sequence', 'ORF Translated Sequence', 'Start Codon Position', 'Stop Codon Position', 'TATA-Box Distance from ORF', 'TATA-Box Position in Sequence', 'Nested ORF'])
        else:
            writer.writerow(['Index', 'ORF Nucleotide Sequence', 'ORF Translated Sequence', 'Start Codon Position', 'Stop Codon Position', 'Regulatory Element Distance to corresponding ORF', 'Regulatory Element Position', 'Nested ORF'])

        # Find open reading frames in the DNA sequence
        open_reading_frames = find_open_reading_frames(sequence, args.organism_type, regulatory_positions)

        # Write the details of identified ORFs to the CSV file
        for idx, orf in enumerate(open_reading_frames, start=1):
            writer.writerow([idx, orf[1], orf[2], orf[3], orf[4], orf[5], orf[6], orf[7]])

main()
