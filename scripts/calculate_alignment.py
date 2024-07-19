import argparse
import csv
from Bio.Align.Applications import ClustalwCommandline, TCoffeeCommandline, MuscleCommandline
from Bio import AlignIO
from Bio.pairwise2 import format_alignment
from nltk.metrics import edit_distance

def main(input_file, alignment_method, shard_size, shard_position):
    # Parse input CSV file
    with open(input_file) as csvfile:
        reader = csv.DictReader(csvfile)

        # Create a list of scaffold sequences
        scaffolds = []

        # Loop over rows in the CSV
        for i, row in enumerate(reader):
            scaffolds.append(row['scaffold'])

            # If we have reached the shard size, process the current shard of sequences
            if i % shard_size == shard_position:
                # Create a FASTA file with the noisy repeat sequences
                with open('noisy_repeats.fasta', 'w') as fasta:
                    for repeat in noisy_repeats:
                        fasta.write(f">noisy_repeat\n{repeat}")

                # Run the specified multi-sequence alignment
                if alignment_method == "clustal":
                    clustal_cline = ClustalwCommandline("clustalw", infile="noisy_repeats.fasta")
                    clustal_cline()

                    # Parse the alignment
                    align = AlignIO.read("noisy_repeats.aln", "clustal")
                elif alignment_method == "tcoffee":
                    tcoffee_cline = TCoffeeCommandline(infile="noisy_repeats.fasta")
                    tcoffee_cline()

                    # Parse the alignment
                    align = AlignIO.read("noisy_repeats.aln", "stockholm")
                elif alignment_method == "muscle":
                    muscle_cline = MuscleCommandline(input="noisy_repeats.fasta")
                    muscle_cline()

                    # Parse the alignment
                    align = AlignIO.read("noisy_repeats.aln", "fasta")

                # Compare the edit distance of the alignment to the scaffold sequences
                for i, scaffold in enumerate(scaffolds):
                    # Get the i-th noisy repeat sequence from the alignment
                    repeat = align[i].seq

                    # Calculate the edit distance between the repeat and scaffold
                    edit_distance = calculate_edit_distance(repeat, scaffold)

                    # Print the edit distance
                    print(f"Edit distance for sequence {i+1}: {edit_distance}")

                # Clear the list of scaffold sequences and noisy repeat sequences
                scaffolds = []
                noisy_repeats = []
            else:
                # Add the noisy repeat sequence to the list
                noisy_repeats.append(row['noisy_repeats'])

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="Perform multi-sequence alignment on a CSV of noisy repeat sequences and compare the edit distance to the scaffold sequence.")
    parser.add_argument("-i", "--input", required=True, help="The input CSV file")
    parser.add_argument("-a", "--alignment", choices=["clustal", "tcoffee", "muscle"], default="clustal", help="The alignment method to use (default: clustal)")
    parser.add_argument("-s", "--shard-size", type=int, default=1000, help="The number of sequences to process at a time (default: 1000)")
    parser.add_argument("-p", "--shard-position", type=int, default=1, help="The position of the current shard of sequences (default: 1)")

    # Parse command line arguments
    args = parser.parse_args()

    # Run main function
    main(args.input, args.alignment, args.shard_size, args.shard_position)


