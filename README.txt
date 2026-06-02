This script analyzes a FASTA file and searches for a very specific 11-nucleotide motif with a defined palindromic structure:
NNN1 = positions 0-2
position 3 = C
NNN2 = positions 4-6
position 7 = G
NNN3 = positions 8-10

Where:
NNN2 = reverse-complement of NNN1
NNN3 = reverse-complement of NNN2
NNN1 = NNN3

Example of searched 11-nucleotide motif:
ATG  C CAT  G ATG
NNN1 C NNN2 G NNN3

Usage:
./NNNCNNNGNNN input.fasta output_summary.txt output_motifs.txt
With
input.fasta = your input FASTA file;
output_summary.txt = summary file (number of motifs per sequence);
output_motifs.txt = file containing each detected motif with its position.