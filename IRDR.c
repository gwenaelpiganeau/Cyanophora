/// This code extracts the longest Inverted Repeat (IR), Directe Repeat (DR) and Homopolymer Repeat (HR) from fasta sequences up to 10kb
/// Gwenael Piganeau 16/01/2025
///
///
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN_REPEAT_SIZE 4
#define MAX_SEQUENCE_LENGTH 10000
#define MAX_NAME_LENGTH 100
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

// Function to compute reverse complement of a sequence
void reverse_complement(const char *seq, char *rev_comp) {
    int len = strlen(seq);
    for (int i = 0; i < len; i++) {
        switch (seq[len - 1 - i]) {
            case 'A': rev_comp[i] = 'T'; break;
            case 'T': rev_comp[i] = 'A'; break;
            case 'C': rev_comp[i] = 'G'; break;
            case 'G': rev_comp[i] = 'C'; break;
            default: rev_comp[i] = 'N'; break;
        }
    }
    rev_comp[len] = '\0';
}

// Function to perform Needleman-Wunsch alignment and find the longest ungapped IR
void find_longest_IR(const char *sequence, int *lgIR, int *loop, char *forward_IR, char *reverse_IR) {
    int len = strlen(sequence);
    *lgIR = 0;
    *loop = 0;
    strcpy(forward_IR, "-");
    strcpy(reverse_IR, "-");

    char rev_comp[MAX_SEQUENCE_LENGTH];
    reverse_complement(sequence, rev_comp);

    int **score_matrix = (int **)malloc((len + 1) * sizeof(int *));
    for (int i = 0; i <= len; i++) {
        score_matrix[i] = (int *)malloc((len + 1) * sizeof(int));
        memset(score_matrix[i], 0, (len + 1) * sizeof(int));
    }

    int max_length = 0, max_i = 0, max_j = 0;

    // Fill the score matrix
    for (int i = 1; i <= len; i++) {
        for (int j = 1; j <= len; j++) {
            if (sequence[i - 1] == rev_comp[j - 1]) {
                score_matrix[i][j] = score_matrix[i - 1][j - 1] + MATCH_SCORE;
            } else {
                score_matrix[i][j] = 0; // No mismatches allowed for IR
            }

            if (score_matrix[i][j] > max_length) {
                max_length = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    if (max_length >= MIN_REPEAT_SIZE) {
        *lgIR = max_length;

        // Extract the forward IR and reverse IR
        strncpy(forward_IR, &sequence[max_i - max_length], max_length);
        forward_IR[max_length] = '\0';
        reverse_complement(forward_IR, reverse_IR);

        // Calculate the loop size: number of nucleotides between the forward IR and reverse IR
        int forward_end = max_i - 1;
        int reverse_start = len - max_j;
//add condition here to get lopp=0 if (max_i - max_length)== reverse_start && (forward_end== reverse_start + max_length - 1)
        if ( (max_i - max_length)== reverse_start && forward_end== (reverse_start + max_length - 1)) {
            *loop=0;
        }
            
       else {
           *loop = reverse_start - forward_end - 1;
        }

        // Debugging output for positions
        printf("Forward IR starts at: %d, ends at: %d\n", max_i - max_length, forward_end);
        printf("Reverse IR starts at: %d, ends at: %d\n", reverse_start, reverse_start + max_length - 1);
    }

    // Cleanup
    for (int i = 0; i <= len; i++) {
        free(score_matrix[i]);
    }
    free(score_matrix);

    if (*lgIR < MIN_REPEAT_SIZE) {
        *lgIR = 0;
        strcpy(forward_IR, "-");
        strcpy(reverse_IR, "-");
    }
}

// Function to find the longest direct repeat (DR)
void find_longest_DR(const char *sequence, int *lgDR, int *nbrepeat, char *DR) {
    int len = strlen(sequence);
    *lgDR = 0;
    *nbrepeat = 0;
    strcpy(DR, "-");

    // Start with the longest possible repeat length and work downward
    for (int sub_len = len / 2; sub_len >= MIN_REPEAT_SIZE; sub_len--) {
        for (int start = 0; start + sub_len <= len; start++) {
            int repeat_count = 1; // Initialize with 1 occurrence for the starting segment
            int pos = start + sub_len; // Start searching after the current segment

            // Check for non-overlapping occurrences
            while (pos + sub_len <= len) {
                if (strncmp(sequence + start, sequence + pos, sub_len) == 0) {
                    repeat_count++;
                    pos += sub_len; // Move past the current match (non-overlapping)
                } else {
                    pos++; // Slide the search window by one if no match
                }
            }

            // Update if this repeat is the longest valid DR
            if (repeat_count >= 2) { // At least 2 occurrences
                if (sub_len > *lgDR || (sub_len == *lgDR && repeat_count > *nbrepeat)) {
                    *lgDR = sub_len;
                    *nbrepeat = repeat_count;
                    strncpy(DR, sequence + start, sub_len);
                    DR[sub_len] = '\0';

                    // Stop further search since we're working from largest sub_len
                    return;
                }
            }
        }
    }

    // Reset results if no valid DR is found
    if (*lgDR < MIN_REPEAT_SIZE || *nbrepeat < 2) {
        *lgDR = 0;
        *nbrepeat = 0;
        strcpy(DR, "-");
    }
}

void find_longest_HR(const char *sequence, int *lgHR, char *HR) {
    int len = strlen(sequence);
    *lgHR = 0;
    strcpy(HR, "-");

    int max_length = 0;       // Longest homopolymer length found
    char max_char = '-';      // Character of the longest homopolymer
    int current_length = 1;   // Length of the current homopolymer

    // Traverse the sequence to detect homopolymer regions
    for (int i = 1; i < len; i++) {
        if (sequence[i] == sequence[i - 1]) {
            current_length++; // Extend the current homopolymer
        } else {
            // Check if the completed homopolymer meets the minimum length
            if (current_length >= MIN_REPEAT_SIZE && current_length > max_length) {
                max_length = current_length;
                max_char = sequence[i - 1];
            }
            current_length = 1; // Reset for the next homopolymer
        }
    }

    // Check the final homopolymer at the end of the sequence
    if (current_length >= MIN_REPEAT_SIZE && current_length > max_length) {
        max_length = current_length;
        max_char = sequence[len - 1];
    }

    // Store the results if a valid HR is found
    if (max_length >= MIN_REPEAT_SIZE) {
        *lgHR = max_length;
        memset(HR, max_char, max_length); // Fill HR with the homopolymer character
        HR[max_length] = '\0';           // Null-terminate the HR string
    } else {
        *lgHR = 0;
        strcpy(HR, "-");
    }
}


int main(int argc, char *argv[]) {
    int lgIR, lgHR, loop, lgDR, nbrepeat;
    char forward_IR[MAX_SEQUENCE_LENGTH], reverse_IR[MAX_SEQUENCE_LENGTH], DR[MAX_SEQUENCE_LENGTH], HR[MAX_SEQUENCE_LENGTH];
    char line[MAX_SEQUENCE_LENGTH];
    char name[MAX_NAME_LENGTH];
    char sequence[MAX_SEQUENCE_LENGTH];
    char outputFile[MAX_NAME_LENGTH];
    sequence[0] = '\0';

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_fasta_file>\n", argv[0]);
        return 1;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        fprintf(stderr, "Error opening input file: %s\n", argv[1]);
        return 1;
    }
    sprintf(outputFile, "%s_repeatoutput.txt", argv[1]);
    FILE *output = fopen(outputFile, "w");
    if (!output) {
        fprintf(stderr, "Error opening output file.\n");
        fclose(input);
        return 1;
    }


    fprintf(output, "Name_sequence\tlgIR\tl\tforward-IR\treverse-IR\tlgDR\tnbrepeat\tDR\tHR_bp\tHR\n");

    while (fgets(line, sizeof(line), input)) {
        if (line[0] == '>') {
            if (strlen(sequence) > 0) {

                find_longest_IR(sequence, &lgIR, &loop, forward_IR, reverse_IR);
                find_longest_DR(sequence, &lgDR, &nbrepeat, DR);
                find_longest_HR(sequence, &lgHR, HR);

                //fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, lgIR, loop, forward_IR, reverse_IR, lgDR, nbrepeat, DR, lgHR, HR);
                if (loop == 0) {
                    // For perfect palindromes, split the sequence into forward and reverse parts
                    int half_length = lgIR / 2;
                    char forward_half[half_length + 1], reverse_half[half_length + 1];
                    strncpy(forward_half, forward_IR, half_length);
                    forward_half[half_length] = '\0';
                    strncpy(reverse_half, forward_IR + half_length, half_length);
                    reverse_half[half_length] = '\0';

                    fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, half_length, loop,
                            forward_half, reverse_half, lgDR, nbrepeat, DR, lgHR, HR);
                } else {
                    // Normal case for inverted repeats
                    fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, lgIR, loop,
                            forward_IR, reverse_IR, lgDR, nbrepeat, DR, lgHR, HR);
                }

                sequence[0] = '\0';
            }
            sscanf(line, ">%s", name);
        } else {
            strcat(sequence, line);
            sequence[strcspn(sequence, "\n")] = '\0';
        }
    }

    if (strlen(sequence) > 0) {

        find_longest_IR(sequence, &lgIR, &loop, forward_IR, reverse_IR);
        find_longest_DR(sequence, &lgDR, &nbrepeat, DR);
        find_longest_HR(sequence, &lgHR, HR);
//        fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, lgIR, loop, forward_IR, reverse_IR, lgDR, nbrepeat, DR, lgHR, HR);
        if (loop == 0) {
            // For perfect palindromes, split the sequence into forward and reverse parts
            int half_length = lgIR / 2;
            char forward_half[half_length + 1], reverse_half[half_length + 1];
            strncpy(forward_half, forward_IR, half_length);
            forward_half[half_length] = '\0';
            strncpy(reverse_half, forward_IR + half_length, half_length);
            reverse_half[half_length] = '\0';

            fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, half_length, loop,
                    forward_half, reverse_half, lgDR, nbrepeat, DR, lgHR, HR);
        } else {
            // Normal case for inverted repeats
            fprintf(output, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\n", name, lgIR, loop,
                    forward_IR, reverse_IR, lgDR, nbrepeat, DR, lgHR, HR);
        }

    }

    fclose(input);
    fclose(output);
    printf("Analysis complete. Results written to %s\n",outputFile);
    return 0;
}
