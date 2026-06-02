#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INIT_SEQ_CAPACITY 10000000
#define MAX_LINE_LEN 1000

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

void reverse_complement3(const char *in, char *out) {
    out[0] = complement(in[2]);
    out[1] = complement(in[1]);
    out[2] = complement(in[0]);
    out[3] = '\0';
}

int detect_motifs_and_write(const char *seq, int len, const char *seqname, FILE *summary, FILE *motifs) {
    int count = 0;

    for (int i = 0; i <= len - 11; i++) {
        if (seq[i + 3] != 'C' || seq[i + 7] != 'G') continue;

        char nnn1[4], nnn2[4], nnn3[4], revcomp1[4], revcomp2[4];
        strncpy(nnn1, &seq[i], 3); nnn1[3] = '\0';
        strncpy(nnn2, &seq[i + 4], 3); nnn2[3] = '\0';
        strncpy(nnn3, &seq[i + 8], 3); nnn3[3] = '\0';

        reverse_complement3(nnn1, revcomp1);
        if (strncmp(nnn2, revcomp1, 3) != 0) continue;

        reverse_complement3(nnn2, revcomp2);
        if (strncmp(nnn3, revcomp2, 3) != 0) continue;

        if (strncmp(nnn1, nnn3, 3) != 0) continue;  // ⚠️ condition ajoutée

        count++;
        fprintf(motifs, "%s\t%d\t", seqname, i);
        for (int j = 0; j < 11; j++) {
            fputc(seq[i + j], motifs);
        }
        fputc('\n', motifs);
    }

    if (count > 0) {
        fprintf(summary, "%s\t%d\n", seqname, count);
    }

    return count;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input.fasta> <output_summary.txt> <output_motifs.txt>\n", argv[0]);
        return 1;
    }

    FILE *fin = fopen(argv[1], "r");
    if (!fin) {
        perror("Erreur ouverture fichier FASTA");
        return 1;
    }

    FILE *fsummary = fopen(argv[2], "w");
    FILE *fmotifs = fopen(argv[3], "w");
    if (!fsummary || !fmotifs) {
        perror("Erreur ouverture fichiers sortie");
        fclose(fin);
        return 1;
    }

    char *sequence = malloc(INIT_SEQ_CAPACITY);
    if (!sequence) {
        perror("Échec malloc");
        return 1;
    }

    int capacity = INIT_SEQ_CAPACITY;
    int seqlen = 0;
    char line[MAX_LINE_LEN];
    char seqname[MAX_LINE_LEN] = "";
    int first = 1;

    while (fgets(line, sizeof(line), fin)) {
        if (line[0] == '>') {
            if (!first && seqlen >= 11) {
                sequence[seqlen] = '\0';
                detect_motifs_and_write(sequence, seqlen, seqname, fsummary, fmotifs);
            }
            sscanf(line, ">%s", seqname);
            seqlen = 0;
            first = 0;
        } else {
            size_t len = strcspn(line, "\r\n");
            if (seqlen + len + 1 > capacity) {
                capacity *= 2;
                char *new_seq = realloc(sequence, capacity);
                if (!new_seq) {
                    perror("Échec realloc");
                    free(sequence);
                    fclose(fin);
                    fclose(fsummary);
                    fclose(fmotifs);
                    return 1;
                }
                sequence = new_seq;
            }
            strncpy(&sequence[seqlen], line, len);
            seqlen += len;
        }
    }

    if (seqlen >= 11) {
        sequence[seqlen] = '\0';
        detect_motifs_and_write(sequence, seqlen, seqname, fsummary, fmotifs);
    }

    free(sequence);
    fclose(fin);
    fclose(fsummary);
    fclose(fmotifs);
    return 0;
}

