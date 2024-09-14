// benchmark.c
#include "closed_syncmers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NUM_SEQUENCES 7000
#define SEQUENCE_LENGTH 10000
#define K 2000
#define S 40

int main() {
    char **sequences = malloc(NUM_SEQUENCES * sizeof(char *));
    if (!sequences) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Generate random sequences
    const char nucleotides[] = "ACGT";
    for (int i = 0; i < NUM_SEQUENCES; i++) {
        sequences[i] = malloc(SEQUENCE_LENGTH + 1);
        if (!sequences[i]) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }
        for (int j = 0; j < SEQUENCE_LENGTH; j++) {
            sequences[i][j] = nucleotides[rand() % 4];
        }
        sequences[i][SEQUENCE_LENGTH] = '\0'; // Null-terminate the string
    }

    // Time the computation
    clock_t start = clock();
    for (int i = 0; i < NUM_SEQUENCES; i++) {
        compute_closed_syncmers(sequences[i], SEQUENCE_LENGTH, K, S, 0);
    }
    clock_t end = clock();
    double total_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Compute throughput
    double total_bytes = (double)NUM_SEQUENCES * SEQUENCE_LENGTH;
    double throughput = (total_bytes / (1024 * 1024)) / total_time; // MB/s

    printf("Processed %d sequences of length %d in %f seconds\n", NUM_SEQUENCES, SEQUENCE_LENGTH, total_time);
    printf("Throughput: %f MB/s\n", throughput);

    // Clean up
    for (int i = 0; i < NUM_SEQUENCES; i++) {
        free(sequences[i]);
    }
    free(sequences);
    return 0;
}

