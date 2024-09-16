#ifndef CLOSED_SYNCMERS_H
#define CLOSED_SYNCMERS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/*** Function Declarations ***/

// Struct to hold minimizer and its position
typedef struct {
    __uint128_t minimizer_hash;
    size_t position;
} MinimizerResult;

// Main function to compute closed syncmers
void compute_closed_syncmers(const char *sequence_input, int len, int K, int S, MinimizerResult *result, int *num_results);

/*** Implementation ***/

// Convert nucleotide base to 2-bit representation
static inline uint8_t base_to_bits(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0; // Treat unknown as 'A'
    }
}

// Get complement base in 2-bit encoding
static inline uint8_t complement_base(uint8_t base) {
    return 3 - base; // Complement: A<->T, C<->G
}

// Dynamically resize array for minimizers
static void add_minimizer(MinimizerResult *results, int *size, __uint128_t minimizer_hash, size_t position) {
    results[*size].minimizer_hash = minimizer_hash;
    results[*size].position = position;
    (*size)++;
}

// Compute closed syncmers
void compute_closed_syncmers(const char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results) {
    *num_results = 0;
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    size_t num_s_mers = len - S + 1;
    __uint128_t *s_mer_hashes = (__uint128_t *)malloc(num_s_mers * sizeof(__uint128_t));

    // Precompute all s-mer hashes
    __uint128_t mask = (((__uint128_t)1) << (2 * S)) - 1;
    __uint128_t hash_fwd = 0, hash_rev = 0;
    __uint128_t rc_shift = 2 * (S - 1);

    for(size_t i = 0; i < len; i++) {
        uint8_t base = base_to_bits(sequence_input[i]);
        hash_fwd = ((hash_fwd << 2) | base) & mask;
        uint8_t comp_base = complement_base(base);
        hash_rev = ((hash_rev >> 2) | ((__uint128_t)comp_base << rc_shift)) & mask;
        if(i >= S - 1) {
            size_t s_mer_pos = i - S + 1;
            __uint128_t canonical_hash = (hash_fwd < hash_rev) ? hash_fwd : hash_rev;
            s_mer_hashes[s_mer_pos] = canonical_hash;
        }
    }

    // Iterate over k-mers
    size_t num_k_mers = len - K + 1;
    for(size_t i = 0; i < num_k_mers; i++) {
        size_t start_s_mer_pos = i;
        size_t end_s_mer_pos = i + K - S;
        __uint128_t min_hash = s_mer_hashes[start_s_mer_pos];
        size_t min_pos = start_s_mer_pos;

        // Find minimal s-mer in the k-mer
        for(size_t j = start_s_mer_pos; j <= end_s_mer_pos; j++) {
            if(s_mer_hashes[j] < min_hash) {
                min_hash = s_mer_hashes[j];
                min_pos = j;
            }
        }

        // Check if minimal s-mer is at start or end
        if(min_pos == start_s_mer_pos || min_pos == end_s_mer_pos) {
            //printf("%.*s\n", K, &sequence_input[i]);
            //printf("%ld (%ld) ", kmer_start, (uint64_t)canonical_hash);
            add_minimizer(results, num_results, min_hash, i);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // CLOSED_SYNCMERS_H

