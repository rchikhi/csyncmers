#ifndef CLOSED_SYNCMERS_NAIVE_H
#define CLOSED_SYNCMERS_NAIVE_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include "closed_syncmers.h"

// Naive computation of closed syncmers
void compute_closed_syncmers_naive(const char *sequence, size_t seq_len, int K, int S, MinimizerResult *results, int *num_results) {
    *num_results = 0;
    __uint128_t mask = (((__uint128_t)1) << (2*S)) - 1;

    // For each k-mer in the sequence
    for (size_t i = 0; i <= seq_len - K; i++) {
        __uint128_t min_hash = ~((__uint128_t)0);
        size_t min_pos_in_kmer = 0;

        // For each s-mer within the k-mer
        for (size_t j = 0; j <= K - S; j++) {
            size_t s_mer_pos = i + j;
            __uint128_t hash_fwd = 0;
            __uint128_t hash_rev = 0;

            // Compute the hash of the s-mer at position s_mer_pos
            for (size_t k = 0; k < S; k++) {
                uint8_t base = base_to_bits(sequence[s_mer_pos + k]);
                // Update forward hash
                hash_fwd = ((hash_fwd << 2) | base) & mask;
                // Update reverse complement hash
                uint8_t comp_base = complement_base(base);
                hash_rev |= ((__uint128_t)comp_base << (2*k));
            }

            // Compute canonical hash
            __uint128_t canonical_hash = (hash_fwd < hash_rev) ? hash_fwd : hash_rev;
            if (canonical_hash < min_hash) {
                min_hash = canonical_hash;
                min_pos_in_kmer = j;
            }
        }

        // Check if the minimal s-mer is at the first or last position within the k-mer
        if (min_pos_in_kmer == 0 || min_pos_in_kmer == K - S) {
            // Record the position and minimizer hash
            results[*num_results].position = i;
            results[*num_results].minimizer_hash = min_hash;
            (*num_results)++;
        }
    }
}


#endif // CLOSED_SYNCMERS_NAIVE_H

