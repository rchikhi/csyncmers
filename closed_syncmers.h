// mostly o1-produced
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

// Main function to compute closed syncmers
void compute_closed_syncmers(const char *sequence_input, int len, int K, int S, int debug_mode);

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

// Deque node structure
typedef struct {
    uint64_t hash;
    size_t pos;
} DequeNode;

// Deque data structure for sliding window minimum
typedef struct {
    DequeNode *data;
    int front;
    int back;
    int capacity;
} Deque;

// Initialize deque
static inline void init_deque(Deque *dq, int capacity) {
    dq->data = (DequeNode *)malloc(capacity * sizeof(DequeNode));
    dq->front = 0;
    dq->back = 0;
    dq->capacity = capacity;
}

// Free deque memory
static inline void free_deque(Deque *dq) {
    free(dq->data);
}

// Get deque size
static inline int deque_size(Deque *dq) {
    return dq->back - dq->front;
}

// Push element to back of deque
static inline void deque_push_back(Deque *dq, DequeNode value) {
    dq->data[dq->back % dq->capacity] = value;
    dq->back++;
}

// Pop element from front of deque
static inline void deque_pop_front(Deque *dq) {
    dq->front++;
}

// Pop element from back of deque
static inline void deque_pop_back(Deque *dq) {
    dq->back--;
}

// Get front element of deque
static inline DequeNode deque_front(Deque *dq) {
    return dq->data[dq->front % dq->capacity];
}

// Get back element of deque
static inline DequeNode deque_back(Deque *dq) {
    return dq->data[(dq->back - 1) % dq->capacity];
}

// Compute closed syncmers
void compute_closed_syncmers(const char *sequence_input, int len, int K, int S, int debug_mode) {
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // Sliding window minimum over canonical s-mer hashes
    int window_size = K - S + 1; // Number of s-mers in a k-mer
    Deque dq;
    init_deque(&dq, window_size + 1); // Capacity must be greater than window_size
    uint64_t mask = ((uint64_t)1 << (2 * S)) - 1; // Mask to keep s-mer length
    uint64_t hash_fwd = 0, hash_rev = 0, canonical_hash;
    uint64_t rc_shift = 2 * (S - 1); // Shift amount for reverse complement
    for(size_t i = 0; i < len; i++) {
        uint8_t base = base_to_bits(sequence_input[i]);
        // Update forward hash
        if(i < S) {
            hash_fwd = (hash_fwd << 2) | base;
            hash_fwd &= mask;
        } else {
            hash_fwd = ((hash_fwd << 2) | base) & mask;
        }
        // Update reverse complement hash
        uint8_t comp_base = complement_base(base);
        if(i < S) {
            hash_rev = hash_rev | ((uint64_t)comp_base << (2 * i));
        } else {
            hash_rev = (hash_rev >> 2) | ((uint64_t)comp_base << rc_shift);
            hash_rev &= mask;
        }
        // Proceed only when we have a full s-mer
        if(i >= S - 1) {
            size_t s_mer_pos = i - S + 1;
            // Compute canonical hash
            canonical_hash = (hash_fwd < hash_rev) ? hash_fwd : hash_rev;
            if(debug_mode) {
                // Print the position and hash of the s-mer
                printf("s-mer at position %zu: hash=%" PRIu64 "\n", s_mer_pos, canonical_hash);
            }
            // Remove nodes outside the current window
            while(deque_size(&dq) > 0 && deque_front(&dq).pos <= s_mer_pos - window_size) {
                deque_pop_front(&dq);
            }
            // Remove nodes with greater or equal hash values
            while(deque_size(&dq) > 0 && deque_back(&dq).hash >= canonical_hash) {
                deque_pop_back(&dq);
            }
            // Add current s-mer hash and position to the deque
            DequeNode node = {canonical_hash, s_mer_pos};
            deque_push_back(&dq, node);
            // Now, if we have processed at least window_size s-mers
            if(s_mer_pos >= window_size - 1) {
                size_t window_start = s_mer_pos - (window_size - 1);
                size_t kmer_start = window_start;
                size_t min_pos_in_window = deque_front(&dq).pos;
                if(debug_mode) {
                    printf("Window starting at position %zu: min s-mer at position %zu\n",
                           window_start, min_pos_in_window);
                }
                if(min_pos_in_window == window_start || min_pos_in_window == s_mer_pos) {
                    // Output k-mer starting at position kmer_start
                    if(kmer_start + K <= len) {
                        //printf("%.*s\n", K, &sequence_input[kmer_start]);
                    }
                }
            }
        }
    }
    free_deque(&dq);
}

#ifdef __cplusplus
}
#endif

#endif // CLOSED_SYNCMERS_H

