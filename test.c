#include "closed_syncmers.h"

int main(int argc, char *argv[]) {
    if(argc < 4) {
        fprintf(stderr, "Usage: %s sequence K S [debug_mode]\n", argv[0]);
        return 1;
    }
    char *sequence_input = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);
    int debug_mode = 0;
    if(argc >= 5) {
        debug_mode = atoi(argv[4]);
    }
    if(S >= K) {
        fprintf(stderr, "Error: S must be less than K\n");
        return 1;
    }
    compute_closed_syncmers(sequence_input, strlen(sequence_input), K, S, debug_mode);
    return 0;
}

