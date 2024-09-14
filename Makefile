all:
	gcc -o test test.c
	gcc -O3 -march=native -o benchmark benchmark.c
