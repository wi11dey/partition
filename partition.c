#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define INTEGERS 100
#define MAX_ITERATIONS 25000

#define KARMARKAR_KARP 0
#define REPEATED_RANDOM 1
#define HILL_CLIMBING 2
#define SIMULATED_ANNEALING 3
#define PREPARTITIONED_REPEATED_RANDOM 11
#define PREPARTITIONED_HILL_CLIMBING 12
#define PREPARTITIONED_SIMULATED_ANNEALING 13

static inline double T(int iter) {
	return 10000000000*pow(0.8, iter/300);
}

static inline bool randbool() {
	return rand() < (RAND_MAX/2);
}

static inline int randindex() {
	return (int) ((((double) rand())/((double) RAND_MAX))*INTEGERS);
}

struct heap {
	long* array;
	int capacity;
	int size;
};

static inline struct heap* heap_new(int capacity) {
	long* array = malloc(sizeof(long)*capacity);
	struct heap* heap = malloc(sizeof(struct heap));
	if (!array || !heap) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}

	heap->array = array;
	heap->capacity = capacity;
	heap->size = 0;

	return heap;
}

static inline void heap_free(struct heap* heap) {
	free(heap->array);
	free(heap);
}

static inline int heap_parent(int current) {
	return (current - 1)/2;
}

static inline int heap_left(int current) {
	return (2*current) + 1;
}

static inline int heap_right(int current) {
	return heap_left(current) + 1;
}

static inline void heap_swap(struct heap* heap, int a, int b) {
	long temp = heap->array[a];
	heap->array[a] = heap->array[b];
	heap->array[b] = temp;
}

void heap_insert(struct heap* heap, long value) {
	if (heap->size == heap->capacity) {
		fprintf(stderr, "heap at capacity\n");
		return;
	}
	heap->array[heap->size] = value;
	int current = heap->size;
	heap->size++;
	int parent;
	while (current && heap->array[parent = heap_parent(current)] < heap->array[current]) {
		heap_swap(heap, parent, current);
		current = parent;
	}
}

void heapify(struct heap* heap, int start) {
	int left = heap_left(start),
		right = heap_right(start),
		max = start;
	if (left < heap->size && heap->array[left] > heap->array[max]) {
		max = left;
	}
	if (right < heap->size && heap->array[right] > heap->array[max]) {
		max = right;
	}
	if (max != start) {
		heap_swap(heap, start, max);
		heapify(heap, max);
	}
}

long heap_deletemax(struct heap* heap) {
	if (!heap->size) {
		fprintf(stderr, "empty heap\n");
		exit(1);
	}
	long max = heap->array[0];
	heap->size--;
	heap->array[0] = heap->array[heap->size];
	heapify(heap, 0);
	return max;
}

long karmarkar_karp(long* A) {
	struct heap* heap = heap_new(INTEGERS);
	for (int i = 0; i < INTEGERS; i++) {
		heap_insert(heap, A[i]);
	}
	while (heap->size >= 2) {
		long a = heap_deletemax(heap),
			b = heap_deletemax(heap);
		heap_insert(heap, labs(a - b));
	}
	long result = heap_deletemax(heap);
	heap_free(heap);
	return result;
}

// N.B. Always allocates, caller responsible for freeing.
long* standard_generate() {
	long* S = malloc(sizeof(long)*INTEGERS);
	for (int i = 0; i < INTEGERS; i++) {
		S[i] = randbool() ? 1 : -1;
	}
	return S;
}

static inline long standard_residue(long* A, long* S) {
	long residue = 0;
	for (int i = 0; i < INTEGERS; i++) {
		residue += S[i]*A[i];
	}
	return labs(residue);
}

// N.B. Always allocates, caller responsible for freeing.
long* standard_neighbor(long* S) {
	long* Sprime = malloc(sizeof(long)*INTEGERS);
	memcpy(Sprime, S, sizeof(long)*INTEGERS);
	int i = randindex();
	Sprime[i] = -Sprime[i];

	// Swap with 50% probability:
	if (randbool()) {
		int j;
		do {
			j = randindex();
		} while (j == i);
		Sprime[j] = -Sprime[j];
	}

	return Sprime;
}

long repeated_random(long* A) {
	long* S = standard_generate();
	for (int i = 0; i < MAX_ITERATIONS; i++) {
		long* Sprime = standard_generate();
		if (standard_residue(A, Sprime) < standard_residue(A, S)) {
			free(S);
			S = Sprime;
		} else {
			free(Sprime);
		}
	}
	long residue = standard_residue(A, S);
	free(S);
	return residue;
}

long hill_climbing(long* A) {
	long* S = standard_generate();
	for (int i = 0; i < MAX_ITERATIONS; i++) {
		long* Sprime = standard_neighbor(S);
		if (standard_residue(A, Sprime) < standard_residue(A, S)) {
			free(S);
			S = Sprime;
		} else {
			free(Sprime);
		}
	}
	long residue = standard_residue(A, S);
	free(S);
	return residue;
}

long simulated_annealing(long* A) {
	long* S = standard_generate();
	long* Sprimeprime = S;
	for (int i = 0; i < MAX_ITERATIONS; i++) {
		long* Sprime = standard_neighbor(S);

		long S_residue = standard_residue(A, S);
		long Sprime_residue = standard_residue(A, Sprime);

		if (Sprime_residue < S_residue || rand() < exp(-((double) (Sprime_residue - S_residue))/T(i))*RAND_MAX) {
			if (S != Sprimeprime) {
				free(S);
			}
			S = Sprime;
		} else {
			free(Sprime);
		}

		if (standard_residue(A, S) < standard_residue(A, Sprimeprime)) {
			free(Sprimeprime);
			Sprimeprime = S;
		}
	}
	long residue = standard_residue(A, Sprimeprime);
	if (Sprimeprime != S) {
		free(Sprimeprime);
	}
	free(S);
	return residue;
}

static long (* const ALGORITHMS[])(long* A) = {
	[KARMARKAR_KARP] = karmarkar_karp,
	[REPEATED_RANDOM] = repeated_random,
	[HILL_CLIMBING] = hill_climbing,
	[SIMULATED_ANNEALING] = simulated_annealing
};

int main(int argc, char** argv) {
	// Make sure every line is flushed immediately:
	setvbuf(stdout, NULL, _IOLBF, 0);

	srand(time(NULL));

	if (argc < 4) {
		fprintf(stderr, "usage: %s <flag> <algorithm> <inputfile>\n", (argc > 0) ? argv[0] : "./partition");
		return 1;
	}

	int algorithm = atoi(argv[2]);
	char* inputfile = argv[3];

	long A[INTEGERS];

	FILE* input = fopen(inputfile, "r");
	if (!input) {
		fprintf(stderr, "error opening %s\n", inputfile);
		return 1;
	}
	for (int i = 0; i < INTEGERS; i++) {
		if (fscanf(input, "%ld", &A[i]) != 1) {
			fprintf(stderr, "could only read %d integers from %s\n", i, inputfile);
			return 1;
		}
	}
	fclose(input);

	printf("%ld\n", ALGORITHMS[algorithm](A));

	return 0;
}
