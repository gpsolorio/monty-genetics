/****************************************************************************\
 File: hvh.c
 Author: Ian Korf
 License: Public Domain
\****************************************************************************/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


int cmp_int(const void *a, const void *b) {
	return *(int*)b - *(int*)a;
}

/*
	The counting table is a sparse 3D matrix that maps to a sparser 4d matrix.

	The most common nucleotide is stored first. Then the next most common. Etc.
	There can be at most 4 different letters. For a sequencing depth of 3,
	there are the following possible values.

	3.0.0.0
	2.1.0.0
	1.1.1.0

	Although there are 4 possible letters, you only need to store 3: if the
	values don't add up to the total depth, you can infer the counts of the
	final letter. This means the matrix can be 3 dimensions instead of 4. Here
	are all the mappings for a depth of 4.

	4.0.0.0 stored as 4.0.0
	3.1.0.0 stored as 3.1.0
	2.1.1.0 stored as 2.1.0
	1.1.1.1 stored as 1.1.1 -> implied 1.1.1.1 because depth must be 4

	The matrix becomes increasingly sparse at higher depths. While this is
	somewhat wasteful, the memory isn't that large. However, it might be more
	cache-able if the size was smaller. But that is a problem for another day.

	Depth  Types  Matrix
	  2       2      27
	  3       3      64
	  4       5     125
	  5       6     216
	 10      23    1331
	 20     108    9261
	 30     297   29791
	 40     632   68921
*/

uint32_t ***count_table(int size) {
	uint32_t ***t = malloc(sizeof(uint32_t **) * size);
	for (int i = 0; i < size; i++) {
		t[i] = malloc(sizeof(uint32_t *) * size);
		for (int j = 0; j < size; j++) {
			t[i][j] = calloc(sizeof(uint32_t), size);
		}
	}
	return t;
}

static char *help = "\
usage: hvh [options] <depth>\n\
  -i <float> iterations [0] (0 means preset iterations)\n\
  -e <float> sequencing error rate [0.1]\n\
  -o <file>  optional output file\n\
  -s <int>   optional random seed\n\
  -z         show only non-zero values\n\
  -v         verbose\n\
";


int main(int argc, char **argv) {
	float iterations = 0;
	float err_rate = 0.1;
	FILE *output = stdout;
	int verbose = 0;
	int skip_zeroes = 0;
	int depth;

	// CLI
	srand(time(NULL));
	if (argc == 1) {printf("%s", help); exit(1);}
	int opt;
	while ((opt = getopt(argc, argv, "i:e:o:s:zv")) != -1) {
		switch (opt) {
			case 'i': iterations = atof(optarg); break;
			case 'e': err_rate = atof(optarg); break;
			case 'o': output = fopen(optarg, "w"); break;
			case 's': srand(atoi(optarg)); break;
			case 'v': verbose = 1; break;
			case 'z': skip_zeroes = 1; break;
		}
	}
	depth = atoi(argv[argc-1]);

	// preset iterations based on depth
	if (iterations == 0) {
		iterations = 1e8 * pow(10, (double)depth/10);
	}
	uint32_t limit = iterations;

	// Create the count lists for major and minor alleles
	uint32_t ***hets = count_table(depth +1);
	uint32_t ***homs = count_table(depth +1);

	// Homozygous
	if (verbose) fprintf(stderr, "Homozygous simulations (%g) ", iterations);
	uint32_t i = 0;
	uint32_t update = limit / 20;
	while (1) {
		if (verbose && i % update == 0) fprintf(stderr, ".");
		if (i++ >= limit) break;
		int count[4] = {0};
		for (int j = 0; j < depth; j++) {
			int slot = 0; // mom
			double r = rand()/(double)RAND_MAX;
			if (r < err_rate) slot = 1 + rand() % 3;
			count[slot] += 1;
		}
		qsort(count, 4, sizeof(int), cmp_int);
		homs[ count[0] ][ count[1] ][ count[2] ] += 1;
	}
	if (verbose) fprintf(stderr, "\n");

	// Heterozygous
	if (verbose) fprintf(stderr, "Heterozygous simulations (%g) ", iterations);
	i = 0;
	while (1) {
		if (verbose && i % update == 0) fprintf(stderr, ".");
		if (i++ >= iterations) break;
		int count[4] = {0};
		for (int j = 0; j < depth; j++) {
			int slot = rand() % 2; // 0:mom, 1:dad
			double r = rand()/(double)RAND_MAX;
			if (r < err_rate) {
				if (slot == 0) {
					slot = rand() % 3 + 1;
				} else {
					slot = rand() % 3 + 1;
					if (slot == 4) slot = 0;
				}
			}
			count[slot] += 1;
		}
		qsort(count, 4, sizeof(int), cmp_int);
		hets[ count[0] ][ count[1] ][ count[2] ] += 1;
	}
	if (verbose) fprintf(stderr, "\n");

	// Output
	fprintf(output, "#counts\thet\thom\tp(het)\t%g\n", iterations);
	for (int i = depth; i >= 0; i--) {
		for (int j = depth -i; j >= 0; j--) {
			if (j > i) continue;
			for (int k = depth -i -j; k >= 0; k--) {
				if (k > j) continue;
				int m = depth -i -j -k;
				if (m > k) continue;
				int hom = homs[i][j][k];
				int het = hets[i][j][k];
				if (skip_zeroes && (hom == 0 || het == 0)) continue;
				fprintf(output, "%d.%d.%d.%d\t%d\t%d\t%g\n", i, j, k, m,
					het, hom, (double)het / (double)(het + hom));
			}
		}
	}

	exit(0);
}
