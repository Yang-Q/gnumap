/*
$Id: $

Description: Setup and find profile space centers

$Log: $
*/

#ifndef _CENTERS_H
#define _CENTERS_H

#include <stdio.h>

// DIM_CENTERS is the number of dimensions (nucleotides, etc)
#define DIM_CENTERS 5
#define N_CENTERS ccount
#if DIM_CENTERS < 2
#error DIM_CENTERS must be greater than or equal to 2
#endif
#define BAND (DIM_CENTERS-1)

// LEVELS is the number of _____________
#define LEVELS 5
#if LEVELS < BAND && METHOD_J
#error LEVELS must be greater than or equal to DIM_CENTERS-1
#endif

// RLEVELS is the granularity at each dimension.  Must be a power of 2
#define RLEVELS 32
#if RLEVELS < BAND
#error RLEVELS must be greater than or equal to DIM_CENTERS-1
#endif

typedef float center_t;
typedef unsigned char center_d;	/* center_d(atatype) */

extern size_t ccount;
extern center_t **centers;
extern unsigned int** sum_mat;

void init_centers(void);
void print_centers(FILE *stream);
center_t distance(center_t *pos1, center_t *pos2);
int find_center_exhaustive(center_t *pos);
void init_lookup(void);
void print_stats(FILE *stream);
int find_center(center_t *pos);
void init_sums(void);
center_d adjust_center(float prevV, float newV, center_d d);

#endif /* _CENTERS_H */
