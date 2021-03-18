// license:BSD-3-Clause
// copyright-holders:Jonathan Gevaryahu
// (C) 2020-2021 Jonathan Gevaryahu AKA Lord Nightmare
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

// basic typedefs
typedef int8_t s8;
typedef uint8_t u8;
typedef int16_t s16;
typedef uint16_t u16;
typedef int32_t s32;
typedef uint32_t u32;
typedef int64_t s64;
typedef uint64_t u64;

// defines
// max allowed threshold for absolute error for a match for two floats due to precision loss
#define THRESH_MAXALLOWED 0.001
// same for doubles
#define THRESH_DMAXALLOWED 0.0000000001

// verbosity defines; V_ERR is effectively 'always'
#define V_ERR   (1)
#define V_PARAM (c.verbose & (1<<0))
#define V_GEN   (c.verbose & (1<<1))
#define V_GEN2  (c.verbose & (1<<2))
#define V_FIX   (c.verbose & (1<<3))
#define V_FIX2  (c.verbose & (1<<4))
#define V_LOAD  (c.verbose & (1<<5))
#define V_LOAD2 (c.verbose & (1<<6))
#define V_FREE  (c.verbose & (1<<7))
#define V_MATH  (c.verbose & (1<<8))

// verbose macro
#define eprintf(v, ...) \
	do { if (v) { fprintf(stderr, __VA_ARGS__); fflush(stderr); } } while (0)

// struct definitions
typedef struct f_array
{
	s32 count; // this is the number of this specific letter used in this table, the numerator for each letter
	float cdf_data;
	float pdf_data;
} f_array;

typedef struct cdf_array
{
	f_array* start;
	f_array* middle;
	f_array* end;
	// these are the number of letters with non-zero counts within each of the tables
	s32 start_buckets;
	s32 middle_buckets;
	s32 end_buckets;
	// these are the total counts for number of letters used to generate the 3 cdf tables, the common denominator for all letters within each array
	s32 start_total;
	s32 middle_total;
	s32 end_total;
} cdf_array;

typedef struct ltrfile
{
	char magic[8];
	u8 num_letters;
	cdf_array* singles;
	cdf_array** doubles;
	cdf_array*** triples;
} ltrfile;

typedef struct s_cfg
{
	const char* const letters;
	int printcdf;
	u32 generate;
	u32 seed;
	bool fix;
	u32 verbose;
	bool dumpstart;
} s_cfg;


/* gcd */
u32 gcd(u32 a, u32 b)
{
	while (b)
	{
		a %= b;
		// clever way to do a swap(a,b) with no other variable/register:
		a ^= b;
		b ^= a;
		a ^= b;
	}
	return a;
}

/* lcm */
u32 lcm(u32 a, u32 b)
{
	return abs(a*b)/gcd(a,b);
}

/* msrand implementation for consistency */
#define MSRAND_MAX 0x7fff
static u32 state = 1;

u32 ms_rand()
{
	state = ((state * 214013) + 2531011) % (1 << 31);
	return (state >> 16) & 0x7fff;
}

void ms_srand(u32 seed)
{
	state = seed;
}

static float nrand() { return (float)ms_rand() / MSRAND_MAX; }
/* end msrand */

// for stock rand():
//static float nrand() { return (float)rand() / RAND_MAX; }

float fget_f(FILE *in)
{
	u32 t = 1;
	if ((*(char *) &t) == 1) // little endian
	{
		for (u32 i = 0; i < 4; i++)
		{
			t >>= 8;
			t |= ((u32)fgetc(in))<<24;
		}
	}
	else // big endian, this isn't tested yet!
	{
		t = 0;
		for (u32 i = 0; i < 4; i++)
		{
			t <<= 8;
			t |= fgetc(in);
		}
	}
	return *(float *) &t;
}

f_array* f_alloc(u32 count)
{
	f_array* f = malloc(count * sizeof(f_array));
	if (f == NULL)
	{
		eprintf(V_ERR,"E* Failure to allocate memory for array of size %d, aborting!\n", count);
		return NULL;
	}
	return f;
}

cdf_array* cdf_alloc(u32 count)
{
	cdf_array *c = malloc(sizeof(cdf_array));
	if (c == NULL)
	{
		eprintf(V_ERR,"E* Failure to allocate memory for array, aborting!\n");
		exit(1);
	}
	c->start = f_alloc(count);
	c->middle = f_alloc(count);
	c->end = f_alloc(count);
	c->start_buckets = c->middle_buckets = c->end_buckets = 0;
	c->start_total = c->middle_total = c->end_total = 0;
	return c;
}

void cdf_free(cdf_array *c)
{
	free(c->end);
	free(c->middle);
	free(c->start);
	free(c);
}

ltrfile* ltr_load(FILE *in, u32 len, s_cfg c)
{
	ltrfile *l = malloc(sizeof(ltrfile));
	u32 pos = 0;
	// magic
	{ // scope-limit
		const char* const compare = "LTR V1.0 ";
		for (u32 i = 0; i < 8; i++)
		{
			l->magic[i] = fgetc(in);
			if (l->magic[i] != compare[i])
			{
				eprintf(V_ERR,"E* Incorrect magic number! Exiting!\n");
				fclose(in);
				exit(1);
			}
			pos++;
		}
	}
	// number of letters
	l->num_letters = fgetc(in);
	if ((l->num_letters < 1) || (l->num_letters > 28))
	{
		eprintf(V_ERR,"E* Invalid number of letters %d! Exiting!\n", l->num_letters);
		fclose(in);
		exit(1);
	}
	pos++;
	eprintf(V_LOAD,"D* LTR header read ok, num_letters = %d\n", l->num_letters);

#define LOAD_LTR_FLOATS(x) \
	{ \
		float acc = 0.0; \
		for (u32 i = 0; i < l->num_letters; i++) \
		{ \
			x[i].cdf_data = fget_f(in); \
			x[i].pdf_data = (x[i].cdf_data) ? (x[i].cdf_data - acc) : 0.0; \
			x[i].count = -1; \
			pos += 4; \
			if (x[i].cdf_data) acc = x[i].cdf_data; \
		} \
	}

	// There was a bug in the original code Bioware used to create .ltr files
	// which caused the single.middle and single.end tables to have their CDF
	// values corrupted for all entries past any which have a probability of
	// zero. Fortunately, this can be corrected for in post, which we do here.

	// If the final nonzero value in the table is not 'exactly' 1.0, then they
	// are corrupt.
	// Note that likely due to precision loss sometime during generation by
	// Bioware's utility, the results, even after correction, may not exactly
	// accumulate to 1.000000f, so we give a small bit of leeway.
#define FIX_LTR_FLOATS(x) \
	{ \
		bool iscorrupt = true; \
		for (u32 i = 0; i < l->num_letters; i++) \
		{ \
			if (fabs(x[i].cdf_data - 1.0) <= THRESH_MAXALLOWED) \
				iscorrupt = false; \
		} \
		if (iscorrupt && c.fix) \
		{ \
			float acc = 0.0; \
			float prev = 0.0; \
			float correction = 0.0; \
			float uncorrected = 0.0; \
			eprintf(V_FIX,"Correcting errors in a probability table...\n"); \
			for (u32 i = 0; i < l->num_letters; i++) \
			{ \
				uncorrected = x[i].cdf_data; \
				if (x[i].cdf_data) \
				{ \
					if ((i > 0) && (prev == 0.0)) \
						correction = acc; \
					acc = x[i].cdf_data + correction; \
					x[i].cdf_data = acc; \
				} \
				x[i].pdf_data = (x[i].cdf_data) ? (x[i].cdf_data - correction) : 0.0; \
				x[i].count = -1; \
				pos += 4; \
				eprintf(V_FIX2,"ltr: %c, original: %f, corrected: %f, acc: %f, offset: %f\n", c.letters[i], uncorrected, x[i].cdf_data, acc, correction); \
				prev = uncorrected; \
			} \
			if (fabs(acc - 1.0) > THRESH_MAXALLOWED) \
				eprintf(V_FIX2,"*W during fixing process, accumulator ended up at a potentially incorrect value of %f!\n", acc); \
		} \
	}

	// allocate and fill singles table
	{ // scope-limit
		l->singles = cdf_alloc(l->num_letters);
		if (l->singles == NULL)
		{
			fclose(in);
			exit(1);
		}
		eprintf(V_LOAD2,"D* successfully allocated the singles cdf table\n");
		LOAD_LTR_FLOATS(l->singles->start);
		LOAD_LTR_FLOATS(l->singles->middle);
		FIX_LTR_FLOATS(l->singles->middle);
		LOAD_LTR_FLOATS(l->singles->end);
		FIX_LTR_FLOATS(l->singles->end);
		eprintf(V_LOAD2,"D* successfully filled the singles cdf table\n");
	}
	// allocate and fill doubles tables
	{ // scope-limit
		l->doubles = malloc( sizeof(l->doubles) * l->num_letters ); // allocate an array of pointers
		for (u32 j = 0; j < l->num_letters; j++)
		{
			l->doubles[j] = cdf_alloc(l->num_letters);
			if (l->doubles[j] == NULL)
			{
				fclose(in);
				exit(1);
			}
			eprintf(V_LOAD2,"D* successfully allocated the doubles cdf table %d\n", j);
			LOAD_LTR_FLOATS(l->doubles[j]->start);
			LOAD_LTR_FLOATS(l->doubles[j]->middle);
			LOAD_LTR_FLOATS(l->doubles[j]->end);
			eprintf(V_LOAD2,"D* successfully filled the doubles cdf table %d\n", j);
		}
	}
	// allocate and fill triples tables
	{ // scope-limit
		l->triples = malloc( sizeof(l->triples) * l->num_letters ); // allocate an array of pointers
		for (u32 k = 0; k < l->num_letters; k++)
		{
			l->triples[k] = malloc( sizeof(l->triples) * l->num_letters ); // allocate an array of pointers
			for (u32 j = 0; j < l->num_letters; j++)
			{
				l->triples[k][j] = cdf_alloc(l->num_letters);
				if (l->triples[k][j] == NULL)
				{
					fclose(in);
					exit(1);
				}
				eprintf(V_LOAD2,"D* successfully allocated the triples cdf table %d:%d\n", k, j);
				LOAD_LTR_FLOATS(l->triples[k][j]->start);
				LOAD_LTR_FLOATS(l->triples[k][j]->middle);
				LOAD_LTR_FLOATS(l->triples[k][j]->end);
				eprintf(V_LOAD2,"D* successfully filled the triples cdf table %d:%d\n", k, j);
			}
		}
	}
	eprintf(V_LOAD,"D* all tables loaded, expected size was %d, final size was %d\n", len, pos);
	return l;
}

void ltr_free(ltrfile* l, s_cfg c)
{
	// first clear triples
	for (u32 k = 0; k < l->num_letters; k++)
	{
		for (u32 j = 0; j < l->num_letters; j++)
		{
			cdf_free(l->triples[k][j]);
		}
		free(l->triples[k]);
	}
	free(l->triples);
	// then doubles
	for (u32 j = 0; j < l->num_letters; j++)
	{
		cdf_free(l->doubles[j]);
	}
	free(l->doubles);
	// then singles
	cdf_free(l->singles);
	// then free the ltrfile itself
	free(l);
	eprintf(V_FREE,"D* everything is freed!\n");
}

double get_mean_squared_error(f_array* input, u32 factor, u8 num_letters)
{
	double acc = 0.0;
	for (u8 i = 0; i < num_letters; i++)
	{
		u32 nearest_int = round(input[i].pdf_data * factor);
		acc += pow((input[i].pdf_data * factor) - (double)nearest_int, 2.0);
	}
	return acc/num_letters;
}

u32 f_array_count_buckets(f_array* f, u8 num_letters, s_cfg c)
{
	u32 count = 0;
	for (u32 i = 0; i < num_letters; i++)
	{
		if (f[i].pdf_data != 0.0)
			count++;
	}
	return count;
}

u32 f_array_analyze(f_array* f, u8 num_letters, s_cfg c)
{
	float minimum = 1.1;
	for (u32 i = 0; i < num_letters; i++)
	{
		if ((f[i].pdf_data != 0.0) && (f[i].pdf_data < minimum))
			minimum = f[i].pdf_data;
	}
	if ((minimum == 1.1) || (fabs(1.1 - minimum) < THRESH_MAXALLOWED)) // we know for sure this table is empty.
		return 0;
	if (minimum == 1.0) // we know for sure this table has a single 1.0 value in it.
		return 1;
	if (minimum == 0.5) // we know for sure this table has exactly two 0.5 values in it.
		return 2;
	// beyond that we can't prove anything.

	// based on this, we can take an educated first guess
	double invmin = 1.0 / minimum;
	eprintf(V_ERR,"D* inverse of minimum probability %f is %f\n", minimum, invmin);
	// but if this is not an integer, try multiplying it by a bunch of guesses to see if we can get it close to an integer.

	if (fabs(invmin - round(invmin)) > THRESH_DMAXALLOWED)
	{
		double min_error = 1000000.0;
		u32 best_guess = 0;
		for (u32 guess = 2; guess < 100; guess++)
		{
			double this_error = fabs((invmin*guess) - round(invmin*guess));
			//eprintf(V_ERR,"D* factor: testing guess of %d (error of %f), min error is %f\n", guess, this_error, min_error);
			if (this_error < min_error)
			{
				//eprintf(V_ERR,"D* factor: got a better guess (with an error of %f) of %d\n", this_error, guess);
				best_guess = guess;
				min_error = this_error;
				if (min_error < THRESH_DMAXALLOWED)
					break;
			}
		}
		eprintf(V_ERR,"D* factor: guessed error factor is %d, yielding %f\n", best_guess, invmin*best_guess);
	}

	// find the smallest possible integer that can be multiplied against all of
	// the PDF numbers that results in the results being completely integer with
	// nothing after the decimal point.
	// this is the number of words used to generate this list in the first place.
	// note we could be using the inverse of the minimum as an initial guess, but
	// this occasionally causes problems.
	double min_error = 1000000.0;
	u32 best_guess = 0;
	for (u32 guess = 3; guess < 30000; guess++) // guess could in theory be initially set to the best_guess variable from the function above this one, instead of 3
	{
		double this_error = get_mean_squared_error(f, guess, num_letters);
		if (this_error < min_error) // we have a better guess!
		{
			eprintf(V_ERR,"D* got a better guess (with an error of %f) of %d\n", this_error, guess);
			best_guess = guess;
			min_error = this_error;
			if (min_error < THRESH_DMAXALLOWED)
				break;
		}
	}
	//eprintf(V_ERR,"D* best guess (with an error of %f) was %d\n", min_error, best_guess);
	return best_guess;
}

// fill in the f_array->count values
void f_array_populate(f_array* f, u8 num_letters, u32 count, s_cfg c)
{
	for (u32 i = 0; i < num_letters; i++)
		f[i].count = round(f[i].pdf_data * count);
}

void cdf_analyze(cdf_array* p, u8 num_letters, s_cfg c)
{
	p->start_buckets = f_array_count_buckets(p->start, num_letters, c);
	p->start_total = f_array_analyze(p->start, num_letters, c);
	f_array_populate(p->start, num_letters, p->start_total, c);
	p->middle_buckets = f_array_count_buckets(p->middle, num_letters, c);
	p->middle_total = f_array_analyze(p->middle, num_letters, c);
	f_array_populate(p->middle, num_letters, p->middle_total, c);
	p->end_buckets = f_array_count_buckets(p->end, num_letters, c);
	p->end_total = f_array_analyze(p->end, num_letters, c);
	f_array_populate(p->end, num_letters, p->end_total, c);
}

// figure out whether the two integers are exact multiples
// this is done by dividing the larger by the smaller and seeing
// if the resulting number has a fractional portion that is exactly zero
bool is_exact_multiple(u32 i, u32 j)
{
	if (i == j)
		return true;
	float fi = i;
	float fj = j;
	float fg = fi;
	float fl = fj;
	if (fj > fg)
	{
		fg = fj;
		fl = fi;
	}
	if (fabs((fg / fl) - round(fg / fl)) > THRESH_MAXALLOWED)
		return false;
	else
		return true;
}

void ltr_analyze(ltrfile* l, s_cfg c)
{
	cdf_analyze(l->singles, l->num_letters, c);
	for (u32 j = 0; j < l->num_letters; j++)
	{
		cdf_analyze(l->doubles[j], l->num_letters, c);
	}
	for (u32 k = 0; k < l->num_letters; k++)
	{
		for (u32 j = 0; j < l->num_letters; j++)
		{
			cdf_analyze(l->triples[k][j], l->num_letters, c);
		}
	}

	// the numbers of names should be correct but could be off by some factor, so lets do some heuristics to correct this
	// first: whichever of the counts for the singles->start_total and singles->end_total is higher is automatically correct, if they're not the same AND are an even multiple of one another
	// some files were hand-edited and hence the two singles tables will not be an exact multiple of one another, and if so, don't try to fix them here!
	if (is_exact_multiple(l->singles->start_total,l->singles->end_total))
	{
		if (l->singles->start_total > l->singles->end_total) // start was higher
		{
			u32 c_factor = l->singles->start_total / l->singles->end_total;
			eprintf(V_ERR,"D* fixing singles->end table by factor of %d:\n", c_factor);
			// iterate through the table and correct the numerators
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->singles->end[i].count *= c_factor;
			}
			// correct the denominator
			l->singles->end_total = l->singles->start_total;
		}
		else if (l->singles->start_total < l->singles->end_total) // end was higher
		{
			u32 c_factor = l->singles->end_total / l->singles->start_total;
			eprintf(V_ERR,"D* fixing singles->start table by factor of %d:\n", c_factor);
			// iterate through the table and correct the numerators
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->singles->start[i].count *= c_factor;
			}
			// correct the denominator
			l->singles->start_total = l->singles->end_total;
		}
	}
	else
	{
		eprintf(V_ERR,"D* cannot equalize start and end tables as they are not an even factor of one another!\n");
	}

	// next heuristic: if the singles->start[*] count for letter * doesn't equal the denominator for doubles[*]->start_total but is off by some factor, increase the latter to match
	for (u32 i = 0; i < l->num_letters; i++)
	{
		if (l->singles->start[i].count != l->doubles[i]->start_total)
		{
			eprintf(V_ERR,"D* count mismatch for singles->start[%c] (%d) vs doubles[%c]->start_total (%d)\n", c.letters[i], l->singles->start[i].count, c.letters[i], l->doubles[i]->start_total);
			if (is_exact_multiple(l->singles->start[i].count, l->doubles[i]->start_total))
			{
				if ((l->singles->start[i].count > l->doubles[i]->start_total) && l->doubles[i]->start_total)
				{
					u32 c_factor = l->singles->start[i].count / l->doubles[i]->start_total;
					eprintf(V_ERR,"D* fixing table by factor of %d:\n", c_factor);
					// iterate through the table and correct the numerators
					for (u32 j = 0; j < l->num_letters; j++)
					{
						l->doubles[i]->start[j].count *= c_factor;
					}
					// correct the denominator
					l->doubles[i]->start_total = l->singles->start[i].count;
				}
				else
					eprintf(V_ERR,"D* cannot fix.\n");
			}
			else
				eprintf(V_ERR,"D* cannot fix due to lack of common factor.\n");
		}
	}

	// next heuristic: if the doubles[i]->start[*] count for letter * doesn't equal the denominator for triples[*][i]->start_total but is off by some factor, increase the latter to match
	for (u32 i = 0; i < l->num_letters; i++)
	{
		for (u32 j = 0; j < l->num_letters; j++)
		{
			if (l->doubles[i]->start[j].count != l->triples[i][j]->start_total)
			{
				eprintf(V_ERR,"D* count mismatch for doubles[%c]->start[%c] (%d) vs triples[%c][%c]->start_total (%d)\n", c.letters[i], c.letters[j], l->doubles[i]->start[j].count, c.letters[i], c.letters[j], l->triples[i][j]->start_total);
				if (is_exact_multiple(l->doubles[i]->start[j].count, l->triples[i][j]->start_total))
				{
					if ((l->doubles[i]->start[j].count > l->triples[i][j]->start_total) && l->triples[i][j]->start_total)
					{
						u32 c_factor = l->doubles[i]->start[j].count / l->triples[i][j]->start_total;
						eprintf(V_ERR,"D* fixing table by factor of %d:\n", c_factor);
						// iterate through the table and correct the numerators
						for (u32 k = 0; k < l->num_letters; k++)
						{
							l->triples[i][j]->start[k].count *= c_factor;
						}
						// correct the denominator
						l->triples[i][j]->start_total = l->doubles[i]->start[j].count;
					}
					else
						eprintf(V_ERR,"D* cannot fix.\n");
				}
				else
					eprintf(V_ERR,"D* cannot fix due to lack of common factor.\n");
			}
		}
	}

	// future heuristic: if the 'singles->end[n].count' is nonzero, and there is only one name for all m, where 'doubles[m]->end[n].count' is nonzero, propagate the count from 'singles->end[n].count' to 'doubles[m]->end[n].count'. Note that if the .ltr file is internally inconsistent due to hand-editing (the humanm and humanf files are definitely like this) do not attempt to do this, as it can propagate incorrect values!
	// an interesting caveat of this is the denominator can be propagated this way, which will adjust the numerator and denominator of many other values in the end tables, if successful.
	for (u32 i = 0; i < l->num_letters; i++)
	{
		u32 parents = 0;
		u32 pidx = 0;
		// if singles->end[i] is 0, don't bother.
		if (!(l->singles->end[i].count))
			continue;
		// if doubles[j]->end[i] buckets is not exactly 1, don't bother.
		//if (!(l->doubles[j]->end_buckets == 1))
		//	continue;
		// which bucket is 1?
		for (u32 j = 0; j < l->num_letters; j++)
		{
			if (l->doubles[j]->end[i].count)
			{
				//eprintf(V_ERR,"D* l->doubles[%c]->end[%c].count is %d\n", c.letters[j], c.letters[i], l->doubles[j]->end[i].count);
				parents++;
				pidx = j;
			}
		}
		if (parents == 1)
		{
			if (l->doubles[pidx]->end[i].count == l->singles->end[i].count)
				eprintf(V_ERR,"D* found exactly one parent (doubles[%c]->end[%c], count of %d out of %d) of singles->end[%c] (count of %d out of %d), but the counts already match, so we don't need to do anything here.\n", c.letters[pidx], c.letters[i], l->doubles[pidx]->end[i].count, l->doubles[pidx]->end_total, c.letters[i], l->singles->end[i].count, l->singles->end_total);
			else
			{
				eprintf(V_ERR,"D* found exactly one parent (doubles[%c]->end[%c], count of %d out of %d) of singles->end[%c] (count of %d out of %d), this may be a candidate for migration.\n", c.letters[pidx], c.letters[i], l->doubles[pidx]->end[i].count, l->doubles[pidx]->end_total, c.letters[i], l->singles->end[i].count, l->singles->end_total);
				// attempt to migrate
				//write me!
				if (is_exact_multiple(l->doubles[pidx]->end[i].count,l->singles->end[i].count))
				{
					eprintf(V_ERR,"D* factors are compatible, attempting migration.\n");
					eprintf(V_ERR,"D* ... or would, if this was written yet.\n");
					for (u32 m = 0; m < l->num_letters; m++)
					{
						// not done yet
					}
				}
			}
		}
	}
}

void cdf_print(cdf_array* p, u8 num_letters, u8 k, u8 j, u8 num, s_cfg c)
{
	u8 x = ' ';
	u8 y = ' ';
	u8 z = ' ';
	for (u8 i = 0; i < num_letters; i++) {
		// formatting
		if (num == 0)
		{
			x = c.letters[i];
		}
		else if (num == 1)
		{
			x = c.letters[j];
			y = c.letters[i];
		}
		else // num == 2
		{
			x = c.letters[k];
			y = c.letters[j];
			z = c.letters[i];
		}
		if ((c.printcdf == 2) || !((p->start[i].cdf_data == 0.0) && (p->middle[i].cdf_data == 0.0) && (p->end[i].cdf_data == 0.0)))
		{
			/*
			printf("%c%c%c      |% .5f    % .5f  |% .5f     % .5f   |% .5f  % .5f\n",
				x, y, z,
				p->start[i].cdf_data, p->start[i].pdf_data,
				p->middle[i].cdf_data, p->middle[i].pdf_data,
				p->end[i].cdf_data, p->end[i].pdf_data);
			*/
			printf("%c%c%c      |% .5f %5d /%5d |% .5f   %5d /%5d |% .5f %5d /%5d\n",
				x, y, z,
				p->start[i].cdf_data, p->start[i].count, p->start_total,
				p->middle[i].cdf_data, p->middle[i].count,p->middle_total,
				p->end[i].cdf_data, p->end[i].count,p->end_total);
		}
	}
}

void ltr_print(ltrfile* l, s_cfg c)
{
	if (!c.printcdf) return;
	printf("Number of letters in LTR: %d\n", l->num_letters);
	printf("Sequence | CDF(start)  P(start) | CDF(middle)  P(middle) | CDF(end)  P(end)\n");
	cdf_print(l->singles,l->num_letters,' ',' ',0,c);
	for (u8 j = 0; j < l->num_letters; j++)
	{
		cdf_print(l->doubles[j],l->num_letters,' ',j,1,c);
	}
	for (u8 k = 0; k < l->num_letters; k++)
	{
		for (u8 j = 0; j < l->num_letters; j++)
		{
			cdf_print(l->triples[k][j],l->num_letters,k,j,2,c);
		}
	}
}

void ltr_dumpstart(ltrfile* l, s_cfg c)
{
	if (!c.dumpstart) return;
	printf("D* %d names total:\n", l->singles->start_total);
	for (u32 k = 0; k < l->num_letters; k++)
	{
		for (u32 j = 0; j < l->num_letters; j++)
		{
			for (u32 i = 0; i < l->num_letters; i++)
			{
				for (u32 h = l->triples[k][j]->start[i].count; h > 0; h--)
				{
					printf("%c%c%c\n", c.letters[k], c.letters[j], c.letters[i]);
				}
			}
		}
	}
}

/*
		if (l->singles->start[i].pdf_data != 0.0)
		{
			for (u32 m = 0; m < (l->singles->start[i].count * singles_m); m++)
			{
				printf("%c", c.letters[i]); // first letter
				// deeper...
				u32 doubles_m = (l->singles->start[i].count / l->doubles[i]->start_total);
				for (u32 j = 0; j < l->num_letters; j++)
				{
					if (l->doubles[i]->start[j].pdf_data != 0.0)
					{
						for (u32 n = 0; n < (l->doubles[i]->start[j].count * doubles_m); n++)
						{
							printf("%c", c.letters[j]); // second letter
						}
					}
				}
				printf("\n");
			}
		}
*/

u8 l2offset(u8 in)
{
	u8 ret = in - 'a';
	if (in == '\'') return 26;
	if (in == '-') return 27;
	return ret;
}

void ltr_generate(ltrfile* l, s_cfg c) // generate exactly one name.
{

	char name[64] = {0};
	u32 index = 0;
	bool done = false;
	bool begin = true;
	float rng = 0.0;
	u8 i, j, k;
	s32 failcnt = 0;
	eprintf(V_GEN2,"D* generating name...\n");
	while (!done) // if we're not done yet
	{
		eprintf(V_GEN,"D* Current name state is \"%s\"\n", name);
		// generate the first 3 letters
		if (begin)
		{
			// initialze some variables here
			failcnt = 0;
			index = 0;
			do
			{
				// roll for a starting letter
				for (i = 0, rng = nrand(); i < l->num_letters; i++)
				{
					if (rng < l->singles->start[i].cdf_data)
						break;
				}

				if (i >= l->num_letters) // sanity check
					continue;

				// roll for the second letter
				for (j = 0, rng = nrand(); j < l->num_letters; j++)
				{
					if (rng < l->doubles[i]->start[j].cdf_data)
						break;
				}

				if (j >= l->num_letters) // sanity check
					continue;

				// roll for the third letter
				for (k = 0, rng = nrand(); k < l->num_letters; k++)
				{
					if (rng < l->triples[i][j]->start[k].cdf_data)
						break;
				}

			} while ((i >= l->num_letters) || (j >= l->num_letters) || (k >= l->num_letters)); // sanity check and loop condition in one

			// we did it! shove these 3 letters into a string
			name[index++] = c.letters[i];
			name[index++] = c.letters[j];
			name[index++] = c.letters[k];
			eprintf(V_GEN,"D* generated 3 first characters %c%c%c\n", c.letters[i], c.letters[j], c.letters[k]);
			begin = false;
		}
		// at this point index is at least 3.

		// make sure k was sane before shifting stuff over
		if (k < l->num_letters)
		{
			i = j;
			j = k;
		}

		// roll for another letter for k but don't use it yet
		rng = nrand();

		// roll to see whether the name ends here; names can't be longer than 12+1 letters and should be biased toward shorter names
		if ( (ms_rand() % 12) <= index ) // did our name end?
		{
			for (k = 0; k < l->num_letters; k++)
			{
				if (rng < l->triples[i][j]->end[k].cdf_data) // use the previous letter roll to find an ending triple
				{
					done = true; // no more letters needed, we just use the ending triple we found directly.
					// note there may be an original bug here, if k from this roll wasn't sane, we end abruptly?
					eprintf(V_GEN,"D* rolled to end the name after the next letter\n");
					break;
				}
			}
		}

		if (!done) // if we're not done yet, we still need more letters.
		{
			for (k = 0; k < l->num_letters; k++)
			{
				if (rng < l->triples[i][j]->middle[k].cdf_data) // use the previous letter roll to find an middle triple
					break;
			}
		}

		if (k < l->num_letters) // our roll was sane?
		{
			name[index++] = c.letters[k];
			eprintf(V_GEN2,"D* generated another character %c\n", c.letters[k]);
		}
		else if ((index > 3) && (failcnt < 100)) // no, it wasn't. we may be stuck in an impossible situation, so back up and try again
		{
			eprintf(V_GEN2,"D* backing up 1 character after failing a roll\n");
			// regenerate the old values for i and j
			j = l2offset(name[index-2]);
			i = l2offset(name[index-3]);
			name[index-1] = '\0'; // DEBUG: nuke the character at index-1
			index--;
			failcnt++;
		}
		else // we're definitely stuck in a bad way. just start over.
		{
			eprintf(V_GEN,"D* giving up and starting over\n");
			index = 0; // DEBUG: set index to 0
			begin = true;
		}
	}
	name[index] = '\0'; // add a trailing null
	// capitalize the first letter if it is a-z, leave it alone if it is - or '
	name[0] = toupper(name[0]);
	eprintf(V_GEN2,"D* generated name: %s\n", name);
	printf("%s\n", name);
}

void usage()
{
	printf("Usage: nwn_getname [options] file.ltr\n");
	printf("Generate one or more names from an .ltr file\n");
	printf("Optional parameters:\n");
	printf("-p\t: print the non-zero rows of the ltr CDF tables\n");
	printf("-pp\t: print all the rows of the ltr CDF tables\n");
	printf("-g #\t: generate # names (Default: 100)\n");
	printf("-s #\t: use # as the seed (Default: random)\n");
	printf("-f\t: if the ltr file has corrupt singles tables, do not fix them\n");
	printf("-d\t: dump the starting letters of every possible name\n");
	printf("-v #\t: verbose bitmask:\n");
	printf("\tParams/Seed  1\n");
	printf("\tGeneration   2\n");
	printf("\t  Details    4\n");
	printf("\tFix (default)8\n");
	printf("\t  Details    16\n");
	printf("\tLtr Load     32\n");
	printf("\t  Details    64\n");
	printf("\tLtr Free     128\n");
}

#define MIN_PARAMETERS 1
int main(int argc, char **argv)
{
	// defaults
	s_cfg c =
	{
		"abcdefghijklmnopqrstuvwxyz'-" // letters
		, 0 // printcdf
		, 100 // generate
		, time(NULL) // seed
		, true // fix
		, 8 /*8 == V_FIX*/ // verbose
		, false // dumpstart
	};

	if (argc < MIN_PARAMETERS+1)
	{
		eprintf(V_ERR,"E* Incorrect number of parameters!\n");
		usage();
		return 1;
	}

// handle optional parameters
	u32 paramidx = 1;
	while (paramidx < (argc-1))
	{
		switch (*(argv[paramidx]++))
		{
			case 'p':
				c.printcdf++;
				if (c.printcdf > 2) { eprintf(V_ERR,"E* Too much verbosity specified for print cdf parameter!\n"); usage(); exit(1); }
				break;
			case 'g':
				paramidx++;
				if (paramidx == (argc-1)) { eprintf(V_ERR,"E* Too few arguments for -g parameter!\n"); usage(); exit(1); }
				if (!sscanf(argv[paramidx], "%d", &c.generate)) { eprintf(V_ERR,"E* Unable to parse argument for -g parameter!\n"); usage(); exit(1); }
				paramidx++;
				break;
			case 's':
				paramidx++;
				if (paramidx == (argc-1)) { eprintf(V_ERR,"E* Too few arguments for -s parameter!\n"); usage(); exit(1); }
				if (!sscanf(argv[paramidx], "%d", &c.seed)) { eprintf(V_ERR,"E* Unable to parse argument for -s parameter!\n"); usage(); exit(1); }
				paramidx++;
				break;
			case 'f':
				c.fix = false;
				break;
			case 'd':
				c.dumpstart = true;
				break;
			case 'v':
				paramidx++;
				if (paramidx == (argc-1)) { eprintf(V_ERR,"E* Too few arguments for -v parameter!\n"); usage(); exit(1); }
				if (!sscanf(argv[paramidx], "%d", &c.verbose)) { eprintf(V_ERR,"E* Unable to parse argument for -v parameter!\n"); usage(); exit(1); }
				paramidx++;
				break;
			case '\0':
				paramidx++;
				break;
			case '-':
				// skip this character.
				break;
			default:
				{ eprintf(V_ERR,"E* Invalid option!\n"); usage(); exit(1); }
				break;
		}
	}
	eprintf(V_PARAM,"D* Parameters: generate: %d, seed: %d, print cdf: %s\n", c.generate, c.seed, c.printcdf?((c.printcdf==2)?"full":"brief"):"no");

// input file
	FILE *in = fopen(argv[argc-1], "rb");
	if (!in)
	{
		eprintf(V_ERR,"E* Unable to open input file %s!\n", argv[argc-1]);
		return 1;
	}

	fseek(in, 0, SEEK_END);
	u32 len = ftell(in);
	rewind(in); //fseek(in, 0, SEEK_SET);

	// simple filesize sanity checks
#define MINFILESIZE (8+1+(sizeof(float)*((1*3)+(1*1*3)+(1*1*1*3))))
#define MAXFILESIZE (8+1+(sizeof(float)*((28*3)+(28*28*3)+(28*28*28*3))))
	if ((len < MINFILESIZE) || (len > MAXFILESIZE))
	{
		eprintf(V_ERR,"E* Input file size of %d is too %s!\n", len, (len < MINFILESIZE)?"small":"large");
		fclose(in);
		exit(1);
	}

	// seed it!
	ms_srand(c.seed);

	// load it!
	ltrfile* infile = ltr_load(in, len, c);
	fclose(in);

	// analyze it!
	ltr_analyze(infile, c);

	// print it!
	ltr_print(infile, c);

	// dump it!
	ltr_dumpstart(infile, c);

	// generate some names!
	for (u32 i = 0; i < c.generate; i++)
		ltr_generate(infile, c);

	// free it!
	ltr_free(infile, c);

	return 0;
}


