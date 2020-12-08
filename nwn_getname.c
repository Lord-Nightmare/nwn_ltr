// license:BSD-3-Clause
// copyright-holders:Jonathan Gevaryahu
// (C) 2020 Jonathan Gevaryahu AKA Lord Nightmare
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

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
#define DEBUG_LOAD 1
#define DEBUG_FREE 1

// struct definitions

/*
typedef enum {
	ZERO = 0,
	ONE,
	TWO
} number_type_t;
*/
/*
// biquad filter params used to generate the filter constants a1/2, b0/1/2
typedef struct biquad_params
{
	biquad_type_t type;
	double fc;
	double q;
	double gain;
	uint32_t sr;
} biquad_params;

// generic biquad filter
typedef struct biquad_filter
{
	biquad_params params;
	double w0;
	double w1;
	double w2;
	// double a0; // only used if gain is done on the input side, which requires an extra multiply, so we're not doing it here.
	double a1;
	double a2;
	double b0;
	double b1;
	double b2;
} biquad_filter;
*/

typedef struct f_array
{
	uint32_t len;
	float* data;
} f_array;

typedef struct cdf_array
{
	f_array* start;
	f_array* middle;
	f_array* end;
} cdf_array;

typedef struct ltrfile
{
	char magic[8];
	uint8_t num_letters;
	cdf_array* singles;
	cdf_array** doubles;
	cdf_array*** triples;
} ltrfile;

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
	else // big endian
	{
		t = 0;
		for (u32 i = 0; i < 4; i++)
		{
			t <<= 8;
			t |= fgetc(in);
		}
		fprintf(stderr,"W* We may not handle big endian stuff properly yet!\n"); fflush(stderr);
		//exit(1);
	}
	return *(float *) &t;
}

f_array* f_alloc(u32 count)
{
	f_array *f = malloc(sizeof(f_array));
	f->len = count;
	f->data = (float *) malloc(count * sizeof(float));
	if (f->data == NULL)
	{
		fprintf(stderr,"E* Failure to allocate memory for array of size %d, aborting!\n", f->len); fflush(stderr);
		return NULL;
	}
	return f;
}

void f_free(f_array *f)
{
	free(f->data);
	f->len = 0;
	free(f);
}

cdf_array* cdf_alloc(u32 count)
{
	cdf_array *c = malloc(sizeof(cdf_array));
	if (c == NULL)
	{
		fprintf(stderr,"E* Failure to allocate memory for array, aborting!\n"); fflush(stderr);
		exit(1);
		//return NULL;
	}
	c->start = f_alloc(count);
	c->middle = f_alloc(count);
	c->end = f_alloc(count);
	return c;
}

void cdf_free(cdf_array *c)
{
	f_free(c->end);
	f_free(c->middle);
	f_free(c->start);
	free(c);
}

ltrfile* ltr_load(FILE *in, u32 len)
{
	ltrfile *l = malloc(sizeof(ltrfile));
	u32 pos = 0;
	// magic
	{ // scope-limit
		char* compare = "LTR V1.0 ";
		for (u32 i = 0; i < 8; i++)
		{
			l->magic[i] = fgetc(in);
			if (l->magic[i] != compare[i])
			{
				fprintf(stderr,"E* Incorrect magic number! Exiting!\n"); fflush(stderr);
				fclose(in);
				exit(1);
			}
			pos++;
		}
	}
	// number of letters
	l->num_letters = fgetc(in);
	if ((l->num_letters < 3) || (l->num_letters > 28))
	{
		fprintf(stderr,"E* Invalid number of letters %d! Exiting!\n", l->num_letters); fflush(stderr);
		fclose(in);
		exit(1);
	}
	pos++;
#ifdef DEBUG_LOAD
	fprintf(stderr,"D* LTR header read ok, num_letters = %d\n", l->num_letters); fflush(stderr);
#endif

	// allocate and fill singles table
	{ // scope-limit
		l->singles = cdf_alloc(l->num_letters);
		if (l->singles == NULL)
		{
			fclose(in);
			exit(1);
		}
#ifdef DEBUG_LOAD
		fprintf(stderr,"D* successfully allocated the singles cdf table\n"); fflush(stderr);
#endif
		for (u32 i = 0; i < l->num_letters; i++)
		{
			l->singles->start->data[i] = fget_f(in); pos+=4;
		}
		for (u32 i = 0; i < l->num_letters; i++)
		{
			l->singles->middle->data[i] = fget_f(in); pos+=4;
		}
		for (u32 i = 0; i < l->num_letters; i++)
		{
			l->singles->end->data[i] = fget_f(in); pos+=4;
		}
#ifdef DEBUG_LOAD
		fprintf(stderr,"D* successfully filled the singles cdf table\n"); fflush(stderr);
#endif
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
#ifdef DEBUG_LOAD
			fprintf(stderr,"D* successfully allocated the doubles cdf table %d\n", j); fflush(stderr);
#endif
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->doubles[j]->start->data[i] = fget_f(in); pos+=4;
			}
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->doubles[j]->middle->data[i] = fget_f(in); pos+=4;
			}
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->doubles[j]->end->data[i] = fget_f(in); pos+=4;
			}
#ifdef DEBUG_LOAD
			fprintf(stderr,"D* successfully filled the doubles cdf table %d\n", j); fflush(stderr);
#endif
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
#ifdef DEBUG_LOAD
				fprintf(stderr,"D* successfully allocated the triples cdf table %d:%d\n", k, j); fflush(stderr);
#endif
				for (u32 i = 0; i < l->num_letters; i++)
				{
					l->triples[k][j]->start->data[i] = fget_f(in); pos+=4;
				}
				for (u32 i = 0; i < l->num_letters; i++)
				{
					l->triples[k][j]->middle->data[i] = fget_f(in); pos+=4;
				}
				for (u32 i = 0; i < l->num_letters; i++)
				{
					l->triples[k][j]->end->data[i] = fget_f(in); pos+=4;
				}
#ifdef DEBUG_LOAD
				fprintf(stderr,"D* successfully filled the triples cdf table %d:%d\n", k, j); fflush(stderr);
#endif
			}
		}
	}
#ifdef DEBUG_LOAD
	fprintf(stderr,"D* all tables loaded, expected size was %d, final size was %d\n", len, pos); fflush(stderr);
#endif
	return l;
}

void ltr_free(ltrfile* l)
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
#ifdef DEBUG_FREE
	fprintf(stderr,"D* everything is freed!\n"); fflush(stderr);
#endif
}

void ltr_print(ltrfile* l)
{
	const char* const letters = "abcdefghijklmnopqrstuvwxyz'-";

    printf("Num letters: %d\n", l->num_letters);
    printf("Sequence | CDF(start)  P(start) | CDF(middle)  P(middle) | CDF(end)  P(end)\n");

	//singles
	float s = 0.0, m = 0.0, e = 0.0;
	/*for (u8 i = 0; i < l->num_letters; i++) {
        struct cdf *p = &ltr->data.singles;
        printf("%c        |% .5f    % .5f  |% .5f     % .5f   |% .5f  % .5f\n", letters[i],
                p->start[i],  p->start[i]  == 0.0 ? 0.0 : p->start[i]  - s,
                p->middle[i], p->middle[i] == 0.0 ? 0.0 : p->middle[i] - m,
                p->end[i],    p->end[i]    == 0.0 ? 0.0 : p->end[i]    - e);

        if (p->start[i]  > 0.0) s = p->start[i];
        if (p->middle[i] > 0.0) m = p->middle[i];
        if (p->end[i]    > 0.0) e = p->end[i];*/
    }

}

void usage()
{
	printf("Usage: nwn_getname file.ltr [number of names]\n");
	printf("Generate one or more names from an .ltr file\n");
	printf("\n");
}

#define NUM_PARAMETERS 1

int main(int argc, char **argv)
{
	if ((argc != NUM_PARAMETERS+1) && (argc != NUM_PARAMETERS+2))
	{
		fprintf(stderr,"E* Incorrect number of parameters!\n"); fflush(stderr);
		usage();
		return 1;
	}

// input file
	FILE *in = fopen(argv[1], "rb");
	if (!in)
	{
		fprintf(stderr,"E* Unable to open input file %s!\n", argv[1]); fflush(stderr);
		return 1;
	}

	fseek(in, 0, SEEK_END);
	u32 len = ftell(in);
	rewind(in); //fseek(in, 0, SEEK_SET);

	// simple filesize sanity checks
#define MINFILESIZE (8+1+(sizeof(float)*((3*3)+(3*3*3)+(3*3*3*3))))
#define MAXFILESIZE (8+1+(sizeof(float)*((28*3)+(28*28*3)+(28*28*28*3))))
	if (len < MINFILESIZE)
	{
		fprintf(stderr,"E* input file size of %d is too small!\n", len);
		fclose(in);
		exit(1);
	}
	else if (len > MAXFILESIZE)
	{
		fprintf(stderr,"E* input file size of %d is too large!\n", len);
		fclose(in);
		exit(1);
	}

	// load it!
	ltrfile* infile = ltr_load(in, len);
	fclose(in);

	// print it!
	ltr_print(infile);

	// free it!
	ltr_free(infile);
/*


	{ // scope limiter for temp
		uint32_t temp = fread(dataArray, sizeof(uint8_t), len, in);
		fclose(in);
		if (temp != len)
		{
			fprintf(stderr,"E* Error reading in %d elements, only read in %d, aborting!\n", len, temp);
			free(dataArray);
			dataArray = NULL;
			return 1;
		}
		fprintf(stderr,"D* Successfully read in %d bytes\n", temp);
	}

// prepare output file
	FILE *out = fopen(argv[2], "wb");
	if (!out)
	{
		fprintf(stderr,"E* Unable to open output file %s!\n", argv[2]);
		free(dataArray);
		dataArray = NULL;
		return 1;
	}
	fflush(stderr);
*/



// actual program goes here

//	fclose(out); 
	//free(dataArray);
	//dataArray = NULL;
	return 0;
}


