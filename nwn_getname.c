// license:BSD-3-Clause
// copyright-holders:Jonathan Gevaryahu
// (C) 1998-2019 Jonathan Gevaryahu AKA Lord Nightmare
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
#define DEBUG 1

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
		fprintf(stderr,"E* We don't handle big endian stuff yet! bailing out!\n"); fflush(stderr);
		exit(1);
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
		fprintf(stderr,"E* Failure to allocate memory for array of size %d, aborting!\n", f->len);
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
	c->start = f_alloc(count);
	c->middle = f_alloc(count);
	c->end = f_alloc(count);
}

void cdf_free(cdf_array *c)
{
	f_free(c->end);
	f_free(c->middle);
	f_free(c->start);
	free(c);
}

ltrfile loadltr(FILE *in, u32 len)
{
	ltrfile ret;
	u32 pos = 0;
	// magic
	{ // scope-limit
		char* compare = "LTR V1.0 ";
		for (u32 i = 0; i < 8; i++)
		{
			ret.magic[i] = fgetc(in);
			if (ret.magic[i] != compare[i])
			{
				fprintf(stderr,"E* Incorrect magic number! Exiting!\n"); fflush(stderr);
				fclose(in);
				exit(1);
			}
			pos++;
		}
	}
	// number of letters
	ret.num_letters = fgetc(in);
	if ((ret.num_letters < 3) || (ret.num_letters > 28))
	{
		fprintf(stderr,"E* Invalid number of letters %d! Exiting!\n", ret.num_letters); fflush(stderr);
		fclose(in);
		exit(1);
	}
	pos++;
#ifdef DEBUG
	fprintf(stderr,"D* LTR header read ok, num_letters = %d\n", ret.num_letters); fflush(stderr);
#endif

	// allocate singles table
	{ // scope-limit
		ret.singles = cdf_alloc(ret.num_letters);
		if (ret.singles == NULL)
		{
			fclose(in);
			exit(1);
		}
	}
#ifdef DEBUG
	fprintf(stderr,"D* successfully allocated the singles cdf table\n"); fflush(stderr);
#endif
	// fill singles table
	{ // scope-limit
		for (u32 i = 0; i < ret.num_letters; i++)
		{
			//ret.singles->start->data[i] = fget_f(in);
			fprintf(stdout,"debug: seen value of %f at offset %d\n", fget_f(in), i);
		}
	}
	//todo: write me!
	return ret;
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

	// encapsulation woo
	ltrfile infile = loadltr(in, len);

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


