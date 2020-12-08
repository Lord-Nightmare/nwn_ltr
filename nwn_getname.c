// license:BSD-3-Clause
// copyright-holders:Jonathan Gevaryahu
// (C) 2020 Jonathan Gevaryahu AKA Lord Nightmare
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>

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
#undef DEBUG_LOAD
#undef DEBUG_FREE
#undef DEBUG_GEN
#undef DEBUG_GEN_VERBOSE
#undef DEBUG_PARAM

// struct definitions
typedef struct cdf_array
{
	float* start;
	float* middle;
	float* end;
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
	int printcdf;
	u32 generate;
	u32 seed;
} s_cfg;

/* msrand implementation for consistency */
#define MSRAND_MAX 0x7fff
static u32 state = 1;

u32 ms_rand()
{
	state = ((state*214013) + 2531011) % (1<<31);
	return (state >> 16)&0x7fff;
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

float* f_alloc(u32 count)
{
	float* f = malloc(count * sizeof(float));
	if (f == NULL)
	{
		fprintf(stderr,"E* Failure to allocate memory for array of size %d, aborting!\n", count); fflush(stderr);
		return NULL;
	}
	return f;
}

cdf_array* cdf_alloc(u32 count)
{
	cdf_array *c = malloc(sizeof(cdf_array));
	if (c == NULL)
	{
		fprintf(stderr,"E* Failure to allocate memory for array, aborting!\n"); fflush(stderr);
		exit(1);
	}
	c->start = f_alloc(count);
	c->middle = f_alloc(count);
	c->end = f_alloc(count);
	return c;
}

void cdf_free(cdf_array *c)
{
	free(c->end);
	free(c->middle);
	free(c->start);
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
	if ((l->num_letters < 1) || (l->num_letters > 28))
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
			l->singles->start[i] = fget_f(in); pos+=4;
		}
		for (u32 i = 0; i < l->num_letters; i++)
		{
			l->singles->middle[i] = fget_f(in); pos+=4;
		}
		for (u32 i = 0; i < l->num_letters; i++)
		{
			l->singles->end[i] = fget_f(in); pos+=4;
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
				l->doubles[j]->start[i] = fget_f(in); pos+=4;
			}
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->doubles[j]->middle[i] = fget_f(in); pos+=4;
			}
			for (u32 i = 0; i < l->num_letters; i++)
			{
				l->doubles[j]->end[i] = fget_f(in); pos+=4;
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
					l->triples[k][j]->start[i] = fget_f(in); pos+=4;
				}
				for (u32 i = 0; i < l->num_letters; i++)
				{
					l->triples[k][j]->middle[i] = fget_f(in); pos+=4;
				}
				for (u32 i = 0; i < l->num_letters; i++)
				{
					l->triples[k][j]->end[i] = fget_f(in); pos+=4;
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

void cdf_print(cdf_array* p, u8 num_letters, u8 k, u8 j, u8 num, s_cfg cfg)
{
	const char* const letters = "abcdefghijklmnopqrstuvwxyz'-";
	u8 a = ' ';
	u8 b = ' ';
	u8 c = ' ';
	float s = 0.0, m = 0.0, e = 0.0;
	for (u8 i = 0; i < num_letters; i++) {
		// formatting
		if (num == 0)
		{
			a = letters[i];
		}
		else if (num == 1)
		{
			a = letters[j];
			b = letters[i];
		}
		else // num == 2
		{
			a = letters[k];
			b = letters[j];
			c = letters[i];
		}
		if ((cfg.printcdf == 2) || !((p->start[i] == 0.0) && (p->middle[i] == 0.0) && (p->end[i] == 0.0)))
		{
			printf("%c%c%c      |% .5f    % .5f  |% .5f     % .5f   |% .5f  % .5f\n", a, b, c,
				p->start[i],  p->start[i]  == 0.0 ? 0.0 : p->start[i]  - s,
				p->middle[i], p->middle[i] == 0.0 ? 0.0 : p->middle[i] - m,
				p->end[i],    p->end[i]    == 0.0 ? 0.0 : p->end[i]    - e);
		}
		if (p->start[i]  > 0.0) s = p->start[i];
		if (p->middle[i] > 0.0) m = p->middle[i];
		if (p->end[i]    > 0.0) e = p->end[i];
	}
}

void ltr_print(ltrfile* l, s_cfg cfg)
{
	if (!cfg.printcdf) return;
	printf("Number of letters in LTR: %d\n", l->num_letters);
	printf("Sequence | CDF(start)  P(start) | CDF(middle)  P(middle) | CDF(end)  P(end)\n");
	cdf_print(l->singles,l->num_letters,' ',' ',0,cfg);
	for (u8 j = 0; j < l->num_letters; j++)
	{
		cdf_print(l->doubles[j],l->num_letters,' ',j,1,cfg);
	}
	for (u8 k = 0; k < l->num_letters; k++)
	{
		for (u8 j = 0; j < l->num_letters; j++)
		{
			cdf_print(l->triples[k][j],l->num_letters,k,j,2,cfg);
		}
	}
}

u8 l2offset(u8 in)
{
	u8 ret = in - 'a';
	if (in == '\'') return 26;
	if (in == '-') return 27;
	return ret;
}

void ltr_generate(ltrfile* l) // generate exactly one name.
{
	const char* const letters = "abcdefghijklmnopqrstuvwxyz'-";
	char name[64] = {0};
	u32 index = 0;
	bool done = false;
	bool begin = true;
	float rng = 0.0;
	u8 i, j, k;
	s32 failcnt = 0;
#ifdef DEBUG_GEN_VERBOSE
	fprintf(stdout,"D* generating name...\n"); fflush(stdout);
#endif
	while (!done) // if we're not done yet
	{
#ifdef DEBUG_GEN
	if (index != 0) fprintf(stdout,"D* Current name state is \"%s\"\n", name); fflush(stdout);
#endif
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
					if (rng < l->singles->start[i])
						break;
				}

				if (i >= l->num_letters) // sanity check
					continue;

				// roll for the second letter
				for (j = 0, rng = nrand(); j < l->num_letters; j++)
				{
					if (rng < l->doubles[i]->start[j])
						break;
				}

				if (j >= l->num_letters) // sanity check
					continue;

				// roll for the third letter
				for (k = 0, rng = nrand(); k < l->num_letters; k++)
				{
					if (rng < l->triples[i][j]->start[k])
						break;
				}

			} while ((i >= l->num_letters) || (j >= l->num_letters) || (k >= l->num_letters)); // sanity check and loop condition in one

			// we did it! shove these 3 letters into a string
			name[index++] = letters[i];
			name[index++] = letters[j];
			name[index++] = letters[k];
#ifdef DEBUG_GEN
			fprintf(stdout,"D* generated 3 first characters %c%c%c\n", letters[i], letters[j], letters[k]); fflush(stdout);
#endif
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
				if (rng < l->triples[i][j]->end[k]) // use the previous letter roll to find an ending triple
				{
					done = true; // no more letters needed, we just use the ending triple we found directly.
					// note there may be an original bug here, if k from this roll wasn't sane, we end abruptly?
#ifdef DEBUG_GEN
			fprintf(stdout,"D* rolled to end the name after the next letter\n"); fflush(stdout);
#endif
					break;
				}
			}
		}

		if (!done) // if we're not done yet, we still need more letters.
		{
			for (k = 0; k < l->num_letters; k++)
			{
				if (rng < l->triples[i][j]->middle[k]) // use the previous letter roll to find an middle triple
					break;
			}
		}

		if (k < l->num_letters) // our roll was sane?
		{
			name[index++] = letters[k];
#ifdef DEBUG_GEN_VERBOSE
			fprintf(stdout,"D* generated another character %c\n", letters[k]); fflush(stdout);
#endif
		}
		else if ((index > 3) && (failcnt < 100)) // no, it wasn't. we may be stuck in an impossible situation, so back up and try again
		{
#ifdef DEBUG_GEN_VERBOSE
			fprintf(stdout,"D* backing up 1 character after failing a roll\n"); fflush(stdout);
#endif
			// regenerate the old values for i and j
			j = l2offset(name[index-2]);
			i = l2offset(name[index-3]);
			name[index-1] = '\0'; // DEBUG: nuke the character at index-1
			index--;
			failcnt++;
		}
		else // we're definitely stuck in a bad way. just start over.
		{
#ifdef DEBUG_GEN
			fprintf(stdout,"D* giving up and starting over\n"); fflush(stdout);
#endif
			index = 0; // DEBUG: set index to 0
			begin = true;
		}
	}
	name[index] = '\0'; // add a trailing null
	// capitalize the first letter if it is a-z, leave it alone if it is - or '
	name[0] = toupper(name[0]);
	printf("%s\n", name);
}

void usage()
{
	printf("Usage: nwn_getname [options] file.ltr\n");
	printf("Generate one or more names from an .ltr file\n");
	printf("Optional arguments:\n");
	printf("-p\t: print the non-zero rows of the ltr CDF tables\n");
	printf("-pp\t: print all the rows of the ltr CDF tables\n");
	printf("-g #\t: generate # names (Default: 100)\n");
	printf("-s #\t: use # as the seed (Default: random)\n");
	//printf("-f\t: if the ltr file has corrupt singles tables, do not fix them\n");
}

#define MIN_PARAMETERS 1

int main(int argc, char **argv)
{
	// defaults
	s_cfg cfg;
	cfg.printcdf = 0;
	cfg.generate = 100;
	cfg.seed = time(NULL);

	if (argc < MIN_PARAMETERS+1)
	{
		fprintf(stderr,"E* Incorrect number of parameters!\n"); fflush(stderr);
		usage();
		return 1;
	}

// handle extra args
	if (argc != MIN_PARAMETERS+1) // if we have more than one parameter
	{
		u32 paramidx = 1;
		while (paramidx < (argc-1))
		{
			switch (*(argv[paramidx]++))
			{
				case 'p':
					cfg.printcdf++;
					if (cfg.printcdf > 2) { fprintf(stderr,"E* Too much verbosity specified for print cdf parameter!\n"); fflush(stderr); usage(); exit(1); }
					break;
				case 'g':
					paramidx++;
					if (paramidx == (argc-1)) { fprintf(stderr,"E* Too few arguments for generate number parameter!\n"); fflush(stderr); usage(); exit(1); }
					if (!sscanf(argv[paramidx], "%d", &cfg.generate)) { fprintf(stderr,"E* unable to parse argument for generate number parameter!\n"); fflush(stderr); usage(); exit(1); }
					paramidx++;
					break;
				case 's':
					paramidx++;
					if (paramidx == (argc-1)) { fprintf(stderr,"E* Too few arguments for seed parameter!\n"); fflush(stderr); usage(); exit(1); }
					if (!sscanf(argv[paramidx], "%d", &cfg.seed)) { fprintf(stderr,"E* unable to parse argument for seed parameter!\n"); fflush(stderr); usage(); exit(1); }
					paramidx++;
					break;
				case '\0':
					paramidx++;
					break;
				case '-':
					// skip this character.
					break;
				default:
					usage();
					exit(1);
					break;
			}
		}
	}
#ifdef DEBUG_PARAM
	fprintf(stderr,"D* Parameters: generate: %d, seed: %d, print cdf: %s\n", cfg.generate, cfg.seed, cfg.printcdf?((cfg.printcdf==2)?"full":"brief"):"no"); fflush(stderr);
#endif

// input file
	FILE *in = fopen(argv[argc-1], "rb");
	if (!in)
	{
		fprintf(stderr,"E* Unable to open input file %s!\n", argv[argc-1]); fflush(stderr);
		return 1;
	}

	fseek(in, 0, SEEK_END);
	u32 len = ftell(in);
	rewind(in); //fseek(in, 0, SEEK_SET);

	// simple filesize sanity checks
#define MINFILESIZE (8+1+(sizeof(float)*((1*3)+(1*1*3)+(1*1*1*3))))
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

	// seed it!
	ms_srand(cfg.seed);

	// load it!
	ltrfile* infile = ltr_load(in, len);
	fclose(in);

	// fix it!
	///TODO: ltr_fix(infile);

	// print it!
	ltr_print(infile, cfg);

	// generate some names!
	for (u32 i = 0; i < cfg.generate; i++)
		ltr_generate(infile);

	// free it!
	ltr_free(infile);

	return 0;
}


