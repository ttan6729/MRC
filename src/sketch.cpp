#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <array>
#include "kvec.h"
#include "creads.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/*static inline uint64_t hash64_(uint64_t key)
{
	key = (~key + (key << 21)); // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)); // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)); // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}
*/
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */

//hash value of minimizer of given reads with quality score larger than threshold

uint32_t sketch_prefix(const char *str,int len,int k)
{
	uint32_t kmer = 0;
	uint64_t mask = (1<<k) - 1;
	for (int i = 0; i < k; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
	}
	return kmer;
}

uint32_t sketch_suffix(const char *str,int len,int k)
{
	uint32_t kmer = 0;
	uint64_t mask = (1<<k) - 1;
	for (int i = len-k; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
	}
	return kmer;

}


uint64_t mm_sketch_single(const char *str, int len,int k)
{

	uint64_t  mask = (1ULL<<2*k) - 1, min=(1ULL<<2*k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, pos = k;
	//assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			pos = i;
		}
	}
	//printf(", pos %d, value: %ld %s\n",pos, min, str+pos-k+1);
	return min;
}


int mm_sketch_pos(const char *str, int len,int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}

	return min_pos;

}

//return two values, first is hash value of minimizer, second is position
std::array<int, 2> mm_sketch_value_pos(const char *str, int len, int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}


	std::array<int, 2> result;
	result[0] = min;
	result[1] = min_pos;
	return result;
}

int mm_sketch_pos2(const char *str, int len,int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}
	if ( kmer != 1732 )
		min_pos = -1;

	return min_pos;

}