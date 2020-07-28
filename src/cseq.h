#ifndef C_SEQ_H
#define C_SEQ_H

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "kseq.h"

//struct for sequence file(fastq) and sequence set
extern unsigned char seq_nt4_table[256];
KSEQ_INIT(gzFile, gzread)
typedef struct {
	int rid;//, n_num, m_m, *pos_n;//rid: index of reads; n_num: number of N char; m_m: array size of pos_n; pos_n: store positon of N
	// uint32_v n_pos;
	char *seq;
	//char *qscore;
	//int quality;
	//void *n_pos;
	// void *p;
} c_seq;

typedef struct
{
	int is_eof;
	gzFile fp;
	kseq_t *ks;
} seq_file;

seq_file *seq_open(const char *fn);
void seq_close(seq_file *fp);
c_seq *seq_init(seq_file *fp, int *n_, int *seq_len); //create and read first sequece
void seq_read(c_seq *&seqs, seq_file *fp, int *n_, int seq_len);
int seq_eof(seq_file *fp);

int average_qs(const char *qs, int len); //average quality score of given char * with len
float average_qscore(c_seq *seq, int seq_len, int l);  //length for based with higher weight

double cputime();
double realtime();
#endif
