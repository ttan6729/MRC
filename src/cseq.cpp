#include "cseq.h"


//KSEQ_INIT(gzFile, gzread)





seq_file *seq_open(const char *fn)
{
	printf("open file:%s\n",fn);
	seq_file *fp;
	gzFile f = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if ( f == 0 ) 
	{
		fprintf(stderr,"failed to open file:%s\n",fn);
		return 0;
	}
	fp = (seq_file *)calloc(1,sizeof(seq_file));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void seq_close(seq_file *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

c_seq *seq_init(seq_file *fp, int *n_, int *seq_len)
{
	kseq_t *ks = fp->ks;
	int m = 0, n = 0; 
	c_seq *seqs = 0;
	while (kseq_read(ks) >= 0)
	{
		c_seq *s;
		if (n >= m)
		{
			m = m ? m<<1 : 1<<28;
			seqs = (c_seq*)realloc(seqs, m * sizeof(c_seq));
		}
		s = &seqs[n];
		s->seq = strdup(ks->seq.s);
		if(*seq_len != ks->seq.l)
		{
			*seq_len = ks->seq.l;
			//printf("seq len:%d, k seq len:%d\n",seq_len,ks->seq.l);
			//fprintf(stderr, "Length of reads are different. The program can not compress it.\n");
			//exit(1);
		}
		n++;
	}
	if (n == 0 ) fp->is_eof = 1;
	*n_ = n;
//	kseq_destroy(ks);
	return seqs;
}

void seq_read(c_seq *&seqs, seq_file *fp, int *n_, int seq_len)
{
	int m, n;
	kseq_t *ks = fp->ks;
	n = *n_;
	m = 256;
	while (m <= n) {
		m <<= 1;
	}
	while(kseq_read(ks) >= 0)
	{
		c_seq * s;
		if (n >= m)
		{
			m = m? m <<1 :256;
		}
		s = &seqs[n];
		s->seq = strdup(ks->seq.s);
		//s->qscore = strdup(ks->qual.s); 
		if (seq_len != ks->seq.l)
		{
			fprintf(stderr, "Length of reads are different, failed to compress.\n");
			exit(1);
		}
		n++;
	}
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
}

int bseq_eof(seq_file *fp)
{
	return fp->is_eof;
}

int average_qs(const char *qs,int len)
{
	int sum = 0;
	int base_value = 33;
	for (int i = 0; i < len; i++)
		sum += (int)qs[i] - base_value;
	return sum/len;
}

// float average_qscore(c_seq *seq, int seq_len, int l)  //length for based with higher weight
// {
// 	float q = 0.0;
// 	int base_value = 33;
// 	if( seq_len < 2 * l)
// 	{
// 		for(int i = 0; i < seq_len; i++)
// 			q += ( (int)seq->qscore[i] - base_value );
// 	} 
// 	else 
// 	{
// 		for(int i = 0; i < l; i++)
// 			q += ( (int)seq->qscore[i] - base_value );
// 		for(int i = l; i < seq_len -l; i++)
// 			q += 0.75 * ( (int)seq->qscore[i] - base_value ); 
// 		for(int i = seq_len - l; i < seq_len; i++)
// 			q += ( (int)seq->qscore[i] - base_value );
// 	}
// 	q = q/seq_len;
// 	seq->quality = q;
// 	return q;
// }


