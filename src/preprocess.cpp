#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "kvec.h"
#include "creads.h"
#include "cseq.h"
#include "khash.h"
#include <unistd.h>
#include <ios>



using namespace std;

vector<uint32_t> read_numbers; //number of Non-n reads
vector<int> feature_ids;
ofstream info_fp;
int **m1; 
float **result_data;
Matrix<double> *m2;
uint32_t *mi_data, *prefix_data, *suffix_data;
struct reads_t;
int file_id = 0, num_bucket, seq_len;
int n_nreads; //number of reads contain 'N'
int k,b;
typedef struct {
	struct reads_t *t;
	long i;
} reads_worker_t; //thread worker

typedef struct reads_t {
	int n_threads;
	int seq_len;
	long n;
	reads_worker_t *w;
	c_seq *seqs;

} reads_t;  //reads for thread work

void mem_usage(double& vm_usage, double& resident_set) {
   vm_usage = 0.0;
   resident_set = 0.0;
   ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc
   //create some variables to get info
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; 
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; //
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

static inline long steal_reads_work(reads_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

int element_pos(vector<int> v, int ele)
{
	auto itr = find(v.begin(),v.end(),ele);
	if(itr != v.end())
		return distance(v.begin(),itr);
	return -1;	
}

string strip_slash(const char *c)
{
	string s(c);
	int pos = s.find_last_of("/");
	return s.substr(pos+1);
}

static void process_read(reads_t *t, uint32_t rid) 
{

	uint32_t mask = (1<<b) - 1,lock=0; //b=14 by default

	c_seq *seq = &t->seqs[rid];
	char *cur_seq = seq->seq;
	// fprintf(stderr, "rid: %ld\n", rid);

	int cntN = 0,cntA = 0, cntT = 0, cntG = 0, cntC = 0;
	for (int i = 0; i < seq_len; ++i) 
	{
		if (cur_seq[i] == 'N')
		{
			__sync_fetch_and_add(&n_nreads,1);
			return;
		}
	}


	lock = mm_sketch_single(cur_seq, seq_len,k);//mm_sketch_single(cur_seq, seq_len,b);
	__sync_fetch_and_add(&m1[file_id][lock],1);
	return;
}

static void *reads_worker(void *data) 
{
	reads_worker_t *w = (reads_worker_t*)data;
	reads_t *t = w->t;
	long i;
	for (;;) 
	{
		i = __sync_fetch_and_add(&w->i, t->n_threads);
		if (i >= t->n) break;
		process_read(t, i);
	}
	while ((i = steal_reads_work(t)) >= 0)
		process_read(t, i);
	pthread_exit(0);
}

void process_file(const char *file_name, int length,reads_t t, int n_threads)
{

	int i;
	seq_file *fp = seq_open(file_name);
	//if (fp == 0) return; 
	int n_seq = 0;
	c_seq *seqs = seq_init(fp,&n_seq,&length);

	info_fp << strip_slash(file_name) << " " << n_seq << "\n";
	t.seqs = seqs;
	t.n = n_seq;
	n_nreads = 0;
//	printf("read length:%d, seq number:%d,thread number:%d\n",length,n_seq,t.n_threads);
	pthread_t *tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = i; //
	// for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i)
		pthread_create(&tid[i], 0, reads_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], NULL);

	read_numbers.push_back(n_seq - n_nreads);
	for(int i = 0; i < n_seq; i++)
		free(seqs[i].seq);
	free(seqs);
	seq_close(fp);

	return;	
}


void init_file()
{

}
template <typename T>
void print_data( T **data, int n, int m)
{
	for( int i = 0; i < n; i++ )
	{
		for( int j = 0; j < m; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}

}

float **pre_process(vector<string> file_list,char *dir, int _k, int selected_number)
{
	k = _k;
	b = 2 * k;
	reads_t t;
	pthread_t *tid;

	char info_name[1024];
	sprintf(info_name,"%s/info",dir);
	info_fp.open(info_name,std::fstream::in | std::fstream::out | std::fstream::app);
	t.n_threads = n_threads;//n_threads;
	t.w = (reads_worker_t*)alloca(n_threads * sizeof(reads_worker_t));

    num_bucket = 1<<b;
	m1 = (int**)malloc(sizeof(int *) * file_list.size());
	for(int i = 0; i < file_list.size(); i++)
	{
		m1[i] = (int *)malloc(sizeof(int) *num_bucket);
		memset(m1[i],0,num_bucket);
	}
	m2 = new Matrix<double>(file_list.size(),num_bucket);
	seq_len = 0;
	for(int i = 0; i < file_list.size(); i++)
	{
		fstream file;
	   	file.open(file_list[i],ios::in); 
	   	if (file.is_open())
	   	{   //checking whether the file is open
	      string line;
	      getline(file, line);
	      getline(file, line);
	      seq_len = strlen(line.c_str());
	      if (i == 0)
	      {
	      	t.seq_len = seq_len;
	      }
	      else if ( t.seq_len != seq_len)
	      {
	      	printf("error, the length of dataset doesn't equal\n");
	      	exit(0);
	      }
	      file.close(); //close the file object.
	      process_file(file_list[i].c_str(),seq_len,t,n_threads);
	      file_id++;
	    }
	    else 
	    {
	    	printf("failed to open %s\n",file_list[i].c_str());
	    } 
	}
	//normalization
	for(int i = 0; i < file_list.size(); i++)
	{
		for(int j = 0; j < num_bucket; j++)
			m2->data[i][j] = 1000.0*m1[i][j]/read_numbers[i];
	}



	vector<double> vars;
	for(int i = 0; i <num_bucket; i++)
		vars.push_back(m2->col_var(i));
	vector<double> vars2(vars);
	sort(vars2.begin(), vars2.end());

	result_data = (float**)malloc(sizeof(float *) * file_list.size());
	for(int i = 0; i < file_list.size(); i++)
		result_data[i] = (float *)malloc(sizeof(float) *selected_number);

	ofstream vector_fp;

	for(int i = 0; i < selected_number; i++)
	{

		double value = vars2[vars2.size()-i-1];
		auto itr = find(vars.begin(),vars.end(),value);
		int index = distance(vars.begin(),itr);
		feature_ids.push_back(index);
	//	printf("feature id:%d,vars %.2f,vars check %.2f\n",index,value,m2->col_var(index));
		for(int j = 0; j < file_list.size(); j++)
		{
			result_data[j][i] = m2->data[j][index];
		}

	}

	free(mi_data);

	delete m2;
	return result_data;
}