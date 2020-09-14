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
pthread_mutex_t mutex;
vector<uint32_t> read_numbers; //number of Non-n reads
vector<uint64_t> feature_ids;
ofstream info_fp;
//int **m1; 
uint32_t *loca; //location of each index with id
uint_v *mi_data;
float **result_data;
Matrix<double> *m2;
//buckets **mi_data;
struct reads_t;
int file_id = 0, seq_len, n_nreads; //number of reads contain 'N'
int k, b = 5, file_num;
uint64_t bucket_size, mutex_size; //(1ULL << (2*b)), mask = (1ULL<<(2*b)) - 1;  
uint32_t data_size = 0, feature_num = 0;
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

// void add_minmizer(uint64_t id)
// {
// 	uint64_t lock = id&mask;
// 	pthread_mutex_lock(&mutex[lock]);
// 	for(int i = 0; i < mi_data[lock].n; i++)
// 	{
// 		if( id == mi_data[lock].a[i].id)
// 		{
// 			mi_data[lock].a[i].n+=1;
// 			pthread_mutex_unlock(&mutex[lock]);
// 			return;
// 		}
// 	}
// 	bucket new_bucket; //= (bucket *)malloc(sizeof(bucket));
// 	new_bucket.id = id;
//     new_bucket.n = 1;
// 	data_size++;
// 	kv_push(bucket,mi_data[lock],new_bucket);

// 	pthread_mutex_unlock(&mutex[lock]);
// 	return;
// } 

static void process_read(reads_t *t, uint32_t rid) 
{
	// if(rid <=1800)
	// 	printf("rid %d %d\n", rid, loca[506]);

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

	uint64_t lock = mm_sketch_single(cur_seq, seq_len,k);//mm_sketch_single(cur_seq, seq_len,b);
//	printf("a\n");
	if( k <= 9 )
	{
		__sync_fetch_and_add(&mi_data[file_id].a[lock],1);
		return;
	}


	if(!loca[lock])
	{
		pthread_mutex_lock(&mutex);
		if(!loca[lock])
		{
			feature_ids.push_back(lock);
			feature_num++;
			loca[lock] = feature_num;
			for(int i = 0; i < file_num; i++)
				kv_push(uint32_t,mi_data[i],0U);
		}
		pthread_mutex_unlock(&mutex);
	}

	__sync_fetch_and_add(&mi_data[file_id].a[loca[lock]-1],1);
	return;
}
//an array contains all dimension, sort after

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
	printf("read length:%d, seq number:%d,thread number:%d\n",length,n_seq,t.n_threads);
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
	file_num = file_list.size();
	k = _k;
	reads_t t;
	pthread_t *tid;

	char info_name[1024];
	sprintf(info_name,"%s/info",dir);
	info_fp.open(info_name,std::fstream::in | std::fstream::out | std::fstream::app);
	t.n_threads = n_threads;//n_threads;
	t.w = (reads_worker_t*)alloca(n_threads * sizeof(reads_worker_t));

	bucket_size = (1ULL<<(2*k));
	loca = (uint32_t *)malloc(sizeof(uint32_t)*bucket_size);
	//memset(loca,0,bucket_size*sizeof(loca[0])/sizeof(char));
	for(int i = 0; i <bucket_size; i++)
		loca[i] = 0;

	pthread_mutex_init(&mutex, 0);
	
	mi_data = (uint_v *)malloc(sizeof(uint_v)*file_list.size());
	for(int i = 0; i < file_list.size(); i++)
		kv_init(mi_data[i]);
	if(k<=9)
	{
		for(int i = 0; i < file_num; i++)
		{
			for(int j = 0; j < bucket_size; j++)
				kv_push(uint32_t,mi_data[i],0U);
		}
				
	}

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

	if(k <=9)
		feature_num = bucket_size;

	vector<uint64_t> filiter_ids;
	for(int j = 0; j < feature_num; j++)
	{
		int count = 0;
		for(int i = 0; i < file_list.size(); i++)
			count += mi_data[i].a[j];
		if( count >= file_num/5 && count > 3)
			filiter_ids.push_back(j);
	}
	int filiter_num = filiter_ids.size();
//	printf("feature number: %ld, filitered: %d\n", feature_num,filiter_num);

//	feature_num /= 10;
	m2 = new Matrix<double>(file_list.size(),filiter_num);
	//normalization
	double fact = (double)read_numbers[0]*0.01/(mi_data[0].a[0]);
	for(int i = 0; i < file_list.size(); i++)
	{
		for(int j = 0; j < filiter_num; j++)
			m2->data[i][j] = fact*mi_data[i].a[filiter_ids[j]]/read_numbers[i] ;
	}

	// printf("test matrix\n");
	// for(int i = 0; i < 5; i++)
	// {
	// 	for(int j = 0; j < 10; j++)
	// 		printf("%.4f ",m2->data[i][j] );
	// 	printf("\n");
	// }



	vector<double> vars;
	for(int i = 0; i <filiter_num; i++)
		vars.push_back(m2->col_var(i));
	vector<double> vars2(vars);
	sort(vars2.begin(), vars2.end());

	result_data = (float**)malloc(sizeof(float *) * file_list.size());
	for(int i = 0; i < file_list.size(); i++)
		result_data[i] = (float *)malloc(sizeof(float) *selected_number);


	// ofstream vector_fp;
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
	normalize(result_data,file_list.size(),selected_number);
	// printf("test matrix\n");
	// for(int i = 0; i < 5; i++)
	// {
	// 	for(int j = 0; j < 10; j++)
	// 		printf("%.4f ",m2->data[i][j] );
	// 	printf("\n");
	// }

	free(mi_data);

	//delete m2;
	return result_data;
}

void normalize(float **matrix, int row, int column)
{
	for(int j = 0; j < column; j++)
	{
		float max = 0.0, min = numeric_limits<float>::max();
	 	for(int i = 0; i < row; i++)
	 	{
	 		if(max < matrix[i][j])
	 			max = matrix[i][j];
	 		if(min > matrix[i][j])
	 			min = matrix[i][j];
	 	}
	 	if( min < numeric_limits<float>::max() && (max-min)>0.0 )
	 	{
	 		float deno = max-min;
	 		for(int i = 0; i < row; i++)
	 			matrix[i][j] = (matrix[i][j]-min)/deno;
	 	}
	}

}