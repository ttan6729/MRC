#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include "creads.h"
//#include "config.h"
#include "cseq.h"
#include "kmeans.h"
using namespace std;

char *folder;
void liftrlimit()  // for linux get and set resoruce limits
{
#ifdef __linux__
	struct rlimit r;    
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}
uint32_t num_locks = 0x1000000;


int diff_threshold = 4; //15 for ERR174310_1 10 for ERR174310
int maxthr; //
int cbthreshold;
int first_mininum = 6; // 4 6 10 20 100

int n_threads = 12;
int rw;
int n = 0,selected_number = 50,m;
float **process(char *dir, int k);
void kmeans_clustering(vector<string> _file_list,float **vectors, char *dir);
char *file_list_path;
vector<string> file_list;

int main (int argc, char **argv)
{	
	int k = 8;
	int c;
	char *dir;
	printf("begin para\n");
	while ((c = getopt(argc, argv, "r:o:t:k:e:")) != -1)	
    switch (c)
    {
      	case 'r':
        file_list_path = optarg;
        break;
     	case 'o':	
     	dir = optarg;
     	break;
     	case 't':
        n_threads = atoi(optarg);
        break;
      	case 'k':
        k = atoi(optarg);
        break;
      	case 'e':
      	break;
    }	
   	printf("end para\n");
	float **vectors = process(dir,k);
	kmeans_clustering(file_list,vectors,dir);
	return 0;
}

float **process(char *dir,int k)
{
	printf("test for MRC\n");
	printf("open file: %s\n",file_list_path);
	fstream file;

   	file.open(file_list_path,ios::in); 
   	if (file.is_open())
   	{   //checking whether the file is open
      string line;
      while(getline(file, line))
      { 
        file_list.push_back(line);
        n++;
      }
      file.close(); //close the file object.
    }
	float **vectors = pre_process(file_list,dir,k,selected_number);


	return vectors;
}

void kmeans_clustering(vector<string> _file_list,float **vectors,char *dir)
{
	m = selected_number;
    printf("n: %d  m: %d\n",n,m);
	string value,line;
	vector<string> *file_list = new vector<string>();

	for(int i = 0; i < n; i++)
		file_list->push_back(_file_list[i]);
		//getline(myfile,line);
	// for(int i = 0; i < n; i++)
	// {
	// 	getline(myfile,line,',');
	// 	file_list->push_back(line);
	// 	for(int j = 0; j < m-1; j++)
	// 	{
	// 		getline(myfile,value,',');
	// 		//printf("%d %d,%s\n",i,j,value.c_str());
	// 		vectors[i][j] = stof(value);
	// 	}
	// 	getline(myfile,value);
	// 	ss.str(value);
	// 	getline(ss,value,',');
	// 	//printf("%d %d,%s\n",i,m-1,value.c_str());
	// 	vectors[i][m-1] = stof(value);
	// }
	
	// for(int i = 0; i < n; i++)
	// {
	//  	float sum = 0.0;
	//  	for(int j = 0; j < m; j++)
 // 			sum += vectors[i][j];
	//  	for(int j = 0; j < m; j++)
	//  		vectors[i][j] = 100000*vectors[i][j]/sum;
	// }
	//myfile.close();
	Kmeans *km = new Kmeans(vectors,n,m,dir);
	km->set_list(file_list);
	printf("begin\n");
	km->runKmeans();


	for(int i = 0; i < n; i++)
		free(vectors[i]);
	free(vectors);	
}	

