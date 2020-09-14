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

float threshold = 2.0;
int n_threads = 12;
int rw;
int n = 0,selected_number = 50,m;
float **process(char *dir, int k);
void kmeans_clustering(vector<string> _file_list,float **vectors, char *dir);
char *file_list_path;
vector<string> file_list;

int main (int argc, char **argv)
{	
	int k = 15;
	int c;
	char *dir;
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
        case 's':
        selected_number = atoi(optarg);
        break;
      	case 'e':
      	threshold = atof(optarg);
      	break;
    }	
	float **vectors = process(dir,k);
	kmeans_clustering(file_list,vectors,dir);
	return 0;
}

float **process(char *dir,int k)
{
	printf("input file list path: %s\n",file_list_path);
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
    else
    {
    	printf("failed to open %s\n",file_list_path);
    	exit(1);
    }
	float **vectors = pre_process(file_list,dir,k,selected_number);


	return vectors;
}

void kmeans_clustering(vector<string> _file_list,float **vectors,char *dir)
{
	m = selected_number;
	string value,line;
	vector<string> *file_list = new vector<string>();

	for(int i = 0; i < n; i++)
		file_list->push_back(_file_list[i]);

	Kmeans *km = new Kmeans(vectors,n,m,dir,threshold);
	km->set_list(file_list);
	km->runKmeans();


	for(int i = 0; i < n; i++)
		free(vectors[i]);
	free(vectors);	
}	

