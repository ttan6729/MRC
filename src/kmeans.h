#include <bitset>  
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>  
#include <assert.h>
#include <omp.h>
#include "kvec.h"


typedef struct { int n, m; int *a; } cluster; //n:element number, m:maximize number, a:element,point id  
typedef struct { int n, m; cluster *a; } clusters; //as output of kmeans
typedef struct { int n; int *size; int **ids; } cluster_result;

class Kmeans
{
private:
 	int n,m; //vector number, dimension of vector
	int k;
	float **vectors;
	cluster *single_points; 
	clusters results;
	std::vector<std::string> *file_list;
	char *dir;
public:
	Kmeans(float **_vectors,int _n, int _m, char *_dir): 
	vectors{_vectors}, n{_n}, m{_m}, dir{_dir} { }
	void split(cluster *c);
	cluster_result *Cluster(cluster *c);
	void runKmeans();
	int NearestCenter(int id, std::vector<int> centers);
	std::vector<int> compute_centers(clusters* cs);
	cluster_result *rebalance(clusters *cs);
	double minDistance(int id, std::vector<int> points);
	void set_list(std::vector<std::string> *v) { file_list = v; }
	float Distance(int id1, int id2)
	{
		float result = 0.0;
		for(int i = 0; i < m; i++)
		{
			float value = vectors[id1][i] - vectors[id2][i];
			if( value < 0)
				value = value*(-1);
			result += value;
		}
	//	printf("distance %d %d, %.4f\n",id1,id2,result);
		return result;
	}


};

