/**
kmeans with cluster size limitation
**/
#include "kmeans.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
int count_split = 0;
ofstream fp;
int clustered_number = 0;

cluster *create_cluster(cluster_result *cr,int id)
{
	cluster *c = (cluster *)malloc(sizeof(cluster));
	kv_init(*c);
	for(int i = 0; i < cr->size[id]; i++)
		kv_push(int,*c,cr->ids[id][i]);	
	return c;
}

void Kmeans::runKmeans()
{
	// for(int i = 0; i < n; i++)
	// {
	// 	for(int j = i+1; j < n; j++)
	// 		printf("distance %d %d: %.4f\n",i,j,Distance(i,j));
	// }
	printf("size of file list:%d\n",file_list->size());
	char file_name[1024];
	sprintf(file_name,"%s/cluster",dir);
	fp.open(file_name);

	printf("m:%d n:%d\n",m,n);
	printf("check matrix\n");
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 10; j++)
			printf("%.4f ",vectors[i][j]);
		printf("\n");	
	}

	k = 2;
	single_points = (cluster *)malloc(sizeof(cluster));
	kv_init(*single_points);
	srand(0);
	cluster *c = (cluster *)malloc(sizeof(cluster));
	kv_init(*c);
	for(int i = 0; i < n; i++)
		kv_push(int,*c,i);
	split(c); //recursive split
	printf("single points: ");
	for(int i = 0; i < single_points->n; i++) 
		printf("%d ",single_points->a[i]);
	printf("\n");
	cluster *c1 = (cluster *)malloc(sizeof(cluster));
	kv_init(*c1);
	kv_copy(int,*c1,*single_points);
	kv_init(*single_points);
	printf("single point size:%d, c1 size %d\n",single_points->n,c1->n);

	printf("\n");	
	for(int i = 0; i < c1->n; i++)
		cout << c1->a[i] << " ";
	printf("\n");
	split(c1);
	printf("single point size:%d\n",single_points->n);
	for(int i = 0; i < single_points->n; i++)
		fp << single_points->a[i] << "\n";
	printf("clustered:% d\n,",clustered_number);
	fp.close();
}

string strip(string s)
{
	int pos = s.find_last_of("/") + 1;
	return s.substr(pos, s.length() - pos);
}

void Kmeans::split(cluster *c)
{
	if ( c->n > 3)
	{
		cluster_result *cr = Cluster(c);

		for(int i = 0; i < cr->n; i++)
		{
			split(create_cluster(cr,i));
		}
	}
	else if (c->n == 1)
	{
		printf("single result:%d\n", c->a[0]);
		kv_push(int,*single_points,c->a[0]);
		//fp << c->a[0] << "\n";
	}	
	else if ( c->n ==2 || c->n == 3)
	{
		clustered_number += c->n;
		//kv_push(cluster,results,*c);
		for(int i = 0; i < c->n-1; i++)
		{
			for(int j = i+1; j < c->n; j++)	
			{
				string s1 = strip((*file_list)[c->a[i]]),  s2 = strip((*file_list)[c->a[j]]);
				printf("distance between %d %d, %s %s: %.2f\n",c->a[i],c->a[j],s1.c_str(),s2.c_str(),Distance(c->a[i],c->a[j]));
			}
		}
		
		printf("result:");
		for(int i = 0; i < c->n; i++)
		{	
			fp << c->a[i] << " ";
			printf(" %d",c->a[i]);
		}
		fp << "\n";
		printf("\n");
	}
	else if (c->n == 0)
	{
		printf("error in cluster size\n");
	//	exit(0);
	}
}


cluster_result *Kmeans::Cluster(cluster *c)
{
	//printf("begin cluster\n");
	int max_iter = 50;
	vector<int> prev_centers, new_centers;

	while( prev_centers.size() < k )
	{	
		int center = rand()%(c->n);
		if( !count(prev_centers.begin(), prev_centers.end(), c->a[center]) )
			prev_centers.push_back(c->a[center]);
		
	}

	clusters *cs = (clusters *)malloc(sizeof(clusters));
	kv_init(*cs);
	for(int i = 0; i < k; i++)
	{
		cluster *_c = (cluster *)malloc(sizeof(cluster));
		kv_init(*_c);
		kv_push(cluster,*cs,*_c);
	}
	//printf("clustering, size %d, ",c->n);
	// printf("cluster:");
	// for(int i = 0; i < c->n; i++)
	// 	printf(" %d",c->a[i]);
	// printf("\n");
	for( int iter = 0; iter < max_iter; iter++ )
	{
		//printf("centers %d %d\n",prev_centers[0],prev_centers[1]);
		for(int i = 0; i < c->n; i++)
			kv_push(int,cs->a[NearestCenter(c->a[i],prev_centers)],c->a[i]);
		//printf("b1\n");
		new_centers = compute_centers(cs);
		if(new_centers == prev_centers)
			break;
		prev_centers = new_centers;	
		if( iter < max_iter-1)
		{
			for(int i = 0; i < cs->n; i++)
			{
				kv_destroy(cs->a[i]);
				kv_init(cs->a[i]);
			}
		}
	}

	cluster_result *result = (cluster_result*)malloc(sizeof(cluster_result));
	result->n = cs->n;
	result->size = (int *)malloc(sizeof(int)*n);
	result->ids = (int **)malloc(sizeof(int *) * cs->n);
	for(int i = 0; i < cs->n; i++)
	{
		cluster c = cs->a[i];
		result->ids[i] = (int *)malloc(sizeof(int)*c.n);
		result->size[i] = c.n;
		for(int j = 0; j < c.n; j++)
			result->ids[i][j] = c.a[j];
	}

	for(int i = 0; i < cs->n; i++)
	{
		if (cs->a[i].n ==1)
		{
			cluster_result *_result = rebalance(cs);
			printf("finish rebalance\n");
			if(_result->n != 0)
				result = _result;
			break;
		}
	}
	printf("return result, cluster size %d\n",result->n);

	// for(int i = 0; i < cs->n; i++)
	// {
	// 	printf("result cluster: ");
	// 	cluster c = cs->a[i];
	// 	for(int j = 0; j < c.n; j++)
	// 		printf(" %d",c.a[j]);
	// 	printf("\n");
	// }	
	return result;
}

//when exists cluster that size equals one
cluster_result *Kmeans::rebalance(clusters *cs)
{
	int point_id, cluster_id;
	for(int i = 0; i < cs->n; i++)
	{
		if(cs->a[i].n == 1)
		{ 
			point_id = cs->a[i].a[0];
			cluster_id = !i;
		}
	}

	printf("rebalance id: %d %d\n",point_id,cluster_id);


	vector<int> ids;
	for(int i = 0; i < cs->a[cluster_id].n; i++)
		ids.push_back(cs->a[cluster_id].a[i]);
	int min_id = ids[NearestCenter(point_id,ids)];
	float min_dis = Distance(point_id,min_id);
	
	cluster_result *result = (cluster_result*)malloc(sizeof(cluster_result));
	result->n  = 0;
	if (cs->a[cluster_id].n <=2)
		return result;

	if( min_dis <  2 * minDistance(min_id,ids) )
	{
		printf("add point, new result %d %d\n",min_id,point_id);
		result->n = 2;
		result->size = (int *)malloc(sizeof(int)*2);
		result->ids = (int **)malloc(sizeof(int *) * 2);
		result->ids[0] = (int *)malloc(sizeof(int)*2);
		result->size[0] = 2;
		result->size[1] = cs->a[cluster_id].n-1;
		result->ids[0][0] = point_id;
		result->ids[0][1] = min_id;
		result->ids[1] = (int *)malloc(sizeof(int)* (cs->a[cluster_id].n-1) );
		int count = 0;
		for(int i = 0; i < cs->a[cluster_id].n; i++)
		{
			if( cs->a[cluster_id].a[i] != min_id )
				result->ids[1][count++] = cs->a[cluster_id].a[i];
		}
		return result;
	}
	return result;
}


vector<int> Kmeans::compute_centers(clusters *cs)
{

	vector<int> centers; 
	for(int i = 0; i < cs->n; i++)
	{
		cluster c = cs->a[i];
		// printf("cluster %d:",i);
		// for(int j = 0; j < c.n; j++)
		//    printf(" %d",c.a[j]);
		// printf("\n");

		float minDistance = numeric_limits<float>::max();
		int minId = -1;
		for(int j = 0; j < c.n; j++)
		{
			minId = c.a[0];
			float distance = 0.0;
			for(int l = 0; l < c.n; l++)
				distance += Distance(c.a[j],c.a[l]);
			if( distance < minDistance )
			{
				minId = c.a[j];
				minDistance = distance;
			}
		}
		centers.push_back(minId);
	}

	return centers;
}




double Kmeans::minDistance(int id, vector<int> points)
{
	float minDistance = numeric_limits<float>::max();
	int k_id = -1;
	float dis;	
	for (int i = 0; i < points.size(); i++)
	{
		dis = Distance(id, points[i]);
		if (dis < minDistance && id != points[i])
		{
			minDistance = dis;
			k_id = i;
		}
	}
	return minDistance;	
}

int Kmeans::NearestCenter(int id, vector<int> centers)
{
	float minDistance = numeric_limits<float>::max();
	int k_id = -1;
	float dis;	
	for (int i = 0; i < centers.size(); i++)
	{
		dis = Distance(id, centers[i]);
		if (dis < minDistance)
		{
			minDistance = dis;
			k_id = i;
		}
	}
	return k_id;
}
