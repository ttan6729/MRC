#ifndef CReads_H
#define CReads_H

#include <stdbool.h>
#include <bitset>  //fixed-sze sequecne of N bits
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <array>
#include "khash.h"
#include "cseq.h"

//struct for reads processing and compression

template <typename T>
class Matrix {
    // implementation


public: // interface
    int col, row;
    T** data;
    Matrix(int r, int c): row(r), col(c) 
    {
    	data = (T**)malloc(sizeof(T*)*col);
    	data[0] = (T*)malloc(sizeof(T)*row*col);
        memset(data[0],0,row*col);
    	for(int i = 1; i < r; i++)
    		data[i] = data[0]+i*col;
    }

    ~Matrix()
    {
    	free(data[0]);
    	free(data);
    }
    //allow to use matrix[i][j]
    T* operator[](int i) {
        return data[i];
    }

    double col_var(int id)
    {
        double sum=0.0, sumsq = 0.0;
        for(int i = 0; i < row; i++)
        {
            sum += data[i][id];
            sumsq += data[i][id] *data[i][id];            
        }   
        return (sumsq - (sum*sum/row))/(row-1);
    }
    void max_min(int index, T *_min, T *_max)
    {
        T min = data[0][index], max =data[0][index];
        for(int i = 1; i < row;i++)
        {
            if( min > data[i][index])
                min = data[i][index];
            if( max < data[i][index])
                max = data[i][index];
        }
        (*_min) = min;
        (*_max) = max;
        return;
    }
};

typedef struct { size_t n; uint64_t id;} bucket;
typedef struct { size_t n,m; bucket *a;} buckets;
typedef struct { size_t n,m; uint *a;} uint_v;
//vector of reads I, n: element number, m: maximum number of elements
typedef struct { size_t n, m; uint64_t *a; } uint64_v; 
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; int32_t *a; } int32_v;
typedef struct { size_t n, m; char *a; } char_v;
typedef struct { size_t n, m; bool *a; } bool_v;
typedef struct { size_t n; uint8_t a; } uint8bit_v;


#define getmax(a,b) ((a)>(b)?(a):(b))
#define getmin(a,b) ((a)<(b)?(a):(b))
#define getmax3(a,b,c) getmax(getmax(a,b),c)
#define getmax4(a,b,c,d) getmax(getmax3(a,b,c),d)
extern int n_threads;
extern int selected_num;
// #ifdef __cplusplus
// extern "C" {
// #endif
	
double cputime();
double realtime();

float **pre_process(std::vector<std::string> file_list, char *dir, int _k, int selected_number);
// #ifdef _PE
void pre_process(const char *fn, const char *fn0, int f, int w, int k, int b, int n_threads);
// #else
void pre_process(const char *fn, int f, int w, int k, int b, int n_threads);
// #endif
int cmp(const void *a, const void *b);
// compute minimizers
uint32_t sketch_suffix(const char *str,int len,int k);
uint32_t sketch_prefix(const char *str,int len,int k);
uint64_t mm_sketch_single(const char *str, int len,int k);
int mm_sketch_pos(const char *str, int len,int k);
int mm_sketch_pos2(const char *str, int len,int k);
std::array<int,2> mm_sketch_value_pos(const char *str, int len, int k);
void add_minmizer(uint64_t id);
void normalize(float **matrix, int row, int column);
//void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mmI_v *p);

// #ifdef __cplusplus
// }
// #endif

//v - ; fp is ofstream; x is the value; n is the number, if v.n == n output v.a, and a.n = 0;
#define DNA_push(v, fp, x) do {									\
		(v).a += ((x) << (2*(v).n));										\
		++ (v).n;										\
		if ((v).n == 4) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#define bit_push(v, fp, x) do {									\
		(v).a += ((x) << ((v).n));										\
		++ (v).n;										\
		if ((v).n == 8) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#define cluster_destroy(v) do { \
		(v).n = (v).m = 0;\
		free((v).ref);\
		free((v).a);\
	} while (0)


#endif