#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>
//#include <ifstream>
#include <fstream>

using namespace std;
std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

int main(int argc, char **argv)
{
	char filename[100];
	if(argc < 2)
	{
		printf("error, file num is required");
		exit(0);
	}
	int num = atoi(argv[1]);
	ofstream file;
	file.open("size.csv");
	
	uint32_t sum;
	vector<uint32_t> file_size;
	for(int i = 1; i <= num; i++)
		file << "," << i;
	file << "\n";

	for(int i = 1; i <= num; i++)
	{
		sprintf(filename,"../pgrc/%d.pgrc",i);
		uint32_t buffer = filesize(filename);
		sum+=buffer;
		file_size.push_back(buffer);
	}


	for(int i = 1; i <= num; i++)
	{
		file << i;
		for (int j = 1; j <= i; j++)
			file << ",";
		for(int j = i+1; j <= num; j++)
		{
			sprintf(filename,"../pgrc/%d_%d.pgrc",i,j);
			file << "," << double(filesize(filename))/double((file_size[i-1]+file_size[j-1]));
		}
		file << "\n";
	}
	file.close();
	cout << "sum:" << sum << endl;
	return 0;
}