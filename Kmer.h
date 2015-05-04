#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED
#include<set>
#include<map>
#include<list>
#include<cmath>
#include<queue>
#include<stack>
#include<ctime>
#include<cstdio>
#include<string>
#include<vector>
#include<climits>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iostream>
#include<algorithm>
#include<unordered_set>
#include<unordered_map>
using namespace std;
typedef long long LL;
typedef unsigned long long ULL;
#define fo(i,a,b) for(int i=a;i<=b;i++)
#define mp make_pair
#define pb push_back
#define CLR(a,x) memset(a,x,sizeof(a))
#define Max 1000000000
#define Min 0
#define epsilon 1e-7

void ReadDataSet(string InputFile);
inline void writeClusterResult(const int *cls, string outfile);
int GetIndex(char S);
string GetTetra(int n, int length);
double* GetFreqs(string seq, int length);
void PrepareSeq(string inputfile,string outfile,int length);

#endif // KMER_H_INCLUDED
