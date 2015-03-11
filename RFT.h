#ifndef RFT_H_INCLUDED
#define RFT_H_INCLUDED
#include<cstdio>
#include<iostream>
#include<string>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<queue>
#include<cstdlib>
#include<vector>
#include<set>
#include<map>
#include <complex>
#include <fstream>
#include<stack>
#include<climits>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
using namespace std;

typedef long long LL;
typedef unsigned long long ULL;


class RFT
{
public:
    RFT(int N):Len(N)
    {
        Phi=new int[Len+1];//1-based
        RS=new double*[Len+1];
        for(int i=0;i<=Len;i++) RS[i]=new double[Len+1];
        //double RS[maxn][maxn];//RS base function
        RFourier=new double*[4];
        for(int i=0;i<4;i++) RFourier[i]=new double[Len+1];
        //double RFourier[4][maxn];//ramanujan fourier signal

    }
    ~RFT();
    void GetRFourier(deque<bool> *BinV);
    void GetPhi();//get phi function of phi[1...n]
    void GetRamaSum();
    int GetLen(){return Len;}
    void OutRFTFeature(ofstream& OutFile, string SpeciesName, double** RFTFeature, int Seqi);
private:
    int Gcd(int m, int n);
    int* Phi;
    double** RS;
    double** RFourier;
    int Len;
};

//inline void writeClusterResult(const int *cls, int SampleNum, string outfile);

#endif // RFT_H_INCLUDED
