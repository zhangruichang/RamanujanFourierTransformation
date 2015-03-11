#ifndef KMEANS_H_INCLUDED
#define KMEANS_H_INCLUDED
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
typedef long double LD;
#define fo(i,a,b) for(int i=a;i<=b;i++)
#define mp make_pair
#define pb push_back
#define CLR(a,x) memset(a,x,sizeof(a))
#define Max 1000000000
#define Min 0
#define epsilon 1e-7
class Kmeans
{
    private:
        int Dim;
        int SampleNum;
        int ClusterNum;
        double **Data;
        bool **ClusterRes;
        double** Center, **NewCenter;

        double Distance(double *x,double *y);//Euclean distance;
        //bool Converge(double** w,double** new_w,int c,int dim);//wheather is converged;
        bool Converge();
    public:

        Kmeans(int Dim_, int SampleNum_, int ClusterNum_):Dim(Dim_), SampleNum(SampleNum_), ClusterNum(ClusterNum_)
        {
            Data=new double*[SampleNum];
            for(int i=0; i<SampleNum; i++) Data[i]=new double[Dim]();

            ClusterRes=new bool*[ClusterNum];
            for(int i=0; i<ClusterNum; i++) ClusterRes[i]=new bool[SampleNum]();

            Center=new double*[ClusterNum];
            NewCenter=new double*[ClusterNum];
            for(int i=0; i<ClusterNum; i++) Center[i]=new double[Dim](), NewCenter[i]=new double[Dim]();
        }
        ~Kmeans();
        void SetData(double** X);//random partition
        void Clustering();
        void Normalize(string RawStr);
        void ShowResult(string OutFileName);
};



#endif // KMEANS_H_INCLUDED
