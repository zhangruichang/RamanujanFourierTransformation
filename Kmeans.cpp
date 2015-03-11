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
#include "Kmeans.h"
#include <windows.h>
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
#define epsilon 1e-3

double Kmeans:: Distance(double *x,double *y)//Euclean distance
{
    ofstream DistanceFile("Distances.txt");
    //freopen("Distances.txt","w",stdout);
	double distance=0;
	for(int i=0;i<Dim;i++)
		distance+=(LD)(x[i]-y[i])*(x[i]-y[i]);
    DistanceFile<<"distance: "<<distance<<endl;
    //Sleep(2000);
	return sqrt(distance);
}

bool Kmeans:: Converge()//wheather is converged
{
	double sum=0;
	for(int i=0;i<ClusterNum;i++)
		sum+=Distance(Center[i],NewCenter[i]);
	cout<<sum<<" "<<epsilon<<endl;
	return sum<epsilon;
}

void Kmeans:: Normalize(string RawStr)
{
    ofstream OutNormalizeRFT(RawStr+".Normalize");
    for(int i=0;i<SampleNum;i++)
    {
        double sum=0;
        for(int j=0;j<Dim;j++) sum+=Data[i][j];
        for(int j=0;j<Dim;j++) Data[i][j]/=sum,OutNormalizeRFT<<Data[i][j]<<" ";
        OutNormalizeRFT<<endl;
    }
}

Kmeans:: ~Kmeans()
{

    for(int i=0; i<SampleNum; i++) delete[] Data[i];
    for(int i=0; i<ClusterNum; i++)
    {
        delete[] ClusterRes[i];
        delete[] Center[i];
        delete[] NewCenter[i];
    }
    delete[] Data;
    delete[] ClusterRes;
    delete[] Center;
    delete[] NewCenter;
}

void Kmeans:: SetData(double** X)//random partition
{
    for(int i=0; i<SampleNum; i++) for(int j=0; j<Dim; j++) Data[i][j]=X[i][j];
    //memset(ClusterRes, 0, sizeof(ClusterRes));

    for(int i=0; i<ClusterNum; i++)
    {
        int Randi=rand()%SampleNum;
        for(int j=0; j<Dim; j++)
            NewCenter[i][j]=X[Randi][j];
    }
}

void Kmeans::Clustering()
{
    int IterativeCount=0;
    //iterative for clustering
    int MaxIterativeCnt=100;
    while(!Converge() && IterativeCount<MaxIterativeCnt)
    {
        //update w;
        for(int i=0; i<ClusterNum; i++)
            for(int j=0; j<Dim; j++)
                Center[i][j]=NewCenter[i][j];
        //compute d_ij
        for(int i=0; i<SampleNum; i++)
        {
            double mindis=Distance(Data[i],Center[0]);
            int mindis_j=0;
            for(int j=1; j<ClusterNum; j++)
            {
                double Dis=Distance(Data[i],Center[j]);
                if(mindis>Dis)
                {
                    mindis=Dis;
                    mindis_j=j;
                }
            }
            for(int j=0; j<ClusterNum; j++) ClusterRes[j][i]=false;
            ClusterRes[mindis_j][i]=true;
        }

        //modify cluster center


        for(int i=0; i<ClusterNum; i++)
        {
            memset(NewCenter[i], 0, sizeof(NewCenter[i]));
            int TotalNum=0;
            for(int j=0; j<SampleNum; j++)
            {
                if(ClusterRes[i][j])
                {
                    for(int k=0; k<Dim; k++)
                        //total_x[i][k]+=x[j][k];
                        NewCenter[i][k]+=Data[j][k];
                    TotalNum++;
                }
            }
            for(int k=0; k<Dim; k++)
                NewCenter[i][k]/=TotalNum;
        }
        //IterativeCount++;
        cout<<"Iteration Time: "<<++IterativeCount<<endl;
    }//end of iterative
    //cout<<"Iterative: "<< IterativeCount<<"\n";
        //output
}
void Kmeans:: ShowResult(string OutFileName)
{
    ofstream OutFile(OutFileName);
    for(int i=0;i<SampleNum;i++)
    {
        for(int j=0;j<ClusterNum;j++)
            if(ClusterRes[j][i]) OutFile<<j<<endl;
    }
    for(int i=0;i<ClusterNum;i++)
    {
        double InDis=0;
        for(int j=0;j<SampleNum;j++)
        {
            if(ClusterRes[i][j])
                InDis+=Distance(Data[j], NewCenter[i]);
        }
        cout<<"Cluster "<<i<<": "<<InDis<<endl;
    }
    OutFile.close();
}
