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
#include "EvalPer.h"
#include<unordered_set>
#include<unordered_map>
using namespace std;
typedef long long LL;
typedef unsigned long long ULL;



void EvalPer(int SeqNum, int ClusterNum, string SeqFile, string ClusterFile)
{
    //SeqNum=173028, ClusterNum=2;
    //int ClusterNum=10;
    //int SeqNum=50000;
    unordered_map<string, int> SpeciesTable;
    vector<vector<int> > SpeciesToCluster(ClusterNum, vector<int>(ClusterNum, 0));

    ifstream FileSeq(SeqFile),
            FileCluster(ClusterFile);
    int SpeciesId=1, ClusterId;
    string Line;
    while(getline(FileSeq, Line))
    {
        string SpeciesName=Line;
        if(SpeciesTable.find(SpeciesName)==SpeciesTable.end())
            SpeciesTable[SpeciesName]=SpeciesId++;
        getline(FileSeq, Line);
        FileCluster>>ClusterId;
        SpeciesToCluster[SpeciesTable[SpeciesName]-1][ClusterId]++;
    }

    for(int i=0;i<ClusterNum;i++) for(int j=0;j<ClusterNum;j++)
        cout<<SpeciesToCluster[i][j]<<" \n"[j==(ClusterNum-1)];

    double Sensitivity, Precision;
    int sum=0;
    for(int i=0;i<ClusterNum;i++)
        sum+=*max_element(SpeciesToCluster[i].begin(), SpeciesToCluster[i].end());
    Sensitivity=(double)sum/SeqNum;

    sum=0;
    for(int j=0;j<ClusterNum;j++)
    {
        int maxe=INT_MIN;
        for(int i=0;i<ClusterNum;i++)
            maxe=max(maxe, SpeciesToCluster[i][j]);
        sum+=maxe;
    }
    Precision=(double)sum/SeqNum;
    cout<<"Sensitivity: "<<Sensitivity<<" Precision: "<<Precision<<" F1: "<<2*Sensitivity*Precision/(Sensitivity+Precision)<<endl;;
}
