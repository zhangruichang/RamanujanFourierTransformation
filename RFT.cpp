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
#include "RFT.h"
#include<unordered_set>
#include<unordered_map>
using namespace std;

#define M_PI 3.141592653

RFT:: ~RFT()
{
    delete[] Phi;
    for(int i=0;i<=Len;i++) delete[] RS[i];
    delete[] RS;
    for(int i=0;i<4;i++) delete[] RFourier[i];
    delete[] RFourier;
}

int RFT::Gcd(int m, int n)
{
    if(!n) return m;
    return Gcd(n, m%n);
}

void RFT::OutRFTFeature(ofstream& OutFileRFTFeature, string SpeciesName, double** RFTFeature, int Seqi)
{
    OutFileRFTFeature<<SpeciesName<<"\n";
    for(int i=1;i<=Len;i++)
    {
        double PS=0;
        for(int j=0;j<4;j++) PS+=abs(RFourier[j][i])*abs(RFourier[j][i]);

        RFTFeature[Seqi][i-1]=PS;

        OutFileRFTFeature<<PS<<" ";
    }
    OutFileRFTFeature<<'\n';
}


void RFT::GetRFourier(deque<bool> *BinV)//Equation 12 in paper
{
    for(int Seqi=0;Seqi<4;Seqi++)//A T G C four binary sequences
    {
        for(int q=1;q<=Len;q++)
        {
            double c=0;
            for(int n=1;n<=Len;n++)
                c+=BinV[Seqi][n-1]*RS[q][n];
            RFourier[Seqi][q]=c/double(Len*Phi[q]);
        }
    }
}

void RFT::GetPhi()//get Phi function of Phi[1...n]
{
    for(int i=1;i<=Len;i++)
    {
        int cnt=0;
        for(int j=1;j<=i;j++)//[1, i], coprime number cnt
        {
            if(Gcd(j, i)==1) cnt++;
        }
        Phi[i]=cnt;
    }
}
void RFT::GetRamaSum()//Equation 2 in the paper
{
    for(int q=1;q<=Len;q++)
    {
        for(int n=1;n<=Len;n++)
        {
            double c=0;
            for(int p=1;p<=q;p++)
            {
                if(Gcd(p,q)!=1) continue;
                double x=2*p*n*(double)M_PI/q;
                //complex<double> tmp(cos(x), sin(x));
                c+=cos(x);//using Euler equation, e^ix=cosx+ isinx; sinx summation is zero, can be proved
            }
            RS[q][n]=c;
        }
    }
}

