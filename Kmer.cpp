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
#include "Kmer.h"
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
const int DIM = 1221;

struct DataItem {
	string id;
	string famIdx;
	double data[DIM];
};

vector<DataItem> dataSet;

void ReadDataSet(string InputFile)
{
	dataSet.clear();
	//cout<<dataSet.max_size()<<endl;
	ifstream fin(InputFile);
	if (!fin)
	{
		cout << "Cannot open the input file" << endl;
		exit(1);
	}

	int count = 0;
	string Line;
	//getline(fin, line);
	while (getline(fin, Line))
	{
		if (count%1000==0)
			cout<<"Reading reads number"<<count++<<endl;
		DataItem item;
        item.famIdx=Line;
        for(int i=0;i<DIM;i++)
            fin>>item.data[i];
		dataSet.push_back(item);
	}
	fin.close();
}


inline void writeClusterResult(const int *cls, std::string outfile)
{
	ofstream fout(outfile);
	for (int i = 0; i < dataSet.size(); i++)
		fout <<cls[i] << endl;
	fout.close();
}

//TetraFreq
int GetIndex(char S)
{
	switch(S){//A+T=G+C=4
		case 'A':
			return 0;
		case 'G':
			return 1;
		case 'C':
			return 2;
		case 'T':
			return 3;
	}
}
string GetTetra(int n, int length)
{
	char s[4]={'A','G', 'C', 'T'};
	int a[length];
	for (int i = 0; i < length; i++) {
		a[i]=n%4;
		n=n/4;
	}
	string result;
	for (int i = 0; i < length; i++) {
		result.push_back(s[a[i]]);
	}
	//delete []a;
	return result;
}
double* GetFreqs(string seq, int length)
{

	transform(seq.begin(), seq.end(), seq.begin(),::toupper);
	int dimension=(int) pow(4.0, length);
	double* freqs=new double[dimension]();

	int a[length];
	int b[length];
	bool hasN=false;
	string charset="ATCG";
	for (int i = 0; i < seq.length()-length+1; i++) {//for each k-mer
		int index1=0;//dna index i
		int index2=0;//reverse dna index i
		//a=new int[length];//k-mer index seq
		//b=new int[length];//reverse k-mer index seq


		for (int j = 0; j < length; j++) {//for each nucleotide of k-mer
			if (charset.find(seq[i+j])==-1){hasN=true; break;}//if seq has N
			a[j]=GetIndex(seq[i+j]);
			index1+=(int)pow(4.0, j)*a[j];//AAAA 0 GAAA 1 CAAA 2 TAAA 3
			// AGCT 0 0123
		}
		if (hasN){continue;}//AACT reverse is TTGA
		for (int j = 0; j < length; j++) {
			b[length-1-j]=3-GetIndex(seq[i+j]);
			index2+=(int)pow(4.0, length-1-j)*b[length-1-j];//TTTT 0 CTTT 1 GTTT 2 ATTT 3
			//TCGA 01234
		}
		freqs[index1]++;// coresponding fregs count++
		freqs[index2]++;
	}
	for (int i = 0; i < dimension; i++) {//count to frequency(<=1)
		freqs[i]=(double)freqs[i]/2/(seq.length()-length+1);//freqs /(2*k-mer number)
	}
	//delete []b;
	//delete []a;
	return freqs;
}
void PrepareSeq(string inputfile,string outfile,int length)//write k-mer feature to ***.freq files
{
		ifstream fin(inputfile);
		ofstream fout(outfile);
		fout<<"feature\t";
		int FeatureNum=(int)pow(4.0, length);
		for (int i = 0; i < FeatureNum; i++)
			fout<<GetTetra(i,length)<<"\t\n"[i==(FeatureNum-1)];
		string title, line, seq;
		int dimension=(int)pow(4.0,length);
		double* freqs=new double[dimension];
		int cnt=0;
		while(getline(fin, line))
		{
			if(line.size()>0 && line[0]=='>')
			{
				if(title=="")
				{
					title=line;
					continue;
				}
				freqs=GetFreqs(seq,length);
				fout<<title<<'\n';
				for (int i = 0; i < dimension; i++)
					fout<<freqs[i]<<" \n"[i==(dimension-1)];
				title=line;
				//seq="";
			}
			else
				seq=line;
		}
        freqs=GetFreqs(seq,length);
        fout<<title<<'\n';
        for (int i = 0; i < dimension; i++)
            fout<<freqs[i]<<" \n"[i==(dimension-1)];
		delete []freqs;
		fout.close();
		fin.close();
}
