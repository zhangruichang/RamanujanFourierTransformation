#include<set>
#include<map>
#include<list>
#include<cmath>
#include<queue>
#include<stack>
#include<ctime>
#include<bitset>
#include<cstdio>
#include<string>
#include<vector>
#include<climits>
#include<cstdlib>
#include<cstring>
#include<iostream>
#include<algorithm>

#include "Kmeans.h"
#include "RFT.h"
#include "Skwic.h"
#include "Kmer.h"
#include "EvalPer.h"
#include "DFT.h"


#include<unordered_set>
#include<unordered_map>
using namespace std;




int Init(int& MaxLen, int& MinLen, int& ClusterNum, int& SampleNum, string RawStr)//return maxlen
//read source meta seq files, and write processed seqs(RawStr.seq) to files, and get maxlen and sample number
{
    //freopen("out.txt","w",stdout);
    ifstream InFile(RawStr);
    ofstream OutFile(RawStr+".seq");
    string DnaSeq, LineSeq, SpeciesName;
    unordered_set<string> SpeciesHash;
    cout<<"Preprocessing Metagenomic Sequences File..."<<endl;


    //Get MaxLen, MinLen
    while(getline(InFile, LineSeq))
    {
        if(LineSeq.size()>=1 && LineSeq[0]=='>')
        {
            if(SpeciesName!="")
            {
                MaxLen=max(MaxLen, (int)DnaSeq.size());
                MinLen=min(MinLen, (int)DnaSeq.size());
                DnaSeq="";
            }
            SpeciesName=LineSeq;
        }
        else DnaSeq+=LineSeq;
    }
    if(SpeciesName!="")
    {
        MaxLen=max(MaxLen, (int)DnaSeq.size());
        MinLen=min(MinLen, (int)DnaSeq.size());
    }
    InFile.close();
    InFile.open(RawStr, ifstream::in);
    DnaSeq="", LineSeq="", SpeciesName="";
    //OutPut Processed MetaSeq to Files
    while(getline(InFile, LineSeq))
    {
        if(LineSeq.size()>=1 && LineSeq[0]=='>')
        {
            if(SpeciesName!="")
            {
                OutFile<<SpeciesName<<"\n"<<DnaSeq<<endl;
                DnaSeq="";
                SampleNum++;
            }
            istringstream Istr(LineSeq);
            string SplitWord;
            while(getline(Istr, SplitWord, '|'));
            istringstream Istr2(SplitWord);
            int Cnt=0;
            while(Cnt<2 && getline(Istr2,SplitWord,'\"')) Cnt++;
            SpeciesName=">"+SplitWord;
            SpeciesHash.insert(SplitWord);
        }
        else
            DnaSeq+=LineSeq;
    }
    OutFile<<SpeciesName<<"\n"<<DnaSeq<<endl; SampleNum++;
    ClusterNum=SpeciesHash.size();
}

void GetRFTFeature(RFT& rft, double** RFTFeature, string RawStr)
//Calculating RFT feature and Output to Files
{

    cout<<"Calculating Ramanujan Sum..."<<endl;
    //pre computer
    rft.GetRamaSum();

    cout<<"Calculating Euler totient function..."<<endl;
    rft.GetPhi();

    cout<<"Ramanujan Fourier Transformation ...."<<endl;
    ifstream InfileMetaSeq(RawStr+".seq");
    ofstream OutFileRFTFeature(RawStr+".rft");

    string LineSeq, SpeciesName="", DnaSeq;
    deque<bool> BinV[4];//(MaxLen, 0);//A T G C;
    int Seqi=0;

    //store RFT feature
    //0-based

    while(getline(InfileMetaSeq, LineSeq))
    {
        if(LineSeq.size()>0 && LineSeq[0]=='>')
        {

            if(SpeciesName=="")
            {
                SpeciesName=LineSeq;
                continue;
            }

            if(Seqi%100==0)
                cout<<"Processing Seq: "<<Seqi<<" ";

            for(int i=0;i<4;i++) BinV[i].assign(rft.GetLen(), 0);
            for(int i=0;i<DnaSeq.size();i++)
            {
                if(DnaSeq[i]=='A') BinV[0][i]=1;
                else if(DnaSeq[i]=='T') BinV[1][i]=1;
                else if(DnaSeq[i]=='G') BinV[2][i]=1;
                else BinV[3][i]=1;
            }

            rft.GetRFourier(BinV);
            rft.OutRFTFeature(OutFileRFTFeature, SpeciesName, RFTFeature, Seqi);

            Seqi++;
            SpeciesName=LineSeq;
        }
        else
        {
            DnaSeq=LineSeq;
        }
    }

    for(int i=0;i<4;i++) BinV[i].assign(rft.GetLen(), 0);
    for(int i=0;i<DnaSeq.size();i++)
    {
        if(DnaSeq[i]=='A') BinV[0][i]=1;
        else if(DnaSeq[i]=='T') BinV[1][i]=1;
        else if(DnaSeq[i]=='G') BinV[2][i]=1;
        else BinV[3][i]=1;
    }

    rft.GetRFourier(BinV);
    rft.OutRFTFeature(OutFileRFTFeature, SpeciesName, RFTFeature, Seqi);
}


void GetBinSeq(int MaxLen,string RawStr)
//Calculating DNA Bin Seq and Output to Files
{
    ifstream InfileMetaSeq(RawStr+".seq");
    ofstream OutBinA(RawStr+"BinA.txt"),
    OutBinT(RawStr+"BinT.txt"),
    OutBinG(RawStr+"BinG.txt"),
    OutBinC(RawStr+"BinC.txt");
    string LineSeq, SpeciesName="", DnaSeq;
    deque<bool> BinV[4];//(MaxLen, 0);//A T G C;
    int Seqi=0;
    while(getline(InfileMetaSeq, LineSeq))
    {
        if(LineSeq.size()>0 && LineSeq[0]=='>')
        {

            if(SpeciesName=="")
            {
                SpeciesName=LineSeq;
                continue;
            }
            if(Seqi%100==0)
                cout<<"Processing Seq: "<<Seqi<<" ";
            for(int i=0;i<4;i++) BinV[i].assign(MaxLen, 0);
            for(int i=0;i<DnaSeq.size();i++)
            {
                if(DnaSeq[i]=='A') BinV[0][i]=1;
                else if(DnaSeq[i]=='T') BinV[1][i]=1;
                else if(DnaSeq[i]=='G') BinV[2][i]=1;
                else BinV[3][i]=1;
            }
            for(int i=0;i<MaxLen;i++)
            {
                OutBinA<<BinV[0][i]<<" ";
                OutBinT<<BinV[1][i]<<" ";
                OutBinG<<BinV[2][i]<<" ";
                OutBinC<<BinV[3][i]<<" ";
            }
            OutBinA<<endl, OutBinT<<endl, OutBinG<<endl, OutBinC<<endl;
            Seqi++;
            SpeciesName=LineSeq;
        }
        else
        {
            DnaSeq=LineSeq;
        }
    }
    for(int i=0;i<4;i++) BinV[i].assign(MaxLen, 0);
    for(int i=0;i<DnaSeq.size();i++)
    {
        if(DnaSeq[i]=='A') BinV[0][i]=1;
        else if(DnaSeq[i]=='T') BinV[1][i]=1;
        else if(DnaSeq[i]=='G') BinV[2][i]=1;
        else BinV[3][i]=1;
    }
    for(int i=0;i<MaxLen;i++)
    {
        OutBinA<<BinV[0][i]<<" ";
        OutBinT<<BinV[1][i]<<" ";
        OutBinG<<BinV[2][i]<<" ";
        OutBinC<<BinV[3][i]<<" ";
    }
    OutBinA<<endl, OutBinT<<endl, OutBinG<<endl, OutBinC<<endl;
}

void ReadRFTFeature(double **X, int Dim, int SampleNum, string RawStr)
{

    int ClusterId;
    string RFTFeatureStr=RawStr+".seq";
    //PrepareSeq(MetaSeqStr, KmerFreqStr, KmerLen);
    //ReadDataSet(KmerFreqStr);
    ifstream InFile(RFTFeatureStr);
    string SpeciesName;
    //getline(InFile, SpeciesName);
    for(int i=0;i<SampleNum;i++)
    {
        getline(InFile, SpeciesName);
        //cout<<SpeciesName<<endl;
        for(int j=0;j<Dim;j++)
            InFile>>X[i][j];//,cout<<X[i][j];
        getline(InFile, SpeciesName);//cin double, behind space and endofline be processed
        //cout<<SpeciesName<<endl;
    }
}



int main()
{


    string RawFastaStr="E:/FourierTransformation/Data/transorder-5class-abundance-5k-454.f3639058.fna";
    int MaxLen=0,MinLen=INT_MAX, SampleNum=0, ClusterNum=5;
    Init(MaxLen, MinLen, ClusterNum, SampleNum, RawFastaStr);
    cout<<"MaxLen: "<<MaxLen<<endl;
    cout<<"MinLen: "<<MinLen<<endl;
    cout<<"SampleNum: "<<SampleNum<<endl;
    cout<<"ClusterNum: "<<ClusterNum<<endl;

/*  double** RFTFeature=new double*[SampleNum+1];
    for(int i=0;i<=SampleNum;i++) RFTFeature[i]=new double[MinLen+1];
    RFT rft(MinLen);
    GetRFTFeature(rft, RFTFeature, RawFastaStr);
*/
    double** RFTFeature=new double*[SampleNum+1];
    for(int i=0;i<=SampleNum;i++) RFTFeature[i]=new double[MaxLen+1];

    RFT rft(MaxLen);
    GetRFTFeature(rft, RFTFeature, RawFastaStr);
    //int SampleNum=50000, MaxLen=1227;;

    //ReadRFTFeature(RFTFeature, MinLen, SampleNum, RawFastaStr);

    Kmeans kmeans(MinLen, SampleNum, ClusterNum);
    kmeans.SetData(RFTFeature);
    //kmeans.Normalize(RawFastaStr);
    kmeans.Clustering();
    //string ClusterResultFile="E:/FourierTransformation/Data/RFT_Kmeans_result_LongReads.txt";
    kmeans.ShowResult(RawFastaStr+".KmeansRes");


/*
    Skwic skwic;
    skwic.setDistMeasure(Skwic::CITYBLOCK);
    skwic.setData(RFTFeature, MaxLen, SampleNum, 1, 0.01);


	int *cls = skwic.clustering(ClusterNum);
	cout<<"Clustering Done!"<<endl<<"Saving Clustering Result.."<<endl;
	string OutFileClusterResult="E:/FourierTransformation/Data/RFT_SKWIC_Result.txt";
	skwic.writeClusterResult(cls, SampleNum, OutFileClusterResult);*/

	for(int i=0;i<=SampleNum;i++) delete[] RFTFeature[i];
	delete[] RFTFeature;
	EvalPer(SampleNum, ClusterNum, RawFastaStr, RawFastaStr+"_KmeansRes");
	return 0;
}
