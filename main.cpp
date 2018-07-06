//
//  main.cpp
//  kaylailene_july6_2018
//
//  Created by Kayla Hardwick on 7/6/18.
//  Copyright © 2018 Kayla Hardwick. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <algorithm>

//Calling Namespace
using namespace std;
char FileRoot[50],FileName[50];

//class
class individual
{
public:
    //Functions
    individual();
    ~individual();
    //Variables
    int** bdmi;
    int intrinsic;
    double* z;
    double** genome;
private:
};
//Function declarations
void input(int* gSympatryPtr,int* gAllopatryPtr,int *NPtr,int *NhybridsPtr,int *numPopsPtr,int* numTraitsPtr,int* numLociPtr,double* rPtr,double* muPtr,double* mPtr, double** thetaAPtr,double*** gammaAPtr,double*** gammaPopsPtr,double* alphaMeanPtr,double*** thetaPopsPtr,double*** muMtrxPtr,int* numBDMILociPtr, double* BDMImuPtr, double* probpairwisePtr, double* intrinsiceffectPtr,double* alphaPtr,string infilestring);
void writesiminfo(ofstream &out_Out,int gSympatry,int gAllopatry,int N,int Nhybrids,int numPops,int numTraits,int numLoci,double r,double mu,double m, double* thetaA,double** gammaA,double** gammaPops,double alphaMean,double** thetaPops,double** muMtrx,int numBDMILoci, double BDMImu, double probpairwise, double intrinsiceffect,double alpha);
void initalizeArrays( double*** GmatrixPtr,int** numAfterSelectionPtr,int**numAfterMigrationPtr,int** poptrackerPtr,int numPops,int Nhybrids,int*** incomp,int** bdmiNum,int*** bdmiMutOrder,int numBDMILoci,int*** bdmitrackerPtr,ofstream &out_TEST);
void initalizePopulations(individual**ancestorPtr, individual**ancestor2Ptr,individual** migpoolPtr, individual** wrkPopPtr, individual** hybridsPtr, individual*** parentalsPtr, individual***popsPtr,
                          individual***pops2Ptr,individual*** pops3Ptr, int numPops, int N,int numTraits,int numLoci, int Nhybrids);
void initalizeSympatricPhase(individual* ancestor,individual* ancestor2, int N, int numLoci,double alphaMean, int numTraits, int numBDMILoci,ofstream &out_TEST);
double calcFitness(individual thatoneguy,int numTraits,int numLoci,double* theta,double** gamma,double intrinsiceffect);
void selection(individual* populationIn, individual* populationOut,int numIn,int* numOutPtr,int numLoci,double** gamma, double* theta, int numTraits, int numBDMILoci,double intrinsiceffect);
void CalculatePhenotypesOfaPopulation(individual* population,int numInd,int numTraits,int numLoci);

int main(int argc, const char * argv[])
{
    string basestring="/Users/hardwickk/Desktop/";
    string dirstring;
    cout<<"Enter name of main directory \n";
    cin >>dirstring;
    int numtrials;
    int trialstart;
    cout<<"Enter number of trials \n";
    cin >>numtrials;
    cout<<"Enter start of trials \n";
    cin >>trialstart;
    for(int trial=0;trial<numtrials;trial++)
    {
        //Setting up file for output.
        ofstream out_Out;
        ofstream out_Parents;
        ofstream out_Hybrids;
        ofstream out_Parentals;
        ofstream out_G;
        ofstream out_GH;
        ofstream out_GP;
        ofstream out_Test;
        ofstream out_bdmi;
        
        ofstream out_bdmi_Parents;
        ofstream out_bdmi_Hybrids;
        ofstream out_bdmi_Parentals;
        
        ofstream out_TEST;
        
        string infilestring=basestring + dirstring + "/InPut.txt";
        string trialstring;
        int trialforstring=trialstart+trial;
        ostringstream convert;
        convert<<trialforstring;
        trialstring=convert.str();
        string outstring=basestring + dirstring + "/trial_" + trialstring + "/out.txt";
        string parentsstring=basestring + dirstring + "/trial_" + trialstring + "/parents.txt";
        string hybridsstring=basestring + dirstring + "/trial_" + trialstring + "/hybrids.txt";
        string parentalsstring=basestring + dirstring + "/trial_" + trialstring + "/parentals.txt";
        string gstring=basestring + dirstring + "/trial_" + trialstring + "/gmatrix_parents.txt";
        string ghstring=basestring + dirstring + "/trial_" + trialstring + "/gmatrix_hybrids.txt";
        string gpstring=basestring + dirstring + "/trial_" + trialstring + "/gmatrix_parentals.txt";
        string teststring=basestring + dirstring + "/trial_" + trialstring + "/matepref.txt";
        string bdmistring=basestring + dirstring + "/trial_" + trialstring + "/bdmis.txt";
        
        string bdmiparentsstring=basestring + dirstring + "/trial_" + trialstring + "/bdmis_parents.txt";
        string bdmihybridsstring=basestring + dirstring + "/trial_" + trialstring + "/bdmis_hybrids.txt";
        string bdmiparentalsstring=basestring + dirstring + "/trial_" + trialstring + "/bdmis_parentals.txt";
        
        string TESTstring=basestring + dirstring + "/trial_" + trialstring + "/TEST";
        
        out_Out.open(outstring);
        out_Parents.open(parentsstring);
        out_Hybrids.open(hybridsstring);
        out_Parentals.open(parentalsstring);
        out_G.open(gstring);
        out_GH.open(ghstring);
        out_GP.open(gpstring);
        out_Test.open(teststring);
        out_bdmi.open(bdmistring);
        
        out_bdmi_Parents.open(bdmiparentsstring);
        out_bdmi_Hybrids.open(bdmihybridsstring);
        out_bdmi_Parentals.open(bdmiparentalsstring);
        
        out_TEST.open(TESTstring);
        
        //Define variables
        srand(time(NULL));
        int numPops,N,Nhybrids,numTraits,numLoci,gAllopatry,gSympatry,numBDMILoci;
        double r, mu, m, alphaMean,alpha,BDMImu,probpairwise,intrinsiceffect;
        double *thetaA;
        double **muMtrx, **thetaPops, **gammaA, **gammaPops;
        int g;
        int* numAfterSelection,*numAfterMigration,*poptracker,*bdmiNum;
        int** incomp,**bdmiMutOrder,**bdmitracker;
        double **GMtrx;
        individual* ancestor, *ancestor2;
        individual* migpool,* wrkPop,*hybrids,** parentals, ** pops, ** pops2, **pops3;
        
        //Initialize some stuff
        input(&gSympatry,&gAllopatry,&N,&Nhybrids,&numPops,&numTraits,&numLoci,&r,&mu,&m,&thetaA,&gammaA,&gammaPops,&alphaMean,&thetaPops,&muMtrx,&numBDMILoci,&BDMImu,&probpairwise,&intrinsiceffect,&alpha,infilestring);
        writesiminfo(out_Out,gSympatry,gAllopatry,N,Nhybrids,numPops,numTraits,numLoci,r,mu,m,thetaA,gammaA,gammaPops,alphaMean,thetaPops,muMtrx,numBDMILoci,BDMImu,probpairwise,intrinsiceffect,alpha);
        initalizeArrays(&GMtrx,&numAfterSelection,&numAfterMigration,&poptracker,numPops,Nhybrids,&incomp,&bdmiNum,&bdmiMutOrder,numBDMILoci,&bdmitracker,out_TEST);
        initalizePopulations(&ancestor,&ancestor2,&migpool,&wrkPop,&hybrids,&parentals,&pops,&pops2,&pops3,numPops,N,numTraits,numLoci,Nhybrids);
        initalizeSympatricPhase(ancestor,ancestor2,N,numLoci,alphaMean,numTraits,numBDMILoci,out_TEST);
        //Sympatric phase
        numAfterSelection[0]=0;
        g=0;
        cout<<"sympatry"<<endl;
        for(int gS=0;gS<gSympatry;gS++)
        {
            //cout<<"g: "<<g<<endl;
            if(gS%100==0)
            {
                cout<<"g: "<<g<<endl;
            }
            selection(ancestor,ancestor2,N,&numAfterSelection[0],numLoci,gammaA,thetaA,numTraits,numBDMILoci,intrinsiceffect);
            g++;
        }
    }
    return 0;
}

//Class Functions
individual::individual()
{
    
}
individual::~individual()
{
    
}
//Functions
void input(int* gSympatryPtr,int* gAllopatryPtr,int *NPtr,int *NhybridsPtr,int *numPopsPtr,int* numTraitsPtr,int* numLociPtr,double* rPtr,double* muPtr,double* mPtr, double** thetaAPtr,double*** gammaAPtr,double*** gammaPopsPtr,double* alphaMeanPtr,double*** thetaPopsPtr,double*** muMtrxPtr,int* numBDMILociPtr, double* BDMImuPtr, double* probpairwisePtr, double* intrinsiceffectPtr,double* alphaPtr,string infilestring)
{
    fstream inFile;
    inFile.open(infilestring);
    string myString;
    myString="Nothing";
    while(myString != "generations" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*gSympatryPtr);
    myString="Nothing";
    while(myString != "generations" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*gAllopatryPtr);
    while(myString != "individuals" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*NPtr);
    while(myString != "hybrids" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*NhybridsPtr);
    while(myString != "populations" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*numPopsPtr);
    while(myString != "traits" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*numTraitsPtr);
    while(myString != "loci" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*numLociPtr);
    while(myString != "rate" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*rPtr);
    while(myString != "ancestral" && inFile.good())
    {
        inFile>>myString;
    }
    (*gammaAPtr)=new double*[(*numTraitsPtr)];
    for(int t1=0;t1<(*numTraitsPtr);t1++)
    {
        (*gammaAPtr)[t1]=new double[(*numTraitsPtr)];
        for(int t2=0;t2<(*numTraitsPtr);t2++)
        {
            inFile>>(*gammaAPtr)[t1][t2];
        }
    }
    while(myString != "pops" && inFile.good())
    {
        inFile>>myString;
    }
    (*gammaPopsPtr)=new double*[(*numTraitsPtr)];
    for(int t1=0;t1<(*numTraitsPtr);t1++)
    {
        (*gammaPopsPtr)[t1]=new double[(*numTraitsPtr)];
        for(int t2=0;t2<(*numTraitsPtr);t2++)
        {
            inFile>>(*gammaPopsPtr)[t1][t2];
        }
    }
    while(myString != "ancestral" && inFile.good())
    {
        inFile>>myString;
    }
    (*thetaAPtr)= new double[(*numTraitsPtr)];
    for(int t=0;t<(*numTraitsPtr);t++)
    {
        inFile>>((*thetaAPtr)[t]);
    }
    while(myString != "pops" && inFile.good())
    {
        inFile>>myString;
    }
    (*thetaPopsPtr)= new double*[(*numPopsPtr)];
    for(int p=0;p<(*numPopsPtr);p++)
    {
        (*thetaPopsPtr)[p]=new double[(*numTraitsPtr)];
        for(int t=0;t<(*numTraitsPtr);t++)
        {
            inFile>>(*thetaPopsPtr)[p][t];
        }
    }
    myString="Nothing";
    while(myString != "matrix" && inFile.good())
    {
        inFile>>myString;
    }
    (*muMtrxPtr)= new double*[(*numTraitsPtr)];
    for(int t1=0;t1<(*numTraitsPtr);t1++)
    {
        (*muMtrxPtr)[t1]=new double[(*numTraitsPtr)];
        for(int t2=0;t2<(*numTraitsPtr);t2++)
        {
            inFile>>(*muMtrxPtr)[t1][t2];
        }
    }
    while(myString != "mutation" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*muPtr);
    while(myString != "rate" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*mPtr);
    while(myString != "mean" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*alphaMeanPtr);
    while(myString != "alpha" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*alphaPtr);
    while(myString != "loci" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*numBDMILociPtr);
    while(myString != "rate" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*BDMImuPtr);
    while(myString != "incompatibility" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*probpairwisePtr);
    while(myString != "effect" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*intrinsiceffectPtr);
}
void writesiminfo(ofstream &out_Out,int gSympatry,int gAllopatry,int N,int Nhybrids,int numPops,int numTraits,int numLoci,double r,double mu,double m, double* thetaA,double** gammaA,double** gammaPops,double alphaMean,double** thetaPops,double** muMtrx,int numBDMILoci, double BDMImu, double probpairwise, double intrinsiceffect,double alpha)//what else do i want to add?
{
    cout<<"gSympatry: "<<gSympatry<<endl;
    out_Out<<"gSympatry: "<<gSympatry<<endl;
    cout<<"gAllopatry: "<<gAllopatry<<endl;
    out_Out<<"gAllopatry: "<<gAllopatry<<endl;
    cout<<"N: "<<N<<endl;
    out_Out<<"N: "<<N<<endl;
    cout<<"Nhybrids: "<<Nhybrids<<endl;
    out_Out<<"Nhybrids: "<<Nhybrids<<endl;
    cout<<"numPops: "<<numPops<<endl;
    out_Out<<"numPops: "<<numPops<<endl;
    cout<<"numTraits: "<<numTraits<<endl;
    out_Out<<"numTraits: "<<numTraits<<endl;
    cout<<"numLoci: "<<numLoci<<endl;
    out_Out<<"numLoci: "<<numLoci<<endl;
    cout<<"r: "<<r<<endl;
    out_Out<<"r: "<<r<<endl;
    cout<<"gamma ancestral: "<<endl;
    out_Out<<"gamma ancestral: "<<endl;
    for(int t1=0;t1<numTraits;t1++)
    {
        for(int t2=0;t2<numTraits;t2++)
        {
            cout<<gammaA[t1][t2]<<" ";
            out_Out<<gammaA[t1][t2]<<" ";
        }
    }
    cout<<endl;
    out_Out<<endl;
    cout<<"gamma populations: "<<endl;
    out_Out<<"gamma populations: "<<endl;
    for(int t1=0;t1<numTraits;t1++)
    {
        for(int t2=0;t2<numTraits;t2++)
        {
            cout<<gammaPops[t1][t2]<<" ";
            out_Out<<gammaPops[t1][t2]<<" ";
        }
    }
    cout<<endl;
    out_Out<<endl;
    cout<<"theta ancestral"<<endl;
    out_Out<<"theta ancestral"<<endl;
    for(int t=0;t<numTraits;t++)
    {
        cout<<thetaA[t]<<" ";
        out_Out<<thetaA[t]<<" ";
    }
    cout<<endl;
    out_Out<<endl;
    cout<<"thetaPops: "<<endl;
    out_Out<<"thetaPops: "<<endl;
    for(int p=0;p<numPops;p++)
    {
        for(int t=0;t<numTraits;t++)
        {
            cout<<thetaPops[p][t]<<" ";
            out_Out<<thetaPops[p][t]<<" ";
        }
    }
    cout<<endl;
    out_Out<<endl;
    cout<<"muMtrx: "<<endl;
    out_Out<<"muMtrx: "<<endl;
    for(int t1=0;t1<numTraits;t1++)
    {
        for(int t2=0;t2<numTraits;t2++)
        {
            cout<<muMtrx[t1][t2]<<" ";
            out_Out<<muMtrx[t1][t2]<<" ";
        }
    }
    cout<<endl;
    out_Out<<endl;
    cout<<"mu: "<<mu<<endl;
    out_Out<<"mu: "<<mu<<endl;
    cout<<"m: "<<m<<endl;
    out_Out<<"m: "<<m<<endl;
    cout<<"alphaMean: "<<alphaMean<<endl;
    out_Out<<"alphaMean: "<<alphaMean<<endl;
    cout<<"alpha: "<<alpha<<endl;
    out_Out<<"alpha: "<<alpha<<endl;
    cout<<"numBDMILoci: "<<numBDMILoci<<endl;
    out_Out<<"numBDMILoci: "<<numBDMILoci<<endl;
    cout<<"BDMImu: "<<BDMImu<<endl;
    out_Out<<"BDMImu: "<<BDMImu<<endl;
    cout<<"probpairwise: "<<probpairwise<<endl;
    out_Out<<"probpairwise: "<<probpairwise<<endl;
    cout<<"intrinsiceffect: "<<intrinsiceffect<<endl;
    out_Out<<"intrinsiceffect: "<<intrinsiceffect<<endl;
}
void initalizeArrays( double*** GmatrixPtr,int** numAfterSelectionPtr,int**numAfterMigrationPtr,int** poptrackerPtr,int numPops,int Nhybrids,int*** incomp,int** bdmiNum,int*** bdmiMutOrder,int numBDMILoci,int*** bdmitrackerPtr,ofstream &out_TEST)
{
    (*GmatrixPtr)=new double*[2];
    (*GmatrixPtr)[0]=new double[2];
    (*GmatrixPtr)[1]=new double[2];
    (*numAfterSelectionPtr)= new int[numPops];
    (*numAfterMigrationPtr)= new int[numPops];
    (*poptrackerPtr)=new int[Nhybrids*2];
    (*incomp)=new int*[2];
    (*bdmiMutOrder)=new int*[2];
    for(int p=0;p<2;p++)
    {
        (*incomp)[p]=new int[numBDMILoci];
        (*bdmiMutOrder)[p]=new int[numBDMILoci];
        for(int l=0;l<numBDMILoci;l++)
        {
            (*incomp)[p][l]=777;
            (*bdmiMutOrder)[p][l]=0;
        }
    }
    out_TEST<<"incomp"<<endl;
    for(int p=0;p<2;p++)
    {
        for(int l=0;l<numBDMILoci;l++)
        {
        out_TEST<<(*incomp)[p][l]<<" ";
        }
        out_TEST<<endl;
    }
    out_TEST<<"bdmiMutOrder"<<endl;
    for(int p=0;p<2;p++)
    {
        for(int l=0;l<numBDMILoci;l++)
        {
            out_TEST<<(*bdmiMutOrder)[p][l]<<" ";
        }
        out_TEST<<endl;
    }
    (*bdmiNum)=new int[2];
    (*bdmiNum)[0]=0;
    (*bdmiNum)[1]=0;
    out_TEST<<"bdmiNum"<<endl;
    for(int p=0;p<2;p++)
    {
        out_TEST<<(*bdmiNum)[p]<<" ";
    }
    out_TEST<<endl;
    (*bdmitrackerPtr)=new int*[2];
    (*bdmitrackerPtr)[0]=new int[numBDMILoci];
    (*bdmitrackerPtr)[1]=new int[numBDMILoci];
    for(int l=0;l<numBDMILoci;l++)
    {
        (*bdmitrackerPtr)[0][l]=0;
        (*bdmitrackerPtr)[1][l]=0;//why do this differently from incomp and bdmiMutOrder?
    }
    out_TEST<<"bdmitracker"<<endl;
    for(int p=0;p<2;p++)
    {
        for(int l=0;l<numBDMILoci;l++)
        {
            out_TEST<<(*bdmitrackerPtr)[p][l]<<" ";
        }
        out_TEST<<endl;
    }
}
void initalizePopulations(individual**ancestorPtr, individual**ancestor2Ptr,individual** migpoolPtr, individual** wrkPopPtr, individual** hybridsPtr, individual*** parentalsPtr, individual***popsPtr,
                          individual***pops2Ptr,individual*** pops3Ptr, int numPops, int N,int numTraits,int numLoci, int Nhybrids)
{
    (*ancestorPtr)=new individual[N];
    (*ancestor2Ptr)=new individual[N];
    (*migpoolPtr)=new individual[N];
    (*wrkPopPtr)=new individual[Nhybrids*2];
    (*hybridsPtr)=new individual[Nhybrids];
    (*parentalsPtr)=new individual*[2];
    (*popsPtr)= new individual*[numPops];
    (*pops2Ptr)= new individual*[numPops];
    (*pops3Ptr)= new individual*[numPops];
    for(int p=0;p<numPops;p++)
    {
        (*popsPtr)[p]=new individual[N];
        (*pops2Ptr)[p]=new individual[N];
        (*pops3Ptr)[p]=new individual[N];
    }
    (*parentalsPtr)[0]=new individual[(Nhybrids/2)];
    (*parentalsPtr)[1]=new individual[(Nhybrids/2)];
}
void initalizeSympatricPhase(individual* ancestor,individual* ancestor2, int N, int numLoci,double alphaMean, int numTraits, int numBDMILoci,ofstream &out_TEST)
{
    //double alphalocus=alphaMean/(double)numLoci;// mean spread over the loci.
    for(int i=0;i<N;i++)
    {
        ancestor[i].genome=new double*[numTraits];
        ancestor2[i].genome=new double*[numTraits];
        ancestor[i].z=new double[numTraits];
        ancestor2[i].z=new double[numTraits];
        ancestor[i].bdmi=new int*[2];
        ancestor2[i].bdmi=new int*[2];
        ancestor[i].intrinsic=0;
        ancestor2[i].intrinsic=0;
        for(int pset=0;pset<2;pset++)
        {
            ancestor[i].bdmi[pset]=new int[numBDMILoci];
            ancestor2[i].bdmi[pset]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                ancestor[i].bdmi[pset][l]=0;
                ancestor2[i].bdmi[pset][l]=0;
            }
        }
        for(int t=0;t<numTraits;t++)
        {
            ancestor[i].genome[t]=new double[numLoci];
            ancestor2[i].genome[t]=new double[numLoci];
            ancestor[i].z[t]=0.0;
            ancestor2[i].z[t]=0.0;
            for(int l=0;l<numLoci;l++)
            {
                ancestor[i].genome[t][l]=0.0;
                //ancestor[i].genome[1][l]=alphalocus;
                ancestor2[i].genome[t][l]=0.0;
            }
        }
    }
    CalculatePhenotypesOfaPopulation(ancestor,N,numTraits,numLoci);
    out_TEST<<"ancestor z"<<endl;
    for(int i=0;i<N;i++)
    {
        for(int t=0;t<numTraits;t++)
        {
            out_TEST<<ancestor[i].z[t]<<" ";
        }
        out_TEST<<endl;
    }
    out_TEST<<"ancestor genome"<<endl;
    for(int i=0;i<N;i++)
    {
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                out_TEST<<ancestor[i].genome[t][l]<<" ";
            }
            out_TEST<<endl;
        }
        out_TEST<<endl;
    }
    out_TEST<<"ancestor bdmis"<<endl;
    for(int i=0;i<N;i++)
    {
        for(int pset=0;pset<2;pset++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                out_TEST<<ancestor[i].bdmi[pset][l]<<" ";
            }
            out_TEST<<endl;
        }
        out_TEST<<endl;
    }
}
void CalculatePhenotypesOfaPopulation(individual* population,int numInd,int numTraits,int numLoci)
{
    for(int i=0;i<numInd;i++)
    {
        for(int t=0;t<numTraits;t++)
        {
            population[i].z[t]=0.0;//makes sure z is initally zero
            for(int l=0;l<numLoci;l++)
            {
                population[i].z[t]+=population[i].genome[t][l];
            }
        }
    }
}
void selection(individual* populationIn, individual* populationOut,int numIn,int* numOutPtr,int numLoci,double** gamma, double* theta, int numTraits, int numBDMILoci,double intrinsiceffect)
{
    (*numOutPtr)=0;
    double w;
    int numDies=0;
    for(int i=0;i<numIn;i++)//Runs through parents
    {
        //Calculating individuals fitness
        w=calcFitness(populationIn[i],numTraits,numLoci,theta,gamma,intrinsiceffect);
        if(rand()/(double)RAND_MAX<w)//individual survives
        {
            for(int t=0;t<numTraits;t++)
            {
                for(int l=0;l<numLoci;l++)
                {
                    populationOut[(*numOutPtr)].genome[t][l]=populationIn[i].genome[t][l];
                }
            }
            for(int j=0;j<2;j++)
            {
                for(int l=0;l<numBDMILoci;l++)
                {
                    populationOut[(*numOutPtr)].bdmi[j][l]=populationIn[i].bdmi[j][l];
                }
            }
            populationOut[(*numOutPtr)].intrinsic=populationIn[i].intrinsic;
            (*numOutPtr)++;
        }
        else//individual dies
        {
            numDies++;
        }
    }
    CalculatePhenotypesOfaPopulation(populationOut,(*numOutPtr),numTraits,numLoci);
}
double calcFitness(individual thatoneguy,int numTraits,int numLoci,double* theta,double** gamma,double intrinsiceffect)
{
    double w=0.0;
    double* vm=new double[numTraits];
    double* v1m= new double[numTraits];
    double exponentm=0.0;
    for(int t1=0;t1<numTraits;t1++)
    {
        vm[t1]=thatoneguy.z[t1]-theta[t1];
    }
    for(int t1=0;t1<numTraits;t1++)
    {
        v1m[t1]=0.0;
        for(int t2=0;t2<numTraits;t2++)
        {
            //v1m[t1]+=gamma*vm[t1];
            v1m[t1]+=gamma[t1][t2]*vm[t1];
        }
    }
    //Calculating vT*v1=exponent I store the exponent in w;
    exponentm=0.0;
    for(int t1=0;t1<numTraits;t1++)
    {
        exponentm+=vm[t1]*v1m[t1];
    }
    exponentm=-0.5*exponentm;
    //Calculating fitnesses
    w=exp(exponentm);
    cout<<exponentm<<" "<<w<<endl;
    //w=w*(pow((1-intrinsiceffect),thatoneguy.intrinsic));
    delete[] v1m;
    delete[] vm;
    return w;
}

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
