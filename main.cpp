//
//  main.cpp
//  svsims_june1
//
//  Created by Kayla Hardwick on 6/1/16.
//  Copyright (c) 2016 Kayla Hardwick. All rights reserved.
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
void initalizeArrays( double*** GmatrixPtr,int** numAfterSelectionPtr,int**numAfterMigrationPtr,int** poptrackerPtr,int numPops,int Nhybrids,int*** incomp,int** bdmiNum,int*** bdmiMutOrder,int numBDMILoci,int*** bdmitrackerPtr);
void initalizePopulations(individual**ancestorPtr, individual**ancestor2Ptr,individual** migpoolPtr, individual** wrkPopPtr, individual** hybridsPtr, individual*** parentalsPtr, individual***popsPtr,
                          individual***pops2Ptr,individual*** pops3Ptr, int numPops, int N,int numTraits,int numLoci, int Nhybrids);
void initalizeSympatricPhase(individual* ancestor,individual* ancestor2, int N, int numLoci,double alphaMean, int numTraits, int numBDMILoci);
double calcFitness(individual thatoneguy,int numTraits,int numLoci,double* theta,double** gamma,double intrinsiceffect);
void selection(individual* populationIn, individual* populationOut,int numIn,int* numOutPtr,int numLoci,double** gamma, double* theta, int numTraits, int numBDMILoci,double intrinsiceffect);
void reproduction(individual* populationIn, individual* populationOut, int numIn, int numOut, double** muMtrx,int numLoci,int numTraits,double recombination,double muRate,bool prefon,int numBDMILoci,double alpha);
void initaizeAllopatricPhase(individual*** pops1Ptr, individual*** pops2Ptr,individual***pops3Ptr,individual* ancestor, individual** migpoolPtr, int N, int numTraits, int numLoci, int numPops, int numBDMILoci);
void initializeTestPops(individual** hybridsPtr,individual** wrkPopPtr,individual*** parentalsPtr,int** poptrackerPtr,int Nhybrids,int numTraits,int numLoci,int numPops, int numBDMILoci);
void deleteindividual(individual ind, int numTraits);
void MakeHybrids(individual* parentsin0, individual* parentsin1, int numInP0, int numInP1, individual* hybridsout, int Nh, int numLoci, int numTraits, double muRate, double** muMtrx, double recombination,double* theta0,double* theta1, double** gamma, int g,int gphase, int** incomp, int numBDMILoci,double intrinsiceffect,ofstream &out_Hybrids);
void MakeParentals(individual* parentsin, int numIn, individual* parentalsout, int Np, int numLoci, int numTraits, double muRate, double** muMtrx, double recombination,double* theta0,double* theta1, double** gamma, int g,int gphase, int** incomp, int numBDMILoci,double intrinsiceffect,ofstream &out_Parentals);
void findbdmi(individual* popsIn,int numIn, int** incomp,int numBDMILoci);
void printParents(individual* popIn,int g,int numInd,int numTraits,int numLoci,double* theta,double** gamma,int gphase,int pop,double intrinsiceffect,ofstream &out_Parents,int numBDMILoci);
//void bdmiMut(int numPops,int numIn, individual** popsIn,int** bdmiMutOrder,int* bdmiNum,int** incomp,int numBDMILoci,double BDMImu,double probpairwise);
void bdmiMut(int numPops,int numIn, individual** popsIn,int** bdmiMutOrder,int* bdmiNum,int** incomp,int numBDMILoci,double BDMImu,double probpairwise,int** bdmitracker);
void migration(individual** popsIn, individual** popsOut, double m, int numBefore, int *numAfter,int numPops, int numLoci, int numTraits, individual* migpool,int numBDMILoci);
void CalculatePhenotypesOfaPopulation(individual* population,int numInd,int numTraits,int numLoci);
void CalculateGmatrix(individual* population, int numInd,int numTraits,int numLoci,double*** GmatrxPtr);
double CalculateCoVariance(double* variable1,double* variable2, int length);
double* boxMuller(int SizeOfSample, double variance,double mean);
double** choldc(double **CovarianceMtrx, int n);
double** generateMultivariateSample(double** CovarianceMatrix,double* mean,int matrixSize, int SizeOfSample/*samplesize must be even*/);
double MatePref(double alpha,double z1,double z2);
void testChamber(individual* wrkPop,individual* pop0, individual* pop1,int* poptracker,int numInP0,int numInP1,int numOut,int numTraits,int numLoci,int g,int gphase,double alpha,ofstream &out_Test);
void writedmis(ofstream &bdmioutfile,int numBDMILoci,individual* popout,int numOut,int gentowrite,int poptowrite);

int main(int argc, const char * argv[])
{
    string basestring="/Users/kaylamh/sims/";
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
        
        //ofstream out_bdmiParents;
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
        
        //string bdmiparentstring=basestring + dirstring + "/trial_" + trialstring + "/bdmisparents.txt";
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
        
        //out_bdmiParents.open(bdmiparentstring);
        
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
        initalizeArrays(&GMtrx,&numAfterSelection,&numAfterMigration,&poptracker,numPops,Nhybrids,&incomp,&bdmiNum,&bdmiMutOrder,numBDMILoci,&bdmitracker);
        initalizePopulations(&ancestor,&ancestor2,&migpool,&wrkPop,&hybrids,&parentals,&pops,&pops2,&pops3,numPops,N,numTraits,numLoci,Nhybrids);
        initalizeSympatricPhase(ancestor,ancestor2,N,numLoci,alphaMean,numTraits,numBDMILoci);
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
            reproduction(ancestor2,ancestor,numAfterSelection[0],N,muMtrx,numLoci,numTraits,r,mu,0,numBDMILoci,alpha);
            g++;
        }
        //Initialize more stuff and delete pops that are no longer needed
        initializeTestPops(&hybrids,&wrkPop,&parentals,&poptracker,Nhybrids,numTraits,numLoci,numPops,numBDMILoci);
        initaizeAllopatricPhase(&pops,&pops2,&pops3,ancestor,&migpool,N,numTraits,numLoci,numPops,numBDMILoci);
        for(int i=0;i<N;i++)
        {
            deleteindividual(ancestor[i],numTraits);
            deleteindividual(ancestor2[i],numTraits);
        }
        delete[] ancestor;
        delete[] ancestor2;
        //Allopatric Phase
        for(int p=0;p<numPops;p++)
        {
            numAfterSelection[p]=0;
            numAfterMigration[p]=N/numPops;
        }
        cout<<"allopatry"<<endl;
        for(int gA=0;gA<gAllopatry;gA++)
        {
            //cout<<"g: "<<g<<endl;
            if(gA%100==0)
            {
                cout<<"g: "<<g<<endl;
                //Write out some stats on pops periodically throughout allopatric phase
                for(int p=0;p<numPops;p++)
                {
                    CalculateGmatrix(pops[p],numAfterMigration[p],numTraits,numLoci,&GMtrx);
                    out_G<<gA<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" P"<<p<<endl;
                    printParents(pops[p],g,numAfterMigration[p],numTraits,numLoci,thetaPops[p],gammaPops,gA,p,intrinsiceffect,out_Parents,numBDMILoci);
                }
                MakeHybrids(pops[0],pops[1],numAfterMigration[0],numAfterMigration[1],hybrids,Nhybrids,numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,g,gA,incomp,numBDMILoci,intrinsiceffect,out_Hybrids);
                CalculateGmatrix(hybrids,Nhybrids,numTraits,numLoci,&GMtrx);
                out_GH<<gA<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<endl;
                MakeParentals(pops[0],numAfterMigration[0],parentals[0],(Nhybrids/2),numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,g,gA,incomp,numBDMILoci,intrinsiceffect,out_Parentals);
                MakeParentals(pops[1],numAfterMigration[1],parentals[1],(Nhybrids/2),numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,g,gA,incomp,numBDMILoci,intrinsiceffect,out_Parentals);
                CalculateGmatrix(parentals[0],(Nhybrids/2),numTraits,numLoci,&GMtrx);
                out_GP<<gA<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" p0"<<endl;
                CalculateGmatrix(parentals[1],(Nhybrids/2),numTraits,numLoci,&GMtrx);
                out_GP<<gA<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" p1"<<endl;
                testChamber(wrkPop,pops[0],pops[1],poptracker,numAfterMigration[0],numAfterMigration[1],Nhybrids,numTraits,numLoci,g,gA,alpha,out_Test);
            }
            if(gA%1000==0)
            {
                for(int p=0;p<numPops;p++)
                {
                    writedmis(out_bdmi_Parents, numBDMILoci, pops[p], numAfterMigration[p],gA,p);
                }
                writedmis(out_bdmi_Hybrids, numBDMILoci, hybrids, Nhybrids,gA,0);
                writedmis(out_bdmi_Parentals, numBDMILoci, parentals[0], (Nhybrids/2),gA,0);
                writedmis(out_bdmi_Parentals, numBDMILoci, parentals[1], (Nhybrids/2),gA,1);
            }
            for(int p=0;p<numPops;p++)
            {
                selection(pops[p],pops2[p],numAfterMigration[p],&numAfterSelection[p],numLoci,gammaPops,thetaPops[p],numTraits,numBDMILoci,intrinsiceffect);
                reproduction(pops2[p],pops3[p],numAfterSelection[p],(N/numPops),muMtrx,numLoci,numTraits,r,mu,1,numBDMILoci,alpha);
            }
            //bdmiMut(numPops,N/numPops,pops3,bdmiMutOrder,bdmiNum,incomp,numBDMILoci, BDMImu, probpairwise);
            bdmiMut(numPops,N/numPops,pops3,bdmiMutOrder,bdmiNum,incomp,numBDMILoci, BDMImu, probpairwise,bdmitracker);
            findbdmi(pops3[0],N/numPops,incomp,numBDMILoci);
            findbdmi(pops3[1],N/numPops,incomp,numBDMILoci);
            migration(pops3,pops,m,(N/numPops),numAfterMigration,numPops,numLoci,numTraits,migpool,numBDMILoci);
            g++;
        }
        //collect stats at end
        for(int p=0;p<numPops;p++)
        {
            CalculateGmatrix(pops[p],numAfterMigration[p],numTraits,numLoci,&GMtrx);
            out_G<<gAllopatry<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" P"<<p<<endl;
            printParents(pops[p],gAllopatry,numAfterMigration[p],numTraits,numLoci,thetaPops[p],gammaPops,gAllopatry,p,intrinsiceffect,out_Parents,numBDMILoci);
        }
        MakeHybrids(pops[0],pops[1],numAfterMigration[0],numAfterMigration[1],hybrids,Nhybrids,numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,gAllopatry,gAllopatry,incomp,numBDMILoci,intrinsiceffect,out_Hybrids);
        CalculateGmatrix(hybrids,Nhybrids,numTraits,numLoci,&GMtrx);
        out_GH<<gAllopatry<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<endl;
        MakeParentals(pops[0],numAfterMigration[0],parentals[0],(Nhybrids/2),numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,gAllopatry,gAllopatry,incomp,numBDMILoci,intrinsiceffect,out_Parentals);
        MakeParentals(pops[1],numAfterMigration[1],parentals[1],(Nhybrids/2),numLoci,numTraits,mu,muMtrx,r,thetaPops[0],thetaPops[1],gammaPops,gAllopatry,gAllopatry,incomp,numBDMILoci,intrinsiceffect,out_Parentals);
        CalculateGmatrix(parentals[0],(Nhybrids/2),numTraits,numLoci,&GMtrx);
        out_GP<<gAllopatry<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" p0"<<endl;
        CalculateGmatrix(parentals[1],(Nhybrids/2),numTraits,numLoci,&GMtrx);
        out_GP<<gAllopatry<<" "<<GMtrx[0][0]<<" "<<GMtrx[1][1]<<" "<<GMtrx[1][0]<<" p1"<<endl;
        testChamber(wrkPop,pops[0],pops[1],poptracker,numAfterMigration[0],numAfterMigration[1],Nhybrids,numTraits,numLoci,gAllopatry,gAllopatry,alpha,out_Test);
        out_bdmi<<"incompatibility matrix"<<endl;
        for(int p=0;p<numPops;p++)
        {
            writedmis(out_bdmi_Parents, numBDMILoci, pops[p], numAfterMigration[p],gAllopatry,p);
        }
        writedmis(out_bdmi_Hybrids, numBDMILoci, hybrids, Nhybrids,gAllopatry,0);
        writedmis(out_bdmi_Parentals, numBDMILoci, parentals[0], (Nhybrids/2),gAllopatry,0);
        writedmis(out_bdmi_Parentals, numBDMILoci, parentals[1], (Nhybrids/2),gAllopatry,1);
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                out_bdmi<<incomp[j][l]<<" ";
            }
            out_bdmi<<endl;
        }
        out_bdmi<<endl;
        
        out_bdmi<<"bdmiMutOrder"<<endl;
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                out_bdmi<<bdmiMutOrder[j][l]<<" ";
            }
            out_bdmi<<endl;
        }
        out_bdmi<<endl;
        
        out_bdmi<<"bdmiNum"<<endl;
        for(int j=0;j<2;j++)
        {
            out_bdmi<<bdmiNum[j]<<" ";
        }
        out_bdmi<<endl;
        
        /*for(int i=0;i<10;i++)
         {
         for(int j=0;j<2;j++)
         {
         for(int l=0;l<numBDMILoci;l++)
         {
         out_bdmiParents<<pops[0][i].bdmi[j][l]<<" ";
         }
         out_bdmiParents<<endl;
         }
         out_bdmiParents<<endl;
         }
         out_bdmiParents<<endl;*/
        
        //DELETED EVERY TIME (because they are modified throughout course of simulation)
        delete[] thetaA;
        for(int t=0;t<numTraits;t++)
        {
            delete[] GMtrx[t];
            delete[] muMtrx[t];
            delete[] gammaA[t];
            delete[] gammaPops[t];
        }
        delete[] GMtrx;
        delete[] muMtrx;
        delete[] gammaA;
        delete[] gammaPops;
        for(int p=0;p<numPops;p++)
        {
            delete[] incomp[p];
            delete[] bdmiMutOrder[p];
            delete[] thetaPops[p];
            delete[] bdmitracker[p];
        }
        delete[] numAfterSelection;
        delete[] numAfterMigration;
        delete[] poptracker;
        delete[] incomp;
        delete[] bdmiMutOrder;
        delete[] thetaPops;
        delete[] bdmiNum;
        delete[] bdmitracker;
        for(int p=0;p<numPops;p++)
        {
            for(int i=0;i<N;i++)
            {
                deleteindividual(pops[p][i],numTraits);
                deleteindividual(pops2[p][i],numTraits);
                deleteindividual(pops3[p][i],numTraits);
            }
            delete[] pops[p];
            delete[] pops2[p];
            delete[] pops3[p];
        }
        for(int i=0;i<N;i++)
        {
            deleteindividual(migpool[i],numTraits);
        }
        for(int i=0;i<Nhybrids;i++)
        {
            deleteindividual(hybrids[i],numTraits);
        }
        for(int i=0;i<(Nhybrids/2);i++)
        {
            deleteindividual(parentals[0][i],numTraits);
            deleteindividual(parentals[1][i],numTraits);
        }
        delete[] parentals[0];
        delete[] parentals[1];
        for(int i=0;i<(Nhybrids*2);i++)
        {
            deleteindividual(wrkPop[i],numTraits);
        }
        delete[] pops;
        delete[] pops2;
        delete[] pops3;
        delete[] migpool;
        delete[] hybrids;
        delete[] parentals;
        delete[] wrkPop;
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
        cout<<thetaA[t]<<",";
        out_Out<<thetaA[t]<<",";
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
void initalizeArrays( double*** GmatrixPtr,int** numAfterSelectionPtr,int**numAfterMigrationPtr,int** poptrackerPtr,int numPops,int Nhybrids,int*** incomp,int** bdmiNum,int*** bdmiMutOrder,int numBDMILoci,int*** bdmitrackerPtr)
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
    (*bdmiNum)=new int[2];
    (*bdmiNum)[0]=0;
    (*bdmiNum)[1]=0;
    (*bdmitrackerPtr)=new int*[2];
    (*bdmitrackerPtr)[0]=new int[numBDMILoci];
    (*bdmitrackerPtr)[1]=new int[numBDMILoci];
    for(int l=0;l<numBDMILoci;l++)
    {
        (*bdmitrackerPtr)[0][l]=0;
        (*bdmitrackerPtr)[1][l]=0;
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
void initalizeSympatricPhase(individual* ancestor,individual* ancestor2, int N, int numLoci,double alphaMean, int numTraits, int numBDMILoci)
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
    w=w*(pow((1-intrinsiceffect),thatoneguy.intrinsic));
    delete[] v1m;
    delete[] vm;
    return w;
}
void reproduction(individual* populationIn, individual* populationOut, int numIn, int numOut, double** muMtrx,int numLoci,int numTraits,double recombination,double muRate,bool prefon,int numBDMILoci,double alpha)
{
    int babycounter=0;
    int* parents= new int[2];
    int parentRead=0;
    double* meanMu= new double[numTraits];
    for(int t=0;t<numTraits;t++)
    {
        meanMu[t]=0.0;
    }
    while(babycounter<numOut)
    {
        parents[0]=rand()%numIn;
        parents[1]=rand()%numIn;
        double probmating;
        if(prefon==0)
        {
            probmating=1;
        }
        else if(prefon==1)
        {
            probmating=(0.5*MatePref(alpha,populationIn[parents[0]].z[0],populationIn[parents[1]].z[0]))+(0.5*MatePref(alpha,populationIn[parents[0]].z[1],populationIn[parents[1]].z[1]));
        }
        double prob=rand()/(double)RAND_MAX;
        if(prob<probmating)
        {
            //make a baby!
            parentRead=rand()%2;
            for(int l=0;l<numLoci;l++)
            {
                for(int t=0;t<numTraits;t++)
                {
                    populationOut[babycounter].genome[t][l]=populationIn[parents[parentRead]].genome[t][l];
                }
                double mutprob=((double)rand()/(double)RAND_MAX);
                if(mutprob<muRate)//Mutation occurs
                {
                   	double** mutationVec=generateMultivariateSample(muMtrx,meanMu,numTraits,1);
                    for(int t=0;t<numTraits;t++)
                    {
                        populationOut[babycounter].genome[t][l]+=mutationVec[0][t];
                    }
                    delete[] mutationVec[0];
                    delete[] mutationVec;
                }
                double recomprob=((double)rand()/(double)RAND_MAX);
                if(recomprob<recombination)
                {
                    parentRead=1-parentRead;
                }
            }
            parentRead=rand()%2;
            for(int pset=0;pset<2;pset++)
            {
                for(int b=0;b<numBDMILoci;b++)
                {
                    populationOut[babycounter].bdmi[pset][b]=populationIn[parents[parentRead]].bdmi[pset][b];
                    double recomprob=((double)rand()/(double)RAND_MAX);
                    if(recomprob<recombination)
                    {
                        parentRead=1-parentRead;
                    }
                }
            }
            populationOut[babycounter].intrinsic=0;
            babycounter++;
        }
    }
    CalculatePhenotypesOfaPopulation(populationOut,numOut,numTraits,numLoci);
    delete[] parents;
    delete[] meanMu;
}
//void bdmiMut(int numPops,int numIn, individual** popsIn,int** bdmiMutOrder,int* bdmiNum,int** incomp,int numBDMILoci,double BDMImu,double probpairwise)
//{
//    //for each individual, try for a new mutation if there's room in the table (insert and inc counters where appropriate
//    for(int p=0;p<numPops;p++)
//    {
//        for(int i=0;i<numIn;i++)
//        {
//            for(int l=0;l<numBDMILoci;l++)
//            {
//                if(popsIn[p][i].bdmi[p][l]==0)
//                {
//                    double rmut=(double)rand()/(double)RAND_MAX;
//                    if(rmut<BDMImu)
//                    {
//                        //new mutation
//                        popsIn[p][i].bdmi[p][l]=1;
//                        //check if mutation has occured before
//                        if(incomp[p][l]==777)
//                        {
//                            //hasn't occured, generate incompatibility
//                            double probbdmi=bdmiNum[1-p]*probpairwise;
//                            double rnum=(double)rand()/(double)RAND_MAX;
//                            if(rnum<probbdmi)
//                            {
//                                //incomp, fill with position
//                                int spot=rand()%bdmiNum[1-p];
//                                int valspot=bdmiMutOrder[1-p][spot];
//                                incomp[p][l]=valspot;
//                            }
//                            else
//                            {
//                                //no incomp, fill with 555
//                                incomp[p][l]=555;
//                            }
//
//                            bdmiMutOrder[p][bdmiNum[p]]=l;
//                            bdmiNum[p]++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

void bdmiMut(int numPops,int numIn, individual** popsIn,int** bdmiMutOrder,int* bdmiNum,int** incomp,int numBDMILoci,double BDMImu,double probpairwise,int** bdmitracker)
{
    int* bdmitrackercounter=new int[2];
    bdmitrackercounter[0]=0;
    bdmitrackercounter[1]=0;
    int* bdmimutusedcount=new int[2];
    bdmimutusedcount[0]=0;
    bdmimutusedcount[1]=0;
    int whichpop;
    //cout<<"bdmiMut"<<endl;
    //for each individual, try for a new mutation if there's room in the table (insert and inc counters where appropriate
    for(int p=0;p<numPops;p++)
    {
        for(int i=0;i<numIn;i++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                if(popsIn[p][i].bdmi[p][l]==0)
                {
                    double rmut=(double)rand()/(double)RAND_MAX;
                    if(rmut<BDMImu)
                    {
                        //new mutation
                        popsIn[p][i].bdmi[p][l]=1;
                        //check if mutation has occured before
                        if(incomp[p][l]==777)
                        {
                            //record it
                            //cout<<"NEW"<<endl;
                            //cout<<"p "<<p<<" l "<<l<<endl;
                            bdmitracker[p][bdmitrackercounter[p]]=l;
                            incomp[p][l]=555;
                            bdmitrackercounter[p]++;
                            
                        }
                    }
                }
            }
        }
    }
    //cout<<"incomp before"<<endl;
    for(int p=0;p<2;p++)
    {
        for(int l=0;l<numBDMILoci;l++)
        {
            //cout<<incomp[p][l]<<" ";
        }
        //cout<<endl;
    }
    //cout<<endl;
    
    //cout<<"bdmitrackercounter "<<bdmitrackercounter[0]<<" "<<bdmitrackercounter[1]<<endl;
    //cout<<"bdmitracker"<<endl;
    /*for(int p=0;p<2;p++)
     {
     for(int l=0;l<numBDMILoci;l++)
     {
     cout<<bdmitracker[p][l]<<" ";
     }
     cout<<endl;
     }*/
    //for the total number of new mutations
    for(int b=0;b<(bdmitrackercounter[0]+bdmitrackercounter[1]);b++)
    {
        //cout<<"counting through new alleles "<<b<<endl;
        //choose which pop has mutation
        whichpop=rand()%2;
        //cout<<"whichpop "<<whichpop<<endl;
        //cout<<"counters "<<bdmimutusedcount[whichpop]<<" "<<bdmitrackercounter[whichpop]<<endl;
        if(bdmimutusedcount[whichpop]==bdmitrackercounter[whichpop])
        {
            whichpop=1-whichpop;
        }
        //cout<<"whichpop "<<whichpop<<endl;
        double probbdmi=bdmiNum[1-whichpop]*probpairwise;
        //cout<<"bdminum "<<bdmiNum[1-whichpop]<<" probbdmi "<<probbdmi<<endl;
        if(((double)rand()/(double)RAND_MAX)<probbdmi)
        {
            //generate incompatibility
            //cout<<"INCOMPATIBLE"<<endl;
            int spot=rand()%bdmiNum[1-whichpop];
            int valspot=bdmiMutOrder[1-whichpop][spot];
            //cout<<"spot "<<spot<<" valspot "<<valspot<<endl;
            incomp[whichpop][bdmitracker[whichpop][bdmimutusedcount[whichpop]]]=valspot;
        }
        bdmiMutOrder[whichpop][bdmiNum[whichpop]]=bdmitracker[whichpop][bdmimutusedcount[whichpop]];
        bdmiNum[whichpop]++;
        bdmimutusedcount[whichpop]++;
    }
    /*cout<<"bdmiNum"<<endl;
     cout<<bdmiNum[0]<<" "<<bdmiNum[1]<<endl;
     cout<<"bdmitrackercounter"<<endl;
     cout<<bdmitrackercounter[0]<<" "<<bdmitrackercounter[1]<<endl;
     cout<<"bdmimutusedcount"<<endl;
     cout<<bdmimutusedcount[0]<<" "<<bdmimutusedcount[1]<<endl;
     cout<<"bdmitracker"<<endl;
     for(int p=0;p<2;p++)
     {
     for(int l=0;l<numBDMILoci;l++)
     {
     cout<<bdmitracker[p][l]<<" ";
     }
     cout<<endl;
     }
     cout<<endl;
     cout<<"incomp"<<endl;
     for(int p=0;p<2;p++)
     {
     for(int l=0;l<numBDMILoci;l++)
     {
     cout<<incomp[p][l]<<" ";
     }
     cout<<endl;
     }
     cout<<endl;
     cout<<"bdmiMutOrder"<<endl;
     for(int p=0;p<2;p++)
     {
     for(int l=0;l<numBDMILoci;l++)
     {
     cout<<bdmiMutOrder[p][l]<<" ";
     }
     cout<<endl;
     }
     cout<<endl;*/
    delete[] bdmitrackercounter;
    delete[] bdmimutusedcount;
}
void findbdmi(individual* popsIn,int numIn, int** incomp,int numBDMILoci)
{
    //for each individual, check for incompat in entire genome
    for(int i=0;i<numIn;i++)
    {
        popsIn[i].intrinsic=0;
        for(int gset=0;gset<2;gset++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                if(popsIn[i].bdmi[gset][l]==1)
                {
                    //check for incomp
                    //look up position in incomp
                    if(incomp[gset][l]!=555)
                    {
                        if(popsIn[i].bdmi[1-gset][incomp[gset][l]]==1)
                        {
                            //incompatibility!
                            popsIn[i].intrinsic++;
                        }
                    }
                }
            }
        }
    }
}
void initializeTestPops(individual** hybridsPtr,individual** wrkPopPtr,individual*** parentalsPtr,int** poptrackerPtr,int Nhybrids,int numTraits,int numLoci,int numPops, int numBDMILoci)
{
    for(int i=0;i<Nhybrids;i++)
    {
        (*hybridsPtr)[i].genome=new double*[numTraits];
        (*hybridsPtr)[i].z=new double[numTraits];
        (*hybridsPtr)[i].bdmi=new int*[2];
        (*hybridsPtr)[i].intrinsic=0;
        for(int t=0;t<numTraits;t++)
        {
            (*hybridsPtr)[i].genome[t]=new double[numLoci];
            (*hybridsPtr)[i].z[t]=0.0;
            for(int l=0;l<numLoci;l++)
            {
                (*hybridsPtr)[i].genome[t][l]=0.0;
            }
        }
        for(int j=0;j<2;j++)
        {
            (*hybridsPtr)[i].bdmi[j]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                (*hybridsPtr)[i].bdmi[j][l]=0;
            }
        }
    }
    for(int i=0;i<(Nhybrids*2);i++)
    {
        (*wrkPopPtr)[i].genome=new double*[numTraits];
        (*wrkPopPtr)[i].z=new double[numTraits];
        (*wrkPopPtr)[i].bdmi=new int*[2];
        (*wrkPopPtr)[i].intrinsic=0;
        (*poptrackerPtr)[i]=0;
        for(int t=0;t<numTraits;t++)
        {
            (*wrkPopPtr)[i].genome[t]=new double[numLoci];
            (*wrkPopPtr)[i].z[t]=0.0;
            for(int l=0;l<numLoci;l++)
            {
                (*wrkPopPtr)[i].genome[t][l]=0.0;
            }
        }
        for(int j=0;j<2;j++)
        {
            (*wrkPopPtr)[i].bdmi[j]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                (*wrkPopPtr)[i].bdmi[j][l]=0;
            }
        }
    }
    for(int i=0;i<(Nhybrids/2);i++)
    {
        (*parentalsPtr)[0][i].genome=new double*[numTraits];
        (*parentalsPtr)[0][i].z=new double[numTraits];
        (*parentalsPtr)[1][i].genome=new double*[numTraits];
        (*parentalsPtr)[1][i].z=new double[numTraits];
        (*parentalsPtr)[0][i].bdmi=new int*[2];
        (*parentalsPtr)[1][i].bdmi=new int*[2];
        (*parentalsPtr)[0][i].intrinsic=0;
        (*parentalsPtr)[1][i].intrinsic=0;
        for(int t=0;t<numTraits;t++)
        {
            (*parentalsPtr)[0][i].z[t]=0.0;
            (*parentalsPtr)[0][i].genome[t]=new double[numLoci];
            (*parentalsPtr)[1][i].z[t]=0.0;
            (*parentalsPtr)[1][i].genome[t]=new double[numLoci];
            for(int l=0;l<numLoci;l++)
            {
                (*parentalsPtr)[0][i].genome[t][l]=0.0;
                (*parentalsPtr)[1][i].genome[t][l]=0.0;
            }
        }
        for(int j=0;j<2;j++)
        {
            (*parentalsPtr)[0][i].bdmi[j]=new int[numBDMILoci];
            (*parentalsPtr)[1][i].bdmi[j]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                (*parentalsPtr)[0][i].bdmi[j][l]=0;
                (*parentalsPtr)[1][i].bdmi[j][l]=0;
            }
        }
    }
}
void initaizeAllopatricPhase(individual*** pops1Ptr, individual*** pops2Ptr,individual***pops3Ptr,individual* ancestor, individual** migpoolPtr, int N, int numTraits, int numLoci, int numPops, int numBDMILoci)
{
    individual wrkInd;
    //Initalizing populations
    for(int p=0;p<numPops;p++)
    {
        for(int i=0;i<N;i++)
        {
            (*pops1Ptr)[p][i].genome=new double*[numTraits];
            (*pops2Ptr)[p][i].genome=new double*[numTraits];
            (*pops3Ptr)[p][i].genome=new double*[numTraits];
            (*pops1Ptr)[p][i].z=new double[numTraits];
            (*pops2Ptr)[p][i].z=new double[numTraits];
            (*pops3Ptr)[p][i].z=new double[numTraits];
            (*pops1Ptr)[p][i].bdmi=new int*[2];
            (*pops2Ptr)[p][i].bdmi=new int*[2];
            (*pops3Ptr)[p][i].bdmi=new int*[2];
            (*pops1Ptr)[p][i].intrinsic=0;
            (*pops2Ptr)[p][i].intrinsic=0;
            (*pops3Ptr)[p][i].intrinsic=0;
            for(int t=0;t<numTraits;t++)
            {
                (*pops1Ptr)[p][i].genome[t]=new double[numLoci];
                (*pops2Ptr)[p][i].genome[t]=new double[numLoci];
                (*pops3Ptr)[p][i].genome[t]=new double[numLoci];
                (*pops2Ptr)[p][i].z[t]=0.0;
                (*pops3Ptr)[p][i].z[t]=0.0;
                for (int l=0;l<numLoci;l++)
                {
                    (*pops2Ptr)[p][i].genome[t][l]=0.0;
                    (*pops3Ptr)[p][i].genome[t][l]=0.0;
                }
            }
            for(int j=0;j<2;j++)
            {
                (*pops1Ptr)[p][i].bdmi[j]=new int[numBDMILoci];
                (*pops2Ptr)[p][i].bdmi[j]=new int[numBDMILoci];
                (*pops3Ptr)[p][i].bdmi[j]=new int[numBDMILoci];
                for(int l=0;l<numBDMILoci;l++)
                {
                    (*pops1Ptr)[p][i].bdmi[j][l]=0;
                    (*pops2Ptr)[p][i].bdmi[j][l]=0;
                    (*pops3Ptr)[p][i].bdmi[j][l]=0;
                }
            }
        }
    }
    for(int i=0;i<N;i++)
    {
        (*migpoolPtr)[i].genome= new double*[numTraits];
        (*migpoolPtr)[i].z= new double[numTraits];
        (*migpoolPtr)[i].bdmi=new int*[2];
        (*migpoolPtr)[i].intrinsic=0;
        for(int t=0;t<numTraits;t++)
        {
            (*migpoolPtr)[i].genome[t]= new double[numLoci];
            (*migpoolPtr)[i].z[t]=0.0;
            for (int l=0;l<numLoci;l++)
            {
                (*migpoolPtr)[i].genome[t][l]=0.0;
            }
        }
        for(int j=0;j<2;j++)
        {
            (*migpoolPtr)[i].bdmi[j]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                (*migpoolPtr)[i].bdmi[j][l]=0;
            }
        }
    }
    for(int p=0;p<numPops;p++)
    {
        for(int i=0;i<N/numPops;i++)
        {
            wrkInd=ancestor[rand()%N];
            for(int t=0;t<numTraits;t++)
            {
                for(int l=0;l<numLoci;l++)
                {
                    (*pops1Ptr)[p][i].genome[t][l]=wrkInd.genome[t][l];
                }
            }
        }
        CalculatePhenotypesOfaPopulation((*pops1Ptr)[p],N/numPops,numTraits,numLoci);
    }
}
void deleteindividual(individual ind, int numTraits)
{
    for(int t=0;t<numTraits;t++)
    {
        delete[] ind.genome[t];
    }
    for(int gset=0;gset<2;gset++)
    {
        delete[] ind.bdmi[gset];
    }
    delete[] ind.genome;
    delete[] ind.z;
    delete[] ind.bdmi;
}
void printParents(individual* popIn,int g,int numInd,int numTraits,int numLoci,double* theta,double** gamma,int gphase,int pop,double intrinsiceffect,ofstream &out_Parents,int numBDMILoci)
{
    for(int i=0;i<numInd;i++)
    {
        double w=calcFitness(popIn[i],numTraits,numLoci,theta,gamma,intrinsiceffect);
        //write totalw, thatoneguy phenotype
        out_Parents<<gphase<<" "<<popIn[i].intrinsic<<" ";
        for(int t=0;t<numTraits;t++)
        {
            out_Parents<<popIn[i].z[t]<<" ";
        }
        out_Parents<<w<<" 0 ";
        //print genome
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                out_Parents<<popIn[i].genome[t][l]<<" ";
            }
        }
        out_Parents<<"p"<<pop<<endl;
    }
}
void MakeHybrids(individual* parentsin0, individual* parentsin1, int numInP0, int numInP1, individual* hybridsout, int Nh, int numLoci, int numTraits, double muRate, double** muMtrx, double recombination,double* theta0,double* theta1, double** gamma, int g,int gphase, int** incomp, int numBDMILoci,double intrinsiceffect,ofstream &out_Hybrids)
{
    individual* hybridparents= new individual[2];
    for(int hp=0;hp<2;hp++)
    {
        hybridparents[hp].z=new double[numTraits];
        hybridparents[hp].genome=new double*[numTraits];
        for(int t=0;t<numTraits;t++)
        {
            hybridparents[hp].z[t]=0.0;
            hybridparents[hp].genome[t]=new double[numLoci];
            for(int l=0;l<numLoci;l++)
            {
                hybridparents[hp].genome[t][l]=0.0;
            }
        }
        hybridparents[hp].bdmi=new int*[2];
        for(int gset=0;gset<2;gset++)
        {
            hybridparents[hp].bdmi[gset]=new int[numBDMILoci];
            for(int l=0;l<numBDMILoci;l++)
            {
                hybridparents[hp].bdmi[gset][l]=0;
            }
        }
        hybridparents[hp].intrinsic=0;
    }
    int hybridparentRead=rand()%2;
    double* meanMu= new double[numTraits];
    for(int t=0;t<numTraits;t++)
    {
        meanMu[t]=0.0;
    }
    for(int i=0;i<Nh;i++)//Cycles through offspring
    {
        int hybridp0=rand()%numInP0;
        int hybridp1=rand()%numInP1;
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                hybridparents[0].genome[t][l]=parentsin0[hybridp0].genome[t][l];
            }
        }
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                hybridparents[1].genome[t][l]=parentsin1[hybridp1].genome[t][l];
            }
        }
        for(int gset=0;gset<2;gset++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                hybridparents[0].bdmi[gset][l]=parentsin0[hybridp0].bdmi[gset][l];
            }
        }
        for(int gset=0;gset<2;gset++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                hybridparents[1].bdmi[gset][l]=parentsin1[hybridp1].bdmi[gset][l];
            }
        }
        for(int l=0;l<numLoci;l++)
        {
            for(int t=0;t<numTraits;t++)
            {
                hybridsout[i].genome[t][l]=hybridparents[hybridparentRead].genome[t][l];
            }
            double hybridmutprob=((double)rand()/(double)RAND_MAX);
            if(hybridmutprob<muRate)//Mutation occurs
            {
                double** hybridmutationVec=generateMultivariateSample(muMtrx,meanMu,numTraits,1);
                for(int t=0;t<numTraits;t++)
                {
                    hybridsout[i].genome[t][l]+=hybridmutationVec[0][t];
                }
                delete[] hybridmutationVec[0];
                delete[] hybridmutationVec;
            }
            double hybridrecomprob=((double)rand()/(double)RAND_MAX);
            if(hybridrecomprob<recombination)
            {
                hybridparentRead=1-hybridparentRead;
            }
        }
        hybridparentRead=rand()%2;
        for(int pset=0;pset<2;pset++)
        {
            for(int b=0;b<numBDMILoci;b++)
            {
                hybridsout[i].bdmi[pset][b]=hybridparents[hybridparentRead].bdmi[pset][b];
                double hybridrecomprob=((double)rand()/(double)RAND_MAX);
                if(hybridrecomprob<recombination)
                {
                    hybridparentRead=1-hybridparentRead;
                }
            }
        }
    }
    CalculatePhenotypesOfaPopulation(hybridsout,Nh,numTraits,numLoci);
    findbdmi(hybridsout,Nh,incomp,numBDMILoci);
    //calcfitness based on both peaks
    for(int i=0;i<Nh;i++)//Cycles through offspring
    {
        double w0=calcFitness(hybridsout[i],numTraits,numLoci,theta0,gamma,intrinsiceffect);
        double w1=calcFitness(hybridsout[i],numTraits,numLoci,theta1,gamma,intrinsiceffect);
        //double totalw=(0.5*w0)+(0.5*w1);
        //write totalw, thatoneguy phenotype
        out_Hybrids<<gphase<<" "<<hybridsout[i].intrinsic<<" ";
        for(int t=0;t<numTraits;t++)
        {
            out_Hybrids<<hybridsout[i].z[t]<<" ";
        }
        out_Hybrids<<w0<<" "<<w1<<" ";
        //print genome
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                out_Hybrids<<hybridsout[i].genome[t][l]<<" ";
            }
        }
        out_Hybrids<<endl;
    }
    for(int pparent=0;pparent<2;pparent++)
    {
        deleteindividual(hybridparents[pparent],numTraits);
    }
    delete[] hybridparents;
    delete[] meanMu;
}
void MakeParentals(individual* parentsin, int numIn, individual* parentalsout, int Np, int numLoci, int numTraits, double muRate, double** muMtrx, double recombination,double* theta0,double* theta1, double** gamma, int g,int gphase, int** incomp, int numBDMILoci,double intrinsiceffect,ofstream &out_Parentals)
{
    int* parentalparents= new int[2];
    //individual* parentalparents= new individual[2];
    int parentalparentRead=rand()%2;
    double* meanMu= new double[numTraits];
    for(int t=0;t<numTraits;t++)
    {
        meanMu[t]=0.0;
    }
    for(int i=0;i<Np;i++)//Cycles through offspring
    {
        //int parentalp0=rand()%numIn;
        //int parentalp1=rand()%numIn;
        //parentalparents[0]=parentsin[parentalp0];
        //parentalparents[1]=parentsin[parentalp1];
        parentalparents[0]=rand()%numIn;
        parentalparents[1]=rand()%numIn;
        for(int l=0;l<numLoci;l++)
        {
            for(int t=0;t<numTraits;t++)
            {
                //parentalsout[i].genome[t][l]=parentalparents[parentalparentRead].genome[t][l];
                parentalsout[i].genome[t][l]=parentsin[parentalparents[parentalparentRead]].genome[t][l];
            }
            double parentalmutprob=((double)rand()/(double)RAND_MAX);
            if(parentalmutprob<muRate)//Mutation occurs
            {
                double** parentalmutationVec=generateMultivariateSample(muMtrx,meanMu,numTraits,1);
                for(int t=0;t<numTraits;t++)
                {
                    parentalsout[i].genome[t][l]+=parentalmutationVec[0][t];
                }
                delete[] parentalmutationVec[0];
                delete[] parentalmutationVec;
            }
            double parentalrecomprob=((double)rand()/(double)RAND_MAX);
            if(parentalrecomprob<recombination)
            {
                parentalparentRead=1-parentalparentRead;
            }
        }
        parentalparentRead=rand()%2;
        for(int pset=0;pset<2;pset++)
        {
            for(int b=0;b<numBDMILoci;b++)
            {
                //parentalsout[i].bdmi[pset][b]=parentalparents[parentalparentRead].bdmi[pset][b];
                parentalsout[i].bdmi[pset][b]=parentsin[parentalparents[parentalparentRead]].bdmi[pset][b];
                double parentalrecomprob=((double)rand()/(double)RAND_MAX);
                if(parentalrecomprob<recombination)
                {
                    parentalparentRead=1-parentalparentRead;
                }
            }
        }
    }
    CalculatePhenotypesOfaPopulation(parentalsout,Np,numTraits,numLoci);
    findbdmi(parentalsout,Np,incomp,numBDMILoci);
    //calcfitness based on both peaks
    for(int i=0;i<Np;i++)//Cycles through offspring
    {
        double w0=calcFitness(parentalsout[i],numTraits,numLoci,theta0,gamma,intrinsiceffect);
        double w1=calcFitness(parentalsout[i],numTraits,numLoci,theta1,gamma,intrinsiceffect);
        //double totalw=(0.5*w0)+(0.5*w1);
        //write totalw, thatoneguy phenotype
        out_Parentals<<gphase<<" "<<parentalsout[i].intrinsic<<" ";
        for(int t=0;t<numTraits;t++)
        {
            out_Parentals<<parentalsout[i].z[t]<<" ";
        }
        out_Parentals<<w0<<" "<<w1<<" ";
        //print genome
        for(int t=0;t<numTraits;t++)
        {
            for(int l=0;l<numLoci;l++)
            {
                out_Parentals<<parentalsout[i].genome[t][l]<<" ";
            }
        }
        out_Parentals<<endl;
    }
    delete[] parentalparents;
    delete[] meanMu;
}
void migration(individual** popsIn, individual** popsOut, double m, int numBefore, int *numAfter,int numPops, int numLoci, int numTraits, individual* migpool,int numBDMILoci)
{
    int numMig=0;
    int pop=0;
    for(int p=0;p<numPops;p++)
    {
        numAfter[p]=0;
        for(int i=0; i<numBefore; i++)
        {
            if(rand()/(double)RAND_MAX>m)//individual stays
            {
                for(int l=0;l<numLoci;l++)
                {
                    for(int t=0;t<numTraits;t++)
                    {
                        popsOut[p][numAfter[p]].genome[t][l]=popsIn[p][i].genome[t][l];
                    }
                }
                for(int j=0;j<2;j++)
                {
                    for(int l=0;l<numBDMILoci;l++)
                    {
                        popsOut[p][numAfter[p]].bdmi[j][l]=popsIn[p][i].bdmi[j][l];
                    }
                }
                popsOut[p][numAfter[p]].intrinsic=popsIn[p][i].intrinsic;
                numAfter[p]++;
            }
            else
            {
                for(int l=0;l<numLoci;l++)
                {
                    for(int t=0;t<numTraits;t++)
                    {
                        migpool[numMig].genome[t][l]=popsIn[p][i].genome[t][l];
                    }
                }
                for(int j=0;j<2;j++)
                {
                    for(int l=0;l<numBDMILoci;l++)
                    {
                        migpool[numMig].bdmi[j][l]=popsIn[p][i].bdmi[j][l];
                    }
                }
                migpool[numMig].intrinsic=popsIn[p][i].intrinsic;
                numMig++;
            }
        }
    }
    for(int i=0;i<numMig;i++)
    {
        pop=rand()%numPops;
        for(int l=0;l<numLoci;l++)
        {
            for(int t=0;t<numTraits;t++)
            {
                popsOut[pop][numAfter[pop]].genome[t][l]=migpool[i].genome[t][l];
            }
        }
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<numBDMILoci;l++)
            {
                popsOut[pop][numAfter[pop]].bdmi[j][l]=migpool[i].bdmi[j][l];
            }
        }
        popsOut[pop][numAfter[pop]].intrinsic=migpool[i].intrinsic;
        numAfter[pop]++;
    }
    for(int p=0;p<numPops;p++)
    {
        CalculatePhenotypesOfaPopulation(popsOut[p], numAfter[p], numTraits, numLoci);
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
void CalculateGmatrix(individual* population, int numInd,int numTraits,int numLoci,double*** GmatrxPtr)
{
    double** traitList= new double*[numTraits];
    for(int t=0;t<numTraits;t++)
    {
        traitList[t]= new double[numInd];
        for(int i=0;i<numInd;i++)
        {
            traitList[t][i]=population[i].z[t];
        }
    }
    for(int t1=0;t1<numTraits;t1++)
    {
        for(int t2=0;t2<numTraits;t2++)
        {
            (*GmatrxPtr)[t1][t2]=CalculateCoVariance(traitList[t1],traitList[t2],numInd);
            //out_Out<<(*GmatrxPtr)[t1][t2]<<",";
        }
        //out_Out<<endl;
    }
    //out_Out<<endl;
    for(int t=0;t<numTraits;t++)
    {
        delete[] traitList[t];
    }
    delete[] traitList;
}
double CalculateCoVariance(double* variable1,double* variable2, int length)
{
    double Covariance=0.0;
    double average1=0.0;
    double average2=0.0;
    for(int i=0;i<length;i++)
    {
        average1+=variable1[i];
        average2+=variable2[i];
    }
    average1=average1/length;
    average2=average2/length;
    for(int i=0;i<length;i++)
    {
        Covariance+=(variable1[i]-average1)*(variable2[i]-average2);
    }
    Covariance=Covariance/length;
    return Covariance;
}
double* boxMuller(int SizeOfSample, double variance,double mean)
{
    if(SizeOfSample%2==1)SizeOfSample++;
    double* uniformVec=new double[SizeOfSample];
    double* normalVec=new double[SizeOfSample];
    for(int i=0; i<SizeOfSample; i++)
    {
        uniformVec[i]=((double)rand())/RAND_MAX;
    }
    for(int i=0; i<(SizeOfSample/2); i++)
    {
        normalVec[2*i]=sqrt(-2*log(uniformVec[2*i]))*cos(2*3.14159265*uniformVec[2*i+1])*((double)sqrt(variance))+mean;
        normalVec[2*i+1]=sqrt(-2*log(uniformVec[2*i]))*sin(2*3.14159265*uniformVec[2*i+1])*((double)sqrt(variance))+mean;
        //The *((double)sqrt(variance))+mean transforms the variance from 1 to .25=(1/2)^2 and makes the mean .5
    }
    
    //Prints out values of normal distribution
    /*out_Out<<"NormalVec is , "<<endl;
     for(int i=0; i<SizeOfSample; i++)
     {
     out_Out<<normalVec[i]<<endl;
     }
     out_Out<<endl;*/
    delete[] uniformVec;
    return normalVec;
}
double** choldc(double **CovarianceMtrx, int n)
{
    double**a= new double*[n];
    for(int i=0;i<n;i++)
    {
        a[i]=new double[n];
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            a[i][j]=CovarianceMtrx[i][j];
        }
    }
    double* p=new double[n]; //alocatedDub+=n;
    int i, j, k;
    double sum;
    for(i=0;i<n;i++)
    {
        for(j=i;j<n;j++)
        {
            sum=a[i][j];
            for(k=i-1;k>=0;k--)
            {
                sum=sum-a[i][k]*a[j][k];
            }
            if(i==j)
            {
                if(sum<=0.0)
                {
                    cout<<"Epic Fail"<<endl;
                }
                else
                {
                    p[i]=sqrt(sum);
                }
            }
            a[j][i]=sum/p[i];
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(j>i)
            {
                a[i][j]=0.0;
            }
        }
    }
    //Printing the matrix
    /*cout<<"Decomposed Mtrx: "<<endl;
     for(int i=0;i<n; i++)
     {
     for(int j=0;j<n; j++)
     {
     cout<<a[i][j]<<",";
     }
     cout<<endl;
     }*/
    delete[] p;// alocatedDub-=n;
    return a;
}
double** generateMultivariateSample(double** CovarianceMatrix,double* mean,int matrixSize, int SizeOfSample/*samplesize must be even*/)
{
    //double* NormalSample=new double[SizeOfSample];
    double** choldcMtrx= choldc(CovarianceMatrix,matrixSize);
    double** MultiNormalSample=new double*[SizeOfSample];
    for(int s=0; s<SizeOfSample;s++)
    {
        MultiNormalSample[s]=new double[matrixSize];
    }
    //Matrix by
    for(int s=0; s<SizeOfSample;s++)
    {
        double* NormalSample= boxMuller(matrixSize,1.0,0.0);
        double runningsum=0;
        for(int i=0;i<matrixSize;i++)
        {
            runningsum=0;
            for(int j=0;j<matrixSize;j++)
            {
                runningsum+=NormalSample[j]*choldcMtrx[i][j];
            }
            MultiNormalSample[s][i]=runningsum+mean[i];
        }
        delete[] NormalSample;
    }
    /*out_Out<<"MultiNormalSample"<<endl;
     for(int s=0; s<SizeOfSample;s++)
     {
     for(int i=0;i<matrixSize;i++)
     {
     out_Out<<MultiNormalSample[s][i]<<",";
     }
     out_Out<<endl;
     }
     out_Out<<endl;*/
    //delete NormalSample;
    for(int i=0;i<matrixSize; i++)
    {
        delete[] choldcMtrx[i];
    }
    delete[] choldcMtrx;
    return MultiNormalSample;
}
double MatePref(double alpha,double z1,double z2)
{
    //alpha=p0.z[1]
    //traitval=p1.z[0]
    //traitopt=p0.z[0]
    double probmating=exp(-alpha*(pow(z1-z2,2)));
    return probmating;
}
void testChamber(individual* wrkPop,individual* pop0, individual* pop1,int* poptracker,int numInP0,int numInP1,int numOut,int numTraits,int numLoci,int g,int gphase,double alpha,ofstream &out_Test)
{
    //cout<<"testchamber"<<endl;
    //cout<<"numinp0 "<<numInP0<<" numinp1 "<<numInP1<<endl;
    int wrkPopcounter=0;
    int numMatings=0;
    int numHybridizations=0;
    int rInd=0;
    while(wrkPopcounter<(numOut*2))
    {
        //cout<<"position "<<wrkPopcounter<<endl;
        if(wrkPopcounter<(numOut))
        {
            rInd=rand()%numInP0;
            //cout<<"rind p0 "<<rInd<<endl;
            for(int t=0;t<numTraits;t++)
            {
                for(int l=0;l<numLoci;l++)
                {
                    wrkPop[wrkPopcounter].genome[t][l]=pop0[rInd].genome[t][l];
                }
            }
            poptracker[wrkPopcounter]=0;
            wrkPopcounter++;
        }
        else
        {
            rInd=rand()%numInP1;
            //cout<<"rind p1 "<<rInd<<endl;
            for(int t=0;t<numTraits;t++)
            {
                for(int l=0;l<numLoci;l++)
                {
                    wrkPop[wrkPopcounter].genome[t][l]=pop1[rInd].genome[t][l];
                }
            }
            poptracker[wrkPopcounter]=1;
            wrkPopcounter++;
        }
    }
    CalculatePhenotypesOfaPopulation(wrkPop, (numOut*2), numTraits, numLoci);
    //now do repro in pop
    int p0=0;
    int p1=0;
    while(numMatings<numOut)
    {
        p0=rand()%(numOut*2);
        p1=rand()%(numOut*2);
        //find probmating
        double probmating=(0.5*MatePref(alpha,wrkPop[p0].z[0],wrkPop[p1].z[0]))+(0.5*MatePref(alpha,wrkPop[p0].z[1],wrkPop[p1].z[1]));
        double prob=rand()/(double)RAND_MAX;
        //cout<<"probmating: "<<probmating<<", rand num: "<<prob<<endl;
        if(prob<probmating)
        {
            numMatings++;
            //if hybrids, inc numhybridizations
            if(poptracker[p0]!=poptracker[p1])
            {
                numHybridizations++;
            }
        }
    }
    out_Test<<gphase<<" "<<numMatings<<" "<<numHybridizations<<endl;
}
void writedmis(ofstream &bdmioutfile,int numBDMILoci,individual* popout,int numOut,int gentowrite,int poptowrite)
{
    for(int i=0;i<numOut;i++)
    {
        //print bdmi loci
        bdmioutfile<<gentowrite<<" p"<<poptowrite<<" ";
        for(int pset=0;pset<2;pset++)
        {
            for(int b=0;b<numBDMILoci;b++)
            {
                bdmioutfile<<popout[i].bdmi[pset][b]<<" ";
            }
        }
        bdmioutfile<<endl;
    }
    bdmioutfile<<endl;
}
