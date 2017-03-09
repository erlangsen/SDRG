#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include "r250.c"
//#include "mpi.h"

struct JBond{
    int Jleft;
    int Jright;
    double jj; //value
};
typedef struct JBond j;

struct Q3Bond{
    double qq; //value
    int Qk; //1
    int Ql; //2
    int Qm; //3
    int Qn; //4
    int Qo; //5
    int Qp; //6
};
typedef struct Q3Bond q3;

struct Q2Bond{
    double qq; //value
    int Q2k; //1
    int Q2l; //2
    int Q2m; //3
    int Q2n; //4
};
typedef struct Q2Bond q2;

struct Point{
    int left1;
    int left2;
    int left3;
    int left4;
    int left5;
    int right1;
    int right2;
    int right3;
    int right4;
    int right5;
};
typedef struct Point p;


int DJ=0;
const double Jmin =20.0;
const double Bondmin = -1000000000000000000.0;

int L;  //length of chain
int N;  //N samples
int livpoint;
int Check;
int Write;
int Gap;
int Bins;
int EE;
int CCheck;
int Num;
int maxxtype;
double D; //disorder strenth
double localJ,localQ;
double Esum,Emax,maxx;
int *singlet,*act,*Nc,*Nd,*Sublattice,*Lambda,*Oo,*siteA,*loopsize,*Vr,*Vl;
double *cr,*dr,*drsta,*load1,*load2,*drcopy,*lncr,*lndrsta,*dr2,*drcopy2,*drsta2,*dr3,*drcopy3,*drsta3,*Q3;
j *bondj;
q3 *bondq3;
q2 *bondq2;
p *neighbor;
int type,bondindex;
char name[200],name2[200],name3[200],name4[200],name5[200],name6[200],name7[200],name8[200],name9[200],name10[200];
int counter;
int *pair;
int *showpair;

void initialize(){
    
    for (int i=0;i<L;i++){
        
        //set initial Jbond and name
        bondj[i].Jleft=i%L;
        bondj[i].Jright=(i+1)%L;
        
        //if (DJ==0)
           bondj[i].jj=Jmin;
        //else
        //double r=double_r250();      // uniform distribution (0,1)

          //  bondj[i].jj=log(r)*D;
        
        //set initial Qbond and name

        bondq3[i].Qk=i;
        bondq3[i].Ql=(i+1)%L;
        bondq3[i].Qm=(i+2)%L;
        bondq3[i].Qn=(i+3)%L;
        bondq3[i].Qo=(i+4)%L;
        bondq3[i].Qp=(i+5)%L;

        
        double r=double_r250();      // uniform distribution (0,1)
        double QQ=pow(r, D);       // d=1-> uniform; d>1 "stronger disorder"
        bondq3[i].qq=log(r)*D;
        Q3[i]=bondq3[i].qq;
       // bondq3[i].qq=Jmin;

        
        bondq2[i].Q2k=i;
        bondq2[i].Q2l=(i+1)%L;
        bondq2[i].Q2m=(i+2)%L;
        bondq2[i].Q2n=(i+3)%L;
        bondq2[i].qq=Jmin;
        
        
        neighbor[i].left1=(i-1+L)%L;
        neighbor[i].left2=(i-2+L)%L;
        neighbor[i].left3=(i-3+L)%L;
        neighbor[i].left4=(i-4+L)%L;
        neighbor[i].left5=(i-5+L)%L;

        neighbor[i].right1=(i+1)%L;
        neighbor[i].right2=(i+2)%L;
        neighbor[i].right3=(i+3)%L;
        neighbor[i].right4=(i+4)%L;
        neighbor[i].right5=(i+5)%L;

        singlet[i]=0;
        act[i]=1;
        pair[i]=-1;
        showpair[i]=-1;
        
        drcopy[i]=0;
        drcopy2[i]=0;
        drcopy3[i]=0;
        
        Sublattice[i]=0;
        Lambda[i]=0;
        Oo[i]=0;
        siteA[i]=0;
        loopsize[i]=0;
        Vr[i]=0;
        Vl[i]=0;
    }
}

void allocateArray(){
    
    bondj =new j[L];
    bondq3 =new q3[L];
    bondq2= new q2[L];
    neighbor = new p[L];
    load1 = new double[6];
    load2 = new double[6];
    act = new int[L];
    singlet = new int[L];
    cr = new double[L];
    Q3 = new double[L];
    
    dr = new double[L];
    dr2 = new double[L];
    dr3 = new double[L];

    drsta = new double[L];
    drsta2 = new double[L];
    drsta3 = new double[L];

    drcopy = new double[L];
    drcopy2 = new double[L];
    drcopy3 = new double[L];

    pair = new int[L];
    showpair = new int[L];
    Nd = new int[L];
    Nc = new int[L];
    lncr = new double[L];
    lndrsta = new double[L];
    Sublattice = new int[L];
    Lambda = new int[L];
    Oo = new int[L];
    siteA = new int[L];
    loopsize = new int[L];
    Vr = new int[L];
    Vl = new int[L];
}

void resetArray(){
    
    for (int i=0; i<L; i++) {
        dr[i]=0.0;
        dr2[i]=0.0;
        dr3[i]=0.0;
        drsta[i]=0.0;
        drsta2[i]=0.0;
        drsta3[i]=0.0;
        singlet[i]=-1;
        cr[i]=0.0;
        act[i]=0;
        drcopy[i]=0.0;
        drcopy2[i]=0.0;
        drcopy3[i]=0.0;
        lncr[i]=0.0;
        lndrsta[i]=0.0;
        Nd[i]=0;
        Nc[i]=0;
        Sublattice[i]=0;
        Lambda[i]=0;
        Oo[i]=0;
        siteA[i]=0;
        loopsize[i]=0;
        Vr[i]=0;
        Vl[i]=0;
    }
}

void freeArray(){
    
    delete [] bondj;
    delete [] bondq3;
    delete [] bondq2;
    delete [] singlet;
    delete [] cr;
    delete [] dr;
    delete [] load1;
    delete [] load2;
    delete [] drsta;
    delete [] drcopy;
    delete [] act;
    delete [] neighbor;
    delete [] pair;
    delete [] showpair;
    delete [] lncr;
    delete [] lndrsta;
    delete [] Nd;
    delete [] Nc;
    delete [] Q3;

}

void resetVariable(){

     livpoint = L;
     type = 0;
     bondindex = 0;
     localJ = 0.0;
     localQ = 0.0;
     counter= 0;
}

double compare(double a,double b){
    
    if (a>b)
        return a;
    else
        return b;
}

double max(){
    
    double x,y,o,p,s,t;
    int xindex,oindex;
    
    x=Bondmin;
    
    for (int i=0; i<L; i++) {
        if (bondj[i].jj<=0){
            y=bondj[i].jj;
            x=compare(x,y);
        }
    }
    
    o=Bondmin;
    
    for (int i=0; i<L; i++) {
        if (bondq3[i].qq<=0){
            p=bondq3[i].qq;
            o=compare(o,p);
        }
    }
    
    s=Bondmin;
    
    for (int i=0; i<L; i++) {
        if (bondq2[i].qq<=0){
            t=bondq2[i].qq;
            s=compare(t,s);
        }
    }
    
    if (compare(compare(x,o),s)==x){
        for (int i=0; i<L; i++) {
            if (bondj[i].jj==x){
                type = 0;//Q1
                bondindex = i;
            }
        }
    }
    
    else if (compare(compare(x,o),s)==o){
        for (int i=0; i<L; i++) {
            if (bondq3[i].qq==o){
                type = 2;//Q3
                bondindex = i;
            }
        }
    }
    
    else if (compare(compare(x,o),s)==s){
        for (int i=0; i<L; i++) {
            if (bondq2[i].qq==s){
                type = 1;//Q2
                bondindex = i;
            }
        }
    }
    return compare(compare(x,o),s);
}

void AddE(){
    
    int NBonds=0;
    Esum=0.0;
    
    for (int i=0; i<L; i++) {
        if (bondj[i].jj!=Jmin) {
            Esum=Esum+exp(bondj[i].jj);
            NBonds++;
        }
        if (bondq2[i].qq!=Jmin) {
            Esum=Esum+exp(bondq2[i].qq);
            NBonds++;
        }
        if (bondq3[i].qq!=Jmin) {
            Esum=Esum+exp(bondq3[i].qq);
            NBonds++;
        }
    }
    Esum=Esum/NBonds;
    Esum=log(Esum);
}




void rgq1(){
    
    //save Jmax's informations
    double finalJ=0.0;
    int a=bondindex;
    int Jk=bondj[a].Jleft;
    int Jl=bondj[a].Jright;
    double Jvalue=bondj[a].jj;
    
    act[Jk]=act[Jl]=0;
    
    livpoint=livpoint-2;
    
    //load cr
    int r=abs(bondj[a].Jright-bondj[a].Jleft);
    
    //prepare dr
    if (r==1){
        int small = (int)compare(Jk,Jl)-1;
        singlet[small]=1;// small is a r=1 singlet index ,dimer[small]=1 ,else =-1
    }
    
    else if(r==(L-1)) //(0,l-1) or (l-1,0)
        singlet[L-1]=1;
    
    
    pair[Jk]=counter;
    pair[Jl]=counter;
    showpair[Jk]=Jl;
    showpair[Jl]=Jk;
    
    counter=counter+1;
    
    //save neighbor(point)
    int qleft1=neighbor[Jk].left1;
    int qleft2=neighbor[Jk].left2;
    int qleft3=neighbor[Jk].left3;
    int qleft4=neighbor[Jk].left4;
    int qleft5=neighbor[Jk].left5;

    int qright1=neighbor[Jl].right1;
    int qright2=neighbor[Jl].right2;
    int qright3=neighbor[Jl].right3;
    int qright4=neighbor[Jl].right4;
    int qright5=neighbor[Jl].right5;

    
    //save cut Q2
    double jj1=Bondmin;
    double jj2=Bondmin;
    
    double leftq3a=Bondmin;
    double leftq3b=Bondmin;
    double leftq3c=Bondmin;
    double rightq3a=Bondmin;
    double rightq3b=Bondmin;
    double rightq3c=Bondmin;

    int effjleft=-1;
    int effjright=-1;
    
    /*
     ------------------------
            Jk      Jl
            ==========
     ==jj1====      ====jj2==
     ------------------------
    */
    
    //find J's neighborJ
    for (int i=0; i<L; i++) {
        if(bondj[i].Jright==Jk){
            jj1=bondj[i].jj;
            effjleft=i;
        }
        
        if(bondj[i].Jleft==Jl){
            jj2=bondj[i].jj;
            effjright=i;
        }
    }
    
    //final
    int final=0;
    if (jj1==jj2&&effjleft==effjright){
        if (jj1!=Jmin){
            final=1;
            bondj[effjleft].jj=Jmin;
            bondj[effjleft].Jleft=-1;
            bondj[effjleft].Jright=-1;
        }
    }

    if (final) return;
    
    
    /*
    ----------------------------------------------------------------------------------------------------------------------------------
    left5       left4       left3       left2       left1       Jk      Jl      right1      right2      right3      right4      right5
                                                                ==========
    ===========================leftq3c============================      ============================rightq3c==========================
                ========================leftq3b===========================
                                                                =========================rightq3b=========================
                            =======================leftq3a============================
                                                    =======================rightq3a===========================
    ----------------------------------------------------------------------------------------------------------------------------------
    */
    
    
    //kill q3bond and save neighbors
    for (int i=0;i<L;i++) {
        if((bondq3[i].Qk==Jk||bondq3[i].Qk==Jl||bondq3[i].Ql==Jk||bondq3[i].Ql==Jl||bondq3[i].Qm==Jk||bondq3[i].Qm==Jl||
            bondq3[i].Qn==Jk||bondq3[i].Qn==Jl||bondq3[i].Qo==Jk||bondq3[i].Qo==Jl||bondq3[i].Qp==Jk||bondq3[i].Qp==Jl)
           &&
           (bondq3[i].qq!=Jmin))
        {
            if ((bondq3[i].Qo==qleft1)&&(bondq3[i].Qp==Jk)){ //cutQ2&per
                leftq3c=bondq3[i].qq;

                if (bondq2[i].qq==Jmin)
                    bondq2[i].qq=leftq3c-log(4);
                
                else
                    bondq2[i].qq=compare(bondq2[i].qq,leftq3c-log(4));
                
                
            }
            
            if ((bondq3[i].Qo==Jk)&&(bondq3[i].Qp==Jl)){ //cutQ2
                if (bondq2[i].qq==Jmin)
                    bondq2[i].qq=bondq3[i].qq;
                else
                    bondq2[i].qq=compare(bondq2[i].qq,bondq3[i].qq);

                leftq3b=bondq3[i].qq;
            }
            
            if ((bondq3[i].Qo==Jl)&&(bondq3[i].Qp==qright1)) //special case&cutJ
                leftq3a=bondq3[i].qq;
            
            
            if ((bondq3[i].Qk==qleft1)&&(bondq3[i].Ql==Jk)) //special case&cutJ
                rightq3a=bondq3[i].qq;
            
            if ((bondq3[i].Qk==Jk)&&(bondq3[i].Ql==Jl)){ //cutQ2
                if (bondq2[neighbor[i].right2].qq==Jmin)
                    bondq2[neighbor[i].right2].qq=bondq3[i].qq;
                else
                    bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq,bondq3[i].qq);
                
                rightq3b=bondq3[i].qq;
            }

            if ((bondq3[i].Qk==Jl)&&(bondq3[i].Ql==qright1)){ //cutQ2&per
                rightq3c=bondq3[i].qq;

                if (bondq2[neighbor[i].right2].qq==Jmin)
                    bondq2[neighbor[i].right2].qq=rightq3c-log(4);
                else
                    bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq,rightq3c-log(4));
                
            }
            
            bondq3[i].Qk=-1;
            bondq3[i].Ql=-1;
            bondq3[i].Qm=-1;
            bondq3[i].Qn=-1;
            bondq3[i].Qo=-1;
            bondq3[i].Qp=-1;
            bondq3[i].qq=Jmin;
        }
    }
    //printf("left1q1=%lf,left1q2=%lf,left1q3=%lf,left2q=%lf,left3q=%lf,right1q1=%lf,right1q2=%lf,right1q3=%lf,rightq2=%lf,righq3=%lf\n",left1q1,left1q2,left1q3,left2q,left3q,right1q1,right1q2,right1q3,right2q,right3q);
    /*
    ----------------------------------------------------------------------------------
    left3       left2       left1       Jk      Jl      right1      right2      right3
                                        ==========
    ==============leftq2c=================      =================rightq2c=============
                ===============leftq2b============
                                        ==============rightq2b============
                            =======leftq2a====rightq2a========
    ----------------------------------------------------------------------------------
    
    */
    //kill q2bond and save neighbors
    double leftq2a=Bondmin;
    double leftq2b=Bondmin;
    double leftq2c=Bondmin;
    double rightq2a=Bondmin;
    double rightq2b=Bondmin;
    double rightq2c=Bondmin;
    
    for (int i=0;i<L;i++) {
        if((bondq2[i].Q2k==Jk||bondq2[i].Q2k==Jl||
            bondq2[i].Q2l==Jk||bondq2[i].Q2l==Jl||
            bondq2[i].Q2m==Jk||bondq2[i].Q2m==Jl||
            bondq2[i].Q2n==Jk||bondq2[i].Q2n==Jl)&&
           (bondq2[i].qq!=Jmin))
        {
            if((bondq2[i].Q2m==qleft1)&&(bondq2[i].Q2n==Jk)){ //cutJ&per
                leftq2c=bondq2[i].qq;
            }
            
            if((bondq2[i].Q2m==Jk)&&(bondq2[i].Q2n==Jl)){ //cutJ
                leftq2b=bondq2[i].qq;
            }
            
            if((bondq2[i].Q2m==Jl)&&(bondq2[i].Q2n==qright1)){ //special case
                leftq2a=bondq2[i].qq;
            }
            
            if((bondq2[i].Q2k==qleft1)&&(bondq2[i].Q2l==Jk)){ //special case
                rightq2a=bondq2[i].qq;
            }
            
            if((bondq2[i].Q2k==Jk)&&(bondq2[i].Q2l==Jl)){ //cutJ
                rightq2b=bondq2[i].qq;
            }
            
            if((bondq2[i].Q2k==Jl)&&(bondq2[i].Q2l==qright1)){ //cutJ&per
                rightq2c=bondq2[i].qq;
            }
           
            bondq2[i].Q2k=-1;
            bondq2[i].Q2l=-1;
            bondq2[i].Q2m=-1;
            bondq2[i].Q2n=-1;
            bondq2[i].qq=Jmin;
        }
    }
    
    double effqq1=Bondmin;
    double effqq2=Bondmin;
    double effJ=Bondmin;
    double specialJ1=Bondmin; //from Q2
    double specialJ2=Bondmin; //from Q3
    double specialJ=Bondmin;

   
    load1[0]=jj1;
    load1[1]=leftq2c;
    load1[2]=leftq3c;
    
    load2[0]=jj2;
    load2[1]=rightq2c;
    load2[2]=rightq3c;
    
    if ((rightq2a==leftq2a)&&(rightq2a!=Bondmin))
        specialJ1=rightq2a-log(4);
    
    if ((leftq3a!=Bondmin)||(rightq3a!=Bondmin))
        specialJ2=compare(leftq3a,rightq3a)-log(4);
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3;j++){
            if ((load1[i]!=Jmin)&&(load2[j]!=Jmin)){
                effJ=compare(effJ,load1[i]+load2[j]-log(2)-Jvalue);
            }
        }
    }

    //create new bond and kill
    bondj[a].jj=compare(effJ,compare(specialJ1,specialJ2));
    
    //printf("effJ=%lf,specialJ1=%lf,specialJ2=%lf\n",effJ,specialJ1,specialJ2);

    
    if (bondj[a].jj<=Bondmin)
        //bondj[a].jj=Jmin;
        bondj[a].jj=Bondmin;
        
    
    bondj[a].Jleft=bondj[effjleft].Jleft;
    int s1=bondj[effjleft].Jleft;
    bondj[a].Jright=bondj[effjright].Jright;
    int s2=bondj[effjright].Jright;
        
    bondj[effjleft].Jleft=-1;  bondj[effjleft].Jright=-1;  bondj[effjleft].jj=Jmin;
    bondj[effjright].Jleft=-1; bondj[effjright].Jright=-1; bondj[effjright].jj=Jmin;
    
    /*
    //////
    if (bondj[a].jj==specialJ) {
        bondj[neighbor[effjleft].left2].jj=Jmin;
        bondj[neighbor[effjleft].left2].Jleft=-1;
        bondj[neighbor[effjleft].left2].Jright=-1;
        
        bondj[neighbor[effjright].right1].jj=Jmin;
        bondj[neighbor[effjright].right1].Jleft=-1;
        bondj[neighbor[effjright].right1].Jright=-1;
    }
    //////
    */
    
    int ss1=0;
    int ss2=0;
   // printf("s1=%d,s2=%d\n",s1,s2);
    
    //Check s1&s2
    for (int i=0; i<L;i++) {
        if (bondj[i].jj!=Jmin) {
            if(bondj[i].Jright==s1 || bondj[i].Jleft==s1)
                ss1=ss1+1;
            else
                ss1=ss1;
            
            if(bondj[i].Jright==s2 || bondj[i].Jleft==s2)
                ss2=ss2+1;
            else
                ss2=ss2;
        }
        
        if (bondq2[i].qq!=Jmin) {
            if (bondq2[i].Q2k==s1 || bondq2[i].Q2l==s1 ||bondq2[i].Q2m==s1 ||bondq2[i].Q2n==s1 )
                ss1=ss1+1;
            else
                ss1=ss1;
            
            if (bondq2[i].Q2k==s2 || bondq2[i].Q2l==s2 ||bondq2[i].Q2m==s2 ||bondq2[i].Q2n==s2 )
                ss2=ss2+1;
            else
                ss2=ss2;
        }
        
        if (bondq3[i].qq!=Jmin) {
            if (bondq3[i].Qk==s1 || bondq3[i].Ql==s1 ||bondq3[i].Qm==s1 || bondq3[i].Qn==s1 || bondq3[i].Qo==s1 || bondq3[i].Qp==s1 )
                ss1=ss1+1;
            else
                ss1=ss1;
            
            if (bondq3[i].Qk==s2 || bondq3[i].Ql==s2 ||bondq3[i].Qm==s2 || bondq3[i].Qn==s2 || bondq3[i].Qo==s2 || bondq3[i].Qp==s2 )
                ss2=ss2+1;
            else
                ss2=ss2;
        }
    }
    
    if (ss1==0){
        if (act[s1])
            livpoint=livpoint-1;
        act[s1]=0;
    }
    
    if (ss2==0) {
        if (act[s2])
            livpoint=livpoint-1;
        act[s2]=0;
    }
    
    
    

    //add cut Q1
    for (int i=0; i<L; i++) {
       
        if (((bondj[i].Jleft==qleft3)&&(bondj[i].Jright==qleft2))&&((act[qleft3]==act[qleft2])&&(act[qleft2]==1))){
            
          /*  if (leftq3a!=Bondmin){
                leftq3a=leftq3a-log(16);

                if (bondj[i].jj==Jmin)
                    bondj[i].jj=leftq3a;
                else
                    bondj[i].jj=compare(bondj[i].jj,leftq3a);
            }*/
        
            if (leftq2c!=Bondmin) {
                leftq2c=leftq2c-log(4);
                
                if (bondj[i].jj==Jmin)
                    bondj[i].jj=leftq2c;
                else
                    bondj[i].jj=compare(bondj[i].jj,leftq2c);
            }
        }
        
        if (((bondj[i].Jleft==qleft2)&&(bondj[i].Jright==qleft1))&&((act[qleft2]==act[qleft1])&&(act[qleft1]==1))&&(leftq2b!=Bondmin)){
            
            if (bondj[i].jj==Jmin)
                bondj[i].jj=leftq2b;
            else
                bondj[i].jj=compare(bondj[i].jj,leftq2b);
        }
       
        if (((bondj[i].Jleft==qright2)&&(bondj[i].Jright==qright3))&&((act[qright2]==act[qright3])&&(act[qright3]==1))){
            
           /* if(rightq3a!=Bondmin){
                rightq3a=rightq3a-log(16);
               
                if (bondj[i].jj==Jmin)
                    bondj[i].jj=rightq3a;
                else
                    bondj[i].jj=compare(bondj[i].jj,rightq3a);
            }*/
        
            if (rightq2c!=Bondmin){
                rightq2c=rightq2c-log(4);

                if (bondj[i].jj==Jmin)
                    bondj[i].jj=rightq2c;
                else
                    bondj[i].jj=compare(bondj[i].jj,rightq2c);
            }
        }
        
        if (((bondj[i].Jleft==qright1)&&(bondj[i].Jright==qright2))&&((act[qright1]==act[qright2])&&(act[qright2]==1))&&(rightq2b!=Bondmin)){
            
            if (bondj[i].jj==Jmin)
                bondj[i].jj=rightq2b;
            else
                    bondj[i].jj=compare(bondj[i].jj,rightq2b);
        }
    }
  /*
  
    //kill bond like (a,b)&(b,a)
    for (int i=0;i<L;i++){
        if ((bondj[a].Jleft==bondj[i].Jright)&&(bondj[a].Jright==bondj[i].Jleft)&&(bondj[a].jj!=Jmin)) {
            bondj[a].Jright=-1;
            bondj[a].Jleft=-1;
            bondj[a].jj=Jmin;
           
        }
    }
   */
   
    int last=-1;
    int rr=-1;
    
    if(livpoint==2){
        
        for (int i=0;i<L;i++){
            if (bondj[i].jj!=Jmin)
                last=i;
        }

        if(last!=-1){
            rr=abs(bondj[last].Jright-bondj[last].Jleft);
            
            //prepare dr
            if (rr==1){
                int small = (int)compare(bondj[last].Jright,bondj[last].Jleft)-1;
                singlet[small]=1;// small is a r=1 singlet index ,dimer[small]=1 ,else =-1
            }
            else if(rr==(L-1)) //(0,l-1) or (l-1,0)
                singlet[L-1]=1;
            
            pair[bondj[last].Jright]=counter;
            pair[bondj[last].Jleft]=counter;
            
            showpair[bondj[last].Jright]=bondj[last].Jleft;
            showpair[bondj[last].Jleft]=bondj[last].Jright;
            
            counter=counter+1;
        }
    }
}


void rgq2(){
    
    double finalQ=0.0;
    int aa=bondindex;
    int s1,s2,s3,s4;
  
    s1 = bondq2[aa].Q2k;
    s2 = bondq2[aa].Q2l;
    s3 = bondq2[aa].Q2m;
    s4 = bondq2[aa].Q2n;
    //printf("(%d,%d,%d,%d)\n",s1,s2,s3,s4);
    
    double qqq = bondq2[aa].qq;
    
    act[s1]=0; act[s2]=0; act[s3]=0; act[s4]=0;
    
    //load dimer
    singlet[s1]=1; singlet[s3]=1;
    
    //save cr
    pair[s1]=counter; pair[s2]=counter; counter=counter+1;
    pair[s3]=counter; pair[s4]=counter; counter=counter+1;
    
    showpair[s1]=s2;
    showpair[s2]=s1;
    showpair[s3]=s4;
    showpair[s4]=s3;
    
    if ( ((livpoint==4)&&(type==0)) ) {
        for (int i=0; i<L; i++) {
            bondj[i].Jleft=-1;
            bondj[i].Jright=-1;
            bondj[i].jj=Jmin;
        }
        localQ=qqq;
        livpoint=0;
    }
    
    else{
        
        double jj1=Bondmin;
        double jj2=Bondmin;
        double leftq3a=Bondmin;
        double leftq3b=Bondmin;
        double leftq3c=Bondmin;
        double leftq3d=Bondmin;
        double leftq3e=Bondmin;
        double rightq3a=Bondmin;
        double rightq3b=Bondmin;
        double rightq3c=Bondmin;
        double rightq3d=Bondmin;
        double rightq3e=Bondmin;

        //save neighbor(point)
        int qleft1=neighbor[s1].left1;
        int qleft2=neighbor[s1].left2;
        int qleft3=neighbor[s1].left3;
        int qleft4=neighbor[s1].left4;
        int qleft5=neighbor[s1].left5;
        
        int qright1=neighbor[s4].right1;
        int qright2=neighbor[s4].right2;
        int qright3=neighbor[s4].right3;
        int qright4=neighbor[s4].right4;
        int qright5=neighbor[s4].right5;
        
        
        //choosed Q's neighbors
        
        int qleft=-1;
        int qright=-1;
        
        for(int i=0;i<L;i++){
            if((bondj[i].Jleft ==s1||bondj[i].Jleft ==s2||bondj[i].Jleft ==s3||bondj[i].Jleft ==s4||
                bondj[i].Jright==s1||bondj[i].Jright==s2||bondj[i].Jright==s3||bondj[i].Jright==s4))
            {
                if (bondj[i].Jright==s1){
                    qleft=i;
                    jj1=bondj[i].jj;
                }
                
                if (bondj[i].Jleft==s4){
                    qright=i;
                    jj2=bondj[i].jj;
                }
                bondj[i].jj=Jmin;
            }
        }
        
        //printf("qleft=%d,qright=%d\n",qleft,qright);
        //printf("jj1=%lf,jj2=%lf\n",jj1,jj2);
        
        /*
         --------------------------------------------------------------------------------------------------------------------------------------------------------
         left5      left4       left3      left2       left1        s1       s2     s3      s4      right1      right2      right3      right4      right5
                                                                    ==========================
         =========================leftq3e=============================                      ============================rightq3e==========================
                    ========================leftq3d============================     ==========================rightq3d========================
                               ===============leftq3c=================================
                                                                             ===========================rightq3c==================
                                           =======================leftq3b=====================
                                                                    =====================rightq3b=====================
                                                        ================leftq3a=====rightq3a==============
         --------------------------------------------------------------------------------------------------------------------------------------------------------
        */
        
        //kill q3bond and save neighbors
        for (int i=0;i<L;i++) {
            if((bondq3[i].Qk==s1||bondq3[i].Qk==s2||bondq3[i].Qk==s3||bondq3[i].Qk==s4||
                bondq3[i].Ql==s1||bondq3[i].Ql==s2||bondq3[i].Ql==s3||bondq3[i].Ql==s4||
                bondq3[i].Qm==s1||bondq3[i].Qm==s2||bondq3[i].Qm==s3||bondq3[i].Qm==s4||
                bondq3[i].Qn==s1||bondq3[i].Qn==s2||bondq3[i].Qn==s3||bondq3[i].Qn==s4||
                bondq3[i].Qo==s1||bondq3[i].Qo==s2||bondq3[i].Qo==s3||bondq3[i].Qo==s4||
                bondq3[i].Qp==s1||bondq3[i].Qp==s2||bondq3[i].Qp==s3||bondq3[i].Qp==s4
                )&&
               (bondq3[i].qq!=Jmin))
            {
                
                
                if ((bondq3[i].Qk==qleft1)&&(bondq3[i].Ql==s1)) //special case
                    leftq3a=bondq3[i].qq;
                
                if ((bondq3[i].Qk==qleft2)&&(bondq3[i].Ql==qleft1)) //cutJ
                    leftq3b=bondq3[i].qq;
                
                if ((bondq3[i].Qk==qleft3)&&(bondq3[i].Ql==qleft2)) //cutJ&per
                    leftq3c=bondq3[i].qq;
                
                if ( (bondq3[i].Qk==qleft4)&&(bondq3[i].Ql==qleft3)&&(bondq3[i].Qm==qleft2)&&(bondq3[i].Qn==qleft1) ){ //cutQ2
                    if (bondq2[i].qq==Jmin)
                        bondq2[i].qq=bondq3[i].qq;
                    else
                        bondq2[i].qq=compare(bondq2[i].qq,bondq3[i].qq);
                    
                    leftq3d=bondq3[i].qq;
                }
                
                if ( (bondq3[i].Qk==qleft5)&&(bondq3[i].Ql==qleft4)&&(bondq3[i].Qm==qleft3)&&(bondq3[i].Qn==qleft2) ){ //cutQ2&per(maybe)
                    leftq3e=bondq3[i].qq;

                    if (bondq2[i].qq==Jmin)
                        bondq2[i].qq=leftq3e-log(4);
                    else
                        bondq2[i].qq=compare(bondq2[i].qq,leftq3e-log(4));
                    
                }
                
                if ((bondq3[i].Qo==s4)&&(bondq3[i].Qp==qright1)) //special case
                    rightq3a=bondq3[i].qq;
                
                if ((bondq3[i].Qo==qright1)&&(bondq3[i].Qp==qright2)) //cutJ
                    rightq3b=bondq3[i].qq;
                
                if ((bondq3[i].Qo==qright2)&&(bondq3[i].Qp==qright3)) //cutJ&per
                    rightq3c=bondq3[i].qq;
                
                if ( (bondq3[i].Qm==qright1)&&(bondq3[i].Qn==qright2)&&(bondq3[i].Qo==qright3)&&(bondq3[i].Qp==qright4) ){ //cutQ2
                    if (bondq2[neighbor[i].right2].qq==Jmin)
                        bondq2[neighbor[i].right2].qq=bondq3[i].qq;
                    else
                        bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq,bondq3[i].qq);
                    
                    rightq3d=bondq3[i].qq;
                }
                
                if ( (bondq3[i].Qm==qright2)&&(bondq3[i].Qn==qright3)&&(bondq3[i].Qo==qright4)&&(bondq3[i].Qp==qright5) ){ //cutQ2&per(maybe)
                    rightq3e=bondq3[i].qq;

                    if (bondq2[neighbor[i].right2].qq==Jmin)
                        bondq2[neighbor[i].right2].qq= rightq3e-log(4);
                    else
                        bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq, rightq3e-log(4));
                    
                }

                
                
                bondq3[i].Qk=-1;
                bondq3[i].Ql=-1;
                bondq3[i].Qm=-1;
                bondq3[i].Qn=-1;
                bondq3[i].Qo=-1;
                bondq3[i].Qp=-1;
                bondq3[i].qq=Jmin;
            }
        }
        //printf("left1q1=%lf,left1q2=%lf,left1q3=%lf,left2q=%lf,left3q=%lf,right1q1=%lf,right1q2=%lf,right1q3=%lf,rightq2=%lf,righq3=%lf\n",left1q1,left1q2,left1q3,left2q,left3q,right1q1,right1q2,right1q3,right2q,right3q);
        
        /*
        -----------------------------------------------------------------------------------------------
        left3      left2      left1      s1      s2      s3      s4      right1      right2      right3
                                         ==========================
        ============leftq2c================                      ==============rightq2c================
                   ===============leftq2b==========      ==============right2b=============
                              =============leftq2a=========
                                                 =========right2a==============
        -----------------------------------------------------------------------------------------------
         
        */
        //kill q2bond and save neighbors
        double leftq2a=Bondmin;
        double leftq2b=Bondmin;
        double leftq2c=Bondmin;
        double rightq2a=Bondmin;
        double rightq2b=Bondmin;
        double rightq2c=Bondmin;

        for (int i=0;i<L;i++) {
            if((bondq2[i].Q2k==s1||bondq2[i].Q2k==s2||bondq2[i].Q2k==s3||bondq2[i].Q2k==s4||
                bondq2[i].Q2l==s1||bondq2[i].Q2l==s2||bondq2[i].Q2l==s3||bondq2[i].Q2l==s4||
                bondq2[i].Q2m==s1||bondq2[i].Q2m==s2||bondq2[i].Q2m==s3||bondq2[i].Q2m==s4||
                bondq2[i].Q2n==s1||bondq2[i].Q2n==s2||bondq2[i].Q2n==s3||bondq2[i].Q2n==s4)&&
               (bondq2[i].qq!=Jmin))
            {
                if((bondq2[i].Q2k==qleft1)&&(bondq2[i].Q2l==s1)){ //per
                    leftq2a=bondq2[i].qq;
                }
                
                if((bondq2[i].Q2k==qleft2)&&(bondq2[i].Q2l==qleft1)){ //cutJ
                    leftq2b=bondq2[i].qq;
                }
                
                if((bondq2[i].Q2k==qleft3)&&(bondq2[i].Q2l==qleft2)){ //cutJ&per(maybe)
                    leftq2c=bondq2[i].qq;
                }
                
                if((bondq2[i].Q2k==s2)&&(bondq2[i].Q2l==s3)){ //per
                    rightq2a=bondq2[i].qq;
                }
                
                if((bondq2[i].Q2k==s3)&&(bondq2[i].Q2l==s4)){ //cutJ
                    rightq2b=bondq2[i].qq;
                }
                
                if((bondq2[i].Q2k==s4)&&(bondq2[i].Q2l==qright1)){ //cutJ&per(maybe)
                    rightq2c=bondq2[i].qq;
                }
                
                bondq2[i].Q2k=-1;
                bondq2[i].Q2l=-1;
                bondq2[i].Q2m=-1;
                bondq2[i].Q2n=-1;
                bondq2[i].qq=Jmin;
            }
        }

        

        
        double effqq1=Bondmin;
        double effqq2=Bondmin;
        double specialJ=Bondmin;
        double effJ=Bondmin;
        
        if ((leftq3a==rightq3a)&&(rightq3a!=Bondmin))
            specialJ=rightq3a-log(16);
        
        load1[0]=jj1;
        load1[1]=leftq2c;
        load1[2]=leftq3e;
        load1[3]=leftq2a;
        load1[4]=leftq3c;
        
        load2[0]=jj2;
        load2[1]=rightq2c;
        load2[2]=rightq3e;
        load2[3]=rightq2a;
        load2[4]=rightq3c;
        
        for (int i=0; i<5; i++) {
            
            if (i<=2) {
                for (int j=3; j<5; j++){
                    if ((load1[i]!=Jmin)&&(load2[j]!=Jmin))
                        effJ=compare(effJ,load1[i]+load2[j]-log(8)-qqq);
                }
            }
            
            else if(i>2){
                for (int j=0; j<5; j++){
                    if ((load1[i]!=Jmin)&&(load2[j]!=Jmin))
                        effJ=compare(effJ,load1[i]+load2[j]-log(8)-qqq);
                }
            }
        }
        
        //create jbond for pbc
        bondj[qleft].jj=compare(effJ,specialJ);
        //printf("effJ=%lf,specialJ=%lf\n",effJ,specialJ);
        
        if (bondj[qleft].jj<=Bondmin) {
            bondj[qleft].jj=Bondmin;
        }
        
        
        bondj[qleft].Jleft=bondj[qleft].Jleft;
        int s1=bondj[qleft].Jleft;
        bondj[qleft].Jright=bondj[qright].Jright;
        int s2=bondj[qright].Jright;
        
        bondj[qright].Jleft=-1; bondj[qright].Jright=-1; bondj[qright].jj=Jmin;
        
       /* //////
        if (bondj[qleft].jj==specialJ) {
            bondj[neighbor[qleft].left2].jj=Jmin;
            bondj[neighbor[qleft].left2].Jleft=-1;
            bondj[neighbor[qleft].left2].Jright=-1;

            bondj[neighbor[qright].right1].jj=Jmin;
            bondj[neighbor[qright].right1].Jleft=-1;
            bondj[neighbor[qright].right1].Jright=-1;
        }
        
        */
        
        livpoint=livpoint-4;
        
        int ss1=0;
        int ss2=0;
        
        //Check s1&s2
        for (int i=0; i<L;i++) {
            if (bondj[i].jj!=Jmin) {
                if(bondj[i].Jright==s1 || bondj[i].Jleft==s1)
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if(bondj[i].Jright==s2 || bondj[i].Jleft==s2)
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
            
            if (bondq2[i].qq!=Jmin) {
                if (bondq2[i].Q2k==s1 || bondq2[i].Q2l==s1 ||bondq2[i].Q2m==s1 ||bondq2[i].Q2n==s1 )
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if (bondq2[i].Q2k==s2 || bondq2[i].Q2l==s2 ||bondq2[i].Q2m==s2 ||bondq2[i].Q2n==s2 )
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
            
            if (bondq3[i].qq!=Jmin) {
                if (bondq3[i].Qk==s1 || bondq3[i].Ql==s1 ||bondq3[i].Qm==s1 || bondq3[i].Qn==s1 || bondq3[i].Qo==s1 || bondq3[i].Qp==s1 )
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if (bondq3[i].Qk==s2 || bondq3[i].Ql==s2 ||bondq3[i].Qm==s2 || bondq3[i].Qn==s2 || bondq3[i].Qo==s2 || bondq3[i].Qp==s2 )
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
        }
        
        if (ss1==0){
            if (act[s1])
                livpoint=livpoint-1;
            act[s1]=0;
        }
        
        if (ss2==0) {
            if (act[s2])
                livpoint=livpoint-1;
            act[s2]=0;
        }

       
        //add cut Q1
        for (int i=0; i<L; i++) {
           
            if (((bondj[i].Jleft==qleft3)&&(bondj[i].Jright==qleft2))&&((act[qleft3]==act[qleft2])&&(act[qleft2]==1))){
                
                if(leftq3c!=Bondmin){
                    leftq3c=leftq3c-log(16);
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq3c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq3c);
                }
                
                
                if(leftq2c!=Bondmin){
                    leftq2c=leftq2c-log(4);

                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq2c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq2c);
                }
            }
            
            
            if (((bondj[i].Jleft==qleft2)&&(bondj[i].Jright==qleft1))&&((act[qleft2]==act[qleft1])&&(act[qleft1]==1))){
                
                if(leftq3b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq3b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq3b);
                }
                
                if(leftq2b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq2b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq2b);
                }
            }
           
            if (((bondj[i].Jleft==qright2)&&(bondj[i].Jright==qright3))&&((act[qright2]==act[qright3])&&(act[qright3]==1))){
                
                if(rightq3c!=Bondmin){
                    rightq3c=rightq3c-log(16);

                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq3c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq3c);
                }
                
                if(rightq2c!=Bondmin){
                    rightq2c=rightq2c-log(4);

                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq2c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq2c);
                }
            }
            
            if (((bondj[i].Jleft==qright1)&&(bondj[i].Jright==qright2))&&((act[qright1]==act[qright2])&&(act[qright2]==1))){
                
                if(rightq3b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq3b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq3b);
                }
                
                if(rightq2b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq2b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq2b);
                }
            }
        }
        
        
        
        
        
        int last=-1;
        int rr=-1;
    
        if(livpoint==2){
        
            for (int i=0;i<L;i++){
                if (bondj[i].jj!=Jmin)
                    last=i;
            }
            
            //printf("last=%d\n",last);

        
            if(last!=-1){
                rr=abs(bondj[last].Jright-bondj[last].Jleft);
            
                //prepare dr
                if (rr==1){
                    int small = (int)compare(bondj[last].Jright,bondj[last].Jleft)-1;
                    singlet[small]=1;// small is a r=1 singlet index ,dimer[small]=1 ,else =-1
                }
                else if(rr==(L-1)) //(0,l-1) or (l-1,0)
                    singlet[L-1]=1;
            
                pair[bondj[last].Jright]=counter;
                pair[bondj[last].Jleft]=counter;
                
                showpair[bondj[last].Jright]=bondj[last].Jleft;
                showpair[bondj[last].Jleft]=bondj[last].Jright;
            
                counter=counter+1;
            }
        }
  /*
        //kill bond when (a,b)&(b,a)occurs
        for (int i=0; i<L; i++) {
            if((bondj[qleft].Jleft==bondj[i].Jright)&&(bondj[qleft].Jright==bondj[i].Jleft)){
                bondj[i].jj=Jmin;
                bondj[i].Jleft=-1;
                bondj[i].Jright=-1;
            }
        }
   */
    }
}

void rgq3(){
    
    double finalQ=0.0;
    int aa=bondindex;
    
    int kk = bondq3[aa].Qk;
    int ll = bondq3[aa].Ql;
    int mm = bondq3[aa].Qm;
    int nn = bondq3[aa].Qn;
    int oo = bondq3[aa].Qo;
    int pp = bondq3[aa].Qp;
    double qqq = bondq3[aa].qq;
    
    
    act[kk]=0; act[ll]=0; act[mm]=0; act[nn]=0; act[oo]=0; act[pp]=0;

    //load dimer
    singlet[kk]=1; singlet[mm]=1; singlet[oo]=1;

    //save cr
    pair[kk]=counter; pair[ll]=counter; counter=counter+1;
    pair[mm]=counter; pair[nn]=counter; counter=counter+1;
    pair[oo]=counter; pair[pp]=counter; counter=counter+1;
    
    showpair[kk]=ll;
    showpair[ll]=kk;
    showpair[mm]=nn;
    showpair[nn]=mm;
    showpair[oo]=pp;
    showpair[pp]=oo;


    if ( ((livpoint==6)&&(type==0)) ) {
        for (int i=0; i<L; i++) {
            bondj[i].Jleft=-1;
            bondj[i].Jright=-1;
            bondj[i].jj=Jmin;
        }
        localQ=qqq;
        livpoint=0;
    }
    
    else{
        
        double jj1=Bondmin;
        double jj2=Bondmin;
        
        double leftq3a=Bondmin;
        double leftq3b=Bondmin;
        double leftq3c=Bondmin;
        double leftq3d=Bondmin;
        double leftq3e=Bondmin;
        
        double rightq3a=Bondmin;
        double rightq3b=Bondmin;
        double rightq3c=Bondmin;
        double rightq3d=Bondmin;
        double rightq3e=Bondmin;

        //save neighbor(point)
        int qleft1=neighbor[kk].left1;
        int qleft2=neighbor[kk].left2;
        int qleft3=neighbor[kk].left3;
        int qleft4=neighbor[kk].left4;
        int qleft5=neighbor[kk].left5;
        
        int qright1=neighbor[pp].right1;
        int qright2=neighbor[pp].right2;
        int qright3=neighbor[pp].right3;
        int qright4=neighbor[pp].right4;
        int qright5=neighbor[pp].right5;
        
        /*
         -------------------------------------------------------------------------
                        kk      ll      mm      nn      oo      pp
         =======jj1=======                                      =======jj2========
         -------------------------------------------------------------------------
        */
        //choosed Q's neighbors
        int qleft=-1;
        int qright=-1;
        
        for(int i=0;i<L;i++){
            if((bondj[i].Jleft ==kk||bondj[i].Jleft ==ll||bondj[i].Jleft ==mm||bondj[i].Jleft ==nn||bondj[i].Jleft ==oo||bondj[i].Jleft ==pp||
                bondj[i].Jright==kk||bondj[i].Jright==ll||bondj[i].Jright==mm||bondj[i].Jright==nn||bondj[i].Jleft ==oo||bondj[i].Jleft ==pp))
            {
                if (bondj[i].Jright==kk){
                    qleft=i;
                    jj1=bondj[i].jj;
                }
                
                if (bondj[i].Jleft==pp){
                    qright=i;
                    jj2=bondj[i].jj;
                }
                bondj[i].jj=Jmin;
            }
        }
        
        //printf("qleft=%d,qright=%d\n",qleft,qright);
        //printf("jj1=%lf,jj2=%lf\n",jj1,jj2);
        
        /*
         --------------------------------------------------------------------------------------------------------------------------------------------------------
         left5     left4       left3       left2       left1       kk     ll     mm     nn    oo    pp     right1      right2      right3      right4      right5
                                                                   ===================================
         =========================leftq3e============================                               ==========================rightq3e==========================
                    ========================leftq3d=========================                  ========================rightq3d=======================
                                ===============leftq3c=============================     ========================rightq3c=================
                                            ===============leftq3b=========================
                                                                                  =====================rightq3b===============
                                                        ==================leftq3a================
                                                                          ==================rightq3a==============
         --------------------------------------------------------------------------------------------------------------------------------------------------------
        */
        //kill q3bond and save neighbor
        for (int i=0;i<L;i++) {
            if((bondq3[i].Qk==kk||bondq3[i].Qk==ll||bondq3[i].Qk==mm||bondq3[i].Qk==nn||bondq3[i].Qk==oo||bondq3[i].Qk==pp||
                bondq3[i].Ql==kk||bondq3[i].Ql==ll||bondq3[i].Ql==mm||bondq3[i].Ql==nn||bondq3[i].Ql==oo||bondq3[i].Ql==pp||
                bondq3[i].Qm==kk||bondq3[i].Qm==ll||bondq3[i].Qm==mm||bondq3[i].Qm==nn||bondq3[i].Qm==oo||bondq3[i].Qm==pp||
                bondq3[i].Qn==kk||bondq3[i].Qn==ll||bondq3[i].Qn==mm||bondq3[i].Qn==nn||bondq3[i].Qn==oo||bondq3[i].Qn==pp||
                bondq3[i].Qo==kk||bondq3[i].Qo==ll||bondq3[i].Qo==mm||bondq3[i].Qo==nn||bondq3[i].Qo==oo||bondq3[i].Qo==pp||
                bondq3[i].Qp==kk||bondq3[i].Qp==ll||bondq3[i].Qp==mm||bondq3[i].Qp==nn||bondq3[i].Qp==oo||bondq3[i].Qp==pp
                )&&
               (bondq3[i].qq!=Jmin))
            {
                
                    if ((bondq3[i].Qk==qleft1)&&(bondq3[i].Ql==kk))
                        leftq3a=bondq3[i].qq; //per
                    
                    if ((bondq3[i].Qk==qleft2)&&(bondq3[i].Ql==qleft1))
                        leftq3b=bondq3[i].qq; //cutQ1
                    
                    if ((bondq3[i].Qk==qleft3)&&(bondq3[i].Ql==qleft2))
                        leftq3c=bondq3[i].qq; //per & cutQ1
                    
                    if ( (bondq3[i].Qk==qleft4)&&(bondq3[i].Ql==qleft3)&&(bondq3[i].Qm==qleft2)&&(bondq3[i].Qn==qleft1) ){ //cutQ2
                        if (bondq2[i].qq==Jmin)
                            bondq2[i].qq=bondq3[i].qq;
                        else
                            bondq2[i].qq=compare(bondq2[i].qq,bondq3[i].qq);
                        
                        leftq3d=bondq3[i].qq;
                    }
                
                    if ( (bondq3[i].Qk==qleft5)&&(bondq3[i].Ql==qleft4)&&(bondq3[i].Qm==qleft3)&&(bondq3[i].Qn==qleft2) ){  //per & cutQ2
                        leftq3e=bondq3[i].qq;

                        if (bondq2[i].qq==Jmin)
                            bondq2[i].qq=bondq3[i].qq-log(4);
                        else
                            bondq2[i].qq=compare(bondq2[i].qq,bondq3[i].qq-log(4));
                        
                    }
                
                    if ((bondq3[i].Qk==ll)&&(bondq3[i].Ql==mm))
                        rightq3a=bondq3[i].qq; //per
                    
                    if ((bondq3[i].Qk==mm)&&(bondq3[i].Ql==nn))
                        rightq3b=bondq3[i].qq; //cutQ1
                    
                    if ((bondq3[i].Qk==nn)&&(bondq3[i].Ql==oo))
                        rightq3c=bondq3[i].qq; //per & cutQ1
                    
                    if ( (bondq3[i].Qm==qright1)&&(bondq3[i].Qn==qright2)&&(bondq3[i].Qo==qright3)&&(bondq3[i].Qp==qright4) ){ //cutQ2
                        if (bondq2[neighbor[i].right2].qq==Jmin)
                            bondq2[neighbor[i].right2].qq=bondq3[i].qq;
                        else
                            bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq,bondq3[i].qq);

                        rightq3d=bondq3[i].qq; //cut Q2
                    }
                
                    if ( (bondq3[i].Qm==qright2)&&(bondq3[i].Qn==qright3)&&(bondq3[i].Qo==qright4)&&(bondq3[i].Qp==qright5) ){  //per & cutQ2
                        rightq3e=bondq3[i].qq;

                        if (bondq2[neighbor[i].right2].qq==Jmin)
                            bondq2[neighbor[i].right2].qq=bondq3[i].qq-log(4);
                        else
                            bondq2[neighbor[i].right2].qq=compare(bondq2[neighbor[i].right2].qq,bondq3[i].qq-log(4));
                        
                        rightq3e=bondq3[i].qq;
                    }
                
                bondq3[i].Qk=-1;
                bondq3[i].Ql=-1;
                bondq3[i].Qm=-1;
                bondq3[i].Qn=-1;
                bondq3[i].Qo=-1;
                bondq3[i].Qp=-1;
                bondq3[i].qq=Jmin;
            }
        }
        //printf("left1q1=%lf,left1q2=%lf,left1q3=%lf,left2q=%lf,left3q=%lf,right1q1=%lf,right1q2=%lf,right1q3=%lf,rightq2=%lf,righq3=%lf\n",left1q1,left1q2,left1q3,left2q,left3q,right1q1,right1q2,right1q3,right2q,right3q);
        
         /*
         ------------------------------------------------------------------------------------------------------------------
         left3      left2       left1       kk      ll      mm      nn      oo      pp      right1      right2      right3
                                            ==========================================
         ================leftq2c==============                                      ================rightq2c==============
                    ===============leftq2b============                      =============rightq2b=============
                                ==========leftq2a=============
                                                                    =================rightq2a=====
         ------------------------------------------------------------------------------------------------------------------
         */
        
        //kill q2bond and save neighbor
        double leftq2a=Bondmin;
        double leftq2b=Bondmin;
        double leftq2c=Bondmin;
        double rightq2a=Bondmin;
        double rightq2b=Bondmin;
        double rightq2c=Bondmin;

        for (int i=0;i<L;i++) {
            if((bondq2[i].Q2k==kk||bondq2[i].Q2k==ll||bondq2[i].Q2k==mm||bondq2[i].Q2k==nn||bondq2[i].Q2k==oo||bondq2[i].Q2k==pp||
                bondq2[i].Q2l==kk||bondq2[i].Q2l==ll||bondq2[i].Q2l==mm||bondq2[i].Q2l==nn||bondq2[i].Q2l==oo||bondq2[i].Q2l==pp||
                bondq2[i].Q2m==kk||bondq2[i].Q2m==ll||bondq2[i].Q2m==mm||bondq2[i].Q2m==nn||bondq2[i].Q2m==oo||bondq2[i].Q2m==pp||
                bondq2[i].Q2n==kk||bondq2[i].Q2n==ll||bondq2[i].Q2n==mm||bondq2[i].Q2n==nn||bondq2[i].Q2n==oo||bondq2[i].Q2n==pp)&&
               (bondq2[i].qq!=Jmin))
            {
                if((bondq2[i].Q2m==qleft1)&&(bondq2[i].Q2n==kk)){
                    leftq2c=bondq2[i].qq;//cut&per
                }
                
                if((bondq2[i].Q2m==kk)&&(bondq2[i].Q2n==ll)){
                    leftq2b=bondq2[i].qq;//cut
                }
                
                if((bondq2[i].Q2m==ll)&&(bondq2[i].Q2n==mm)){
                    leftq2a=bondq2[i].qq;//per
                }
                
                if((bondq2[i].Q2k==pp)&&(bondq2[i].Q2l==qright1)){
                    rightq2c=bondq2[i].qq;//cut&per
                }
                
                if((bondq2[i].Q2k==oo)&&(bondq2[i].Q2l==pp)){
                    rightq2b=bondq2[i].qq;//cut
                }
                
                if((bondq2[i].Q2k==nn)&&(bondq2[i].Q2l==oo)){
                    rightq2a=bondq2[i].qq;//per
                }
                
                bondq2[i].Q2k=-1;
                bondq2[i].Q2l=-1;
                bondq2[i].Q2m=-1;
                bondq2[i].Q2n=-1;
                bondq2[i].qq=Jmin;
            }
        }
        
        //printf("q1qq=%lf,q2qq=%lf,q3qq=%lf,q4qq=%lf,qqleft1=%lf,qqleft2=%lf,qqright1=%lf,qqright2=%lf\n",q1qq,q2qq,q3qq,q4qq,qqleft1,qqleft2,qqright1,qqright2);
        
        double effqq1=Bondmin;
        double effqq2=Bondmin;
        double effJ=Bondmin;
        
        load1[0]=jj1;
        load1[1]=leftq3e;
        load1[2]=leftq2c;
        load1[3]=leftq2a;
        load1[4]=leftq3c;
        load1[5]=leftq3a;

        load2[0]=jj2;
        load2[1]=rightq3e;
        load2[2]=rightq2c;
        load2[3]=rightq2a;
        load2[4]=rightq3c;
        load2[5]=rightq3a;

        //comparison of which combination does the most good
        for (int i=0; i<6; i++) {
            
            if (load1[i]!=Jmin) {
                
                if(i<3){
                    for (int j=5; j<6; j++){
                        if (load2[j]!=Jmin)
                            effJ=compare(effJ,load1[i]+load2[j]-log(32)-qqq);
                    }
                }
            
                else if ((i>=3)&&(i<=4)) {
                    for (int j=3; j<6; j++){
                        if (load2[j]!=Jmin)
                            effJ=compare(effJ,load1[i]+load2[j]-log(32)-qqq);
                    }
                }
            
                else{
                    for (int j=0; j<6; j++){
                        if (load2[j]!=Jmin)
                            effJ=compare(effJ,load1[i]+load2[j]-log(32)-qqq);
                    }
                }
            }
        }
       
        //create jbond for pbc
        bondj[qleft].jj=effJ;
        
        if (effJ<=Bondmin)//if no effJ
            bondj[qleft].jj=Bondmin;;
        
        bondj[qleft].Jleft=bondj[qleft].Jleft;
        int s1=bondj[qleft].Jleft;
        bondj[qleft].Jright=bondj[qright].Jright;
        int s2=bondj[qright].Jright;
        
        bondj[qright].Jleft=-1; bondj[qright].Jright=-1; bondj[qright].jj=Jmin;
        
        livpoint=livpoint-6;
        
        int ss1=0;
        int ss2=0;
        
        //Check s1&s2
        for (int i=0; i<L;i++) {
            if (bondj[i].jj!=Jmin) {
                if(bondj[i].Jright==s1 || bondj[i].Jleft==s1)
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if(bondj[i].Jright==s2 || bondj[i].Jleft==s2)
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
            
            if (bondq2[i].qq!=Jmin) {
                if (bondq2[i].Q2k==s1 || bondq2[i].Q2l==s1 ||bondq2[i].Q2m==s1 ||bondq2[i].Q2n==s1 )
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if (bondq2[i].Q2k==s2 || bondq2[i].Q2l==s2 ||bondq2[i].Q2m==s2 ||bondq2[i].Q2n==s2 )
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
            
            if (bondq3[i].qq!=Jmin) {
                if (bondq3[i].Qk==s1 || bondq3[i].Ql==s1 ||bondq3[i].Qm==s1 || bondq3[i].Qn==s1 || bondq3[i].Qo==s1 || bondq3[i].Qp==s1 )
                    ss1=ss1+1;
                else
                    ss1=ss1;
                
                if (bondq3[i].Qk==s2 || bondq3[i].Ql==s2 ||bondq3[i].Qm==s2 || bondq3[i].Qn==s2 || bondq3[i].Qo==s2 || bondq3[i].Qp==s2 )
                    ss2=ss2+1;
                else
                    ss2=ss2;
            }
        }
        
        if (ss1==0){
            if (act[s1])
                livpoint=livpoint-1;
            act[s1]=0;
        }
        
        if (ss2==0) {
            if (act[s2])
                livpoint=livpoint-1;
            act[s2]=0;
        }
        
        //add cut Q1
        for (int i=0; i<L; i++) {
          
            if (((bondj[i].Jleft==qleft3)&&(bondj[i].Jright==qleft2))&&((act[qleft3]==act[qleft2])&&(act[qleft2]==1))){
                
                if(leftq3c!=Bondmin){
                    leftq3c=leftq3c-log(16);
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq3c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq3c);
                }
                
                if(leftq2c!=Bondmin){
                    leftq2c=leftq2c-log(4);
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq2c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq2c);
                }
            }
           
            if (((bondj[i].Jleft==qleft2)&&(bondj[i].Jright==qleft1))&&((act[qleft2]==act[qleft1])&&(act[qleft1]==1))){
                
                if(leftq3b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq3b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq3b);
                }
                
                if(leftq2b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=leftq2b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,leftq2b);
                }
            }
           
            if (((bondj[i].Jleft==qright2)&&(bondj[i].Jright==qright3))&&((act[qright2]==act[qright3])&&(act[qright3]==1))){
                
                if(rightq3c!=Bondmin){
                    rightq3c=rightq3c-log(16);
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq3c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq3c);
                }
                
                if(rightq2c!=Bondmin){
                    rightq2c=rightq2c-log(4);
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq2c;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq2c);
                }
            }
            
            if (((bondj[i].Jleft==qright1)&&(bondj[i].Jright==qright2))&&((act[qright1]==act[qright2])&&(act[qright2]==1))){
                
                if(rightq3b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq3b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq3b);
                }
                
                if(rightq2b!=Bondmin){
                    if (bondj[i].jj==Jmin)
                        bondj[i].jj=rightq2b;
                    else
                        bondj[i].jj=compare(bondj[i].jj,rightq2b);
                }
            }
        }

        
        int last=-1;
        int rr=-1;
        
        if(livpoint==2){
            for (int i=0;i<L;i++){
                if (bondj[i].jj!=Jmin)
                    last=i;
                }
            //printf("last=%d\n",last);
            
            if(last!=-1){
                rr=abs(bondj[last].Jright-bondj[last].Jleft);
                
                //prepare dr
                if (rr==1){
                    int small = (int)compare(bondj[last].Jright,bondj[last].Jleft)-1;
                    singlet[small]=1;// small is a r=1 singlet index ,dimer[small]=1 ,else =-1
                    }
        
            else if(rr==(L-1)) //(0,l-1) or (l-1,0)
                singlet[L-1]=1;
                
                pair[bondj[last].Jright]=counter;
                pair[bondj[last].Jleft]=counter;
            
                showpair[bondj[last].Jright]=bondj[last].Jleft;
                showpair[bondj[last].Jleft]=bondj[last].Jright;
                
                counter=counter+1;
            }
        }
      /*
        //kill bond when (a,b)&(b,a)occurs
        for (int i=0; i<L; i++) {
            if((bondj[qleft].Jleft==bondj[i].Jright)&&(bondj[qleft].Jright==bondj[i].Jleft)){
                bondj[i].jj=Jmin;
                bondj[i].Jleft=-1;
                bondj[i].Jright=-1;
            }
        }
       */
    }
}

void Dimer(){
    
    int j,k,m;   //i j k m means point and r=k-i;
    
    for (int i=0; i<L ; i++){
        for (int r=2; r<=L/2; r++) {
            k=(i+r)%L;
            if(singlet[i]&&singlet[k]){ //(i,j) (k,m)
                dr[r]=dr[r]+0.5625;
                drcopy[r]=drcopy[r]+0.5625;
            }
        }
    }
}


void Cc(){
    
    for (int i=0; i<L/2 ; i++){
        for (int r=1; r<=L/2; r++) {
            int x = (i+r)%L;
            if((pair[i]==pair[x])&&(pair[i]!=-1)){ //(i,j) (k,m)
                cr[r]=cr[r]+0.75;
		//printf("(%d,%d)\n",i,x);
            }
        }
    }
}

void transposeloop(){
    
    int i=0;
    int j=0;
    int lll=0;
    int k=0; //number used to label the loop
    
    int lnum=0;
    
    
    
    for (int z=0; z<L/2; z++){
        siteA[z]=0;
        loopsize[z]=0;
    }
    
    for (int z=0; z<L; z++) {
        Lambda[z]=0;
        Oo[z]=0;
        Vr[z]=showpair[z];
        Vl[z]=showpair[z];
    }
    
    
    
    //counting sublattice A;just for 1D case
    
    for (int z=1; z<=L/2; z++)  siteA[z]=2*z-1;
    
    //printf("yes\n");
    
    //for (int z=0; z<L; z++) printf("loop%d,%d\n",z,showpair[z]);
    
    
   // for (int z=1; z<=L/2; z++)  printf("siteA%d,%d\n",z,siteA[z]);
    
    
    
    for (int lll=1; lll<=L/2; lll++) {
        j=siteA[lll];
        if (Lambda[j]!=0)
            continue;
        
        lnum=lnum+1;
        i=j;
        // printf("yes\n");
        
        while(1){
            loopsize[lnum]=loopsize[lnum]+2;
            Lambda[i]=lnum; Oo[i]=k; i=Vl[i]; k=k+1;
            Lambda[i]=lnum; Oo[i]=k; i=Vr[i]; k=k+1;
            
            //printf("(i,j)=(%d,%d)\n",i,j);
            
            if (i==j)
                break;
        }
        
    }
    
}

void Dimer2(){
    
    int j,k,m,Sk,Sl,Op,Ol;   //i j k m means point and r=k-i;
    
    for (int i=0; i<L ; i++){
        j=(i+1)%L;
        
        if (Lambda[i]==Lambda[j]) {
            Ol=Oo[i]+Oo[j]-compare(Oo[i],Oo[j]);
            Op=compare(Oo[i],Oo[j]);
            
            Ol=2*(Ol/2)+1;
            Op=2*(Op/2);
            
        }
        for (int r=2; r<=L/2; r++) {
            k=(i+r)%L;
            m=(k+1)%L;
            
            if(Lambda[i]==Lambda[j]){ //(i,j)
                if(Lambda[k]==Lambda[m]){ //(k,m)
                    if(Lambda[k]==Lambda[i]){ //(i,j,k,m)
                        
                        if ((Oo[k]>=Ol)&&(Oo[k]<=Op))
                            Sk=1;
                        else
                            Sk=2;
                        
                        if ((Oo[m]>=Ol)&&(Oo[m]<=Op))
                            Sl=1;
                        else
                            Sl=2;
                        
                        
                        if (Sk==Sl){
                            dr2[r]=dr2[r]+0.5625;
                            drcopy2[r]=drcopy2[r]+0.5625;
                            dr3[r]=dr3[r]+0.5625;
                            drcopy3[r]=drcopy3[r]+0.5625;

                            //printf("%d\t\t\tothers\n",r);
                        }
                        else{
                            dr2[r]=dr2[r]-0.1875;
                            drcopy2[r]=drcopy2[r]-0.1875;
                            //printf("%d\t\t\tothers\n",r);
                            dr3[r]=dr3[r]-0.1875;
                            drcopy3[r]=drcopy3[r]-0.1875;

                            
                        }
                    }
                    
                    else {//(i,j) (k,m)
                        dr[r]=dr[r]+0.5625;
                        drcopy[r]=drcopy[r]+0.5625;
                        
                        dr3[r]=dr3[r]+0.5625;
                        drcopy3[r]=drcopy3[r]+0.5625;

                        //printf("%d\tcase1\n",r);
                    }
                    
                    
                }
            }
            
            
            else if ((Lambda[i]==Lambda[k])&&(Lambda[j]==Lambda[m])){ //(i,k) (j,m)
                dr2[r]=dr2[r]+0.1875;
                drcopy2[r]=drcopy2[r]+0.1875;
                dr3[r]=dr3[r]+0.1875;
                drcopy3[r]=drcopy3[r]+0.1875;

                //printf("%d\t\t\tothers\n",r);
            }
            
            
            else if ((Lambda[i]==Lambda[m])&&(Lambda[j]==Lambda[k])){ //(i,m) (j,k)
                dr2[r]=dr2[r]+0.1875;
                drcopy2[r]=drcopy2[r]+0.1875;
                dr3[r]=dr3[r]+0.1875;
                drcopy3[r]=drcopy3[r]+0.1875;

                //printf("%d\t%d\t\tcase2\n",i,r);

            }
            
        }
    }
    
}


void Dimerstar(){
    
    double de=0;
    
    de = 0.5*(drcopy[0]-drcopy[1]);
    drsta[0] = drsta[0] + de;
    
    de = 0.5*(drcopy2[0]-drcopy2[1]);
    drsta2[0] = drsta2[0] + de;
    
    de = 0.5*(drcopy3[0]-drcopy3[1]);
    drsta3[0] = drsta3[0] + de;

    
    for (int i=2; i<L/2; i++) {
        de = 0.5*(drcopy[i]-0.5*drcopy[i-1]-0.5*drcopy[i+1])*pow(-1,i);
        drsta[i] = drsta[i] + de;
        
        de = 0.5*(drcopy2[i]-0.5*drcopy2[i-1]-0.5*drcopy2[i+1])*pow(-1,i);
        drsta2[i] = drsta2[i] + de;
        
        de = 0.5*(drcopy3[i]-0.5*drcopy3[i-1]-0.5*drcopy3[i+1])*pow(-1,i);
        drsta3[i] = drsta3[i] + de;
    }
    
    de = 0.5*(drcopy[L/2]-drcopy[(L/2)-1])*pow(-1,L/2);
    drsta[L/2] = drsta[L/2] + de;
    
    de = 0.5*(drcopy2[L/2]-drcopy2[(L/2)-1])*pow(-1,L/2);
    drsta2[L/2] = drsta2[L/2] + de;
    
    de = 0.5*(drcopy3[L/2]-drcopy3[(L/2)-1])*pow(-1,L/2);
    drsta3[L/2] = drsta3[L/2] + de;
    
}



void writesinglet(){
    
    for(int i=0;i<L;i++)
        printf("singlet%d=%d\n",i,singlet[i]);
}

void writepair(){
    
    for(int i=0; i<L; i++)
        if((showpair[i]>i)&&(showpair[i]>0))
            printf("%d\t%d\n",i,showpair[i]);
}

void WriteYR(){
    for(int i=0;i<L;i++)
        printf("%d\t%lf\t%d\n",i,Q3[i],showpair[i]);
}


void writefile(){
    
    if (localQ==0){
        
        FILE *fp=fopen(name, "a");
        if(fp==NULL){
            printf("XX\n");
        }
        else{
                fprintf(fp,"%lf\n",localJ);
                fclose(fp);
            }
    }
    
    else if(localJ==0){
        
        FILE *fp=fopen(name2, "a");
        if(fp==NULL){
            printf("XX\n");
        }
        else{
                fprintf(fp,"%lf\n",localQ);
                fclose(fp);
            }
    }

}

void writeCrr(){
        
        FILE *fpp=fopen(name3, "a");
        if(fpp==NULL)
            printf("YY\n");
            //error++;
        
        else{
            for (int i=1; i<L;i++ ) {
                if(cr[i]>0)
                    fprintf(fpp,"%d\t%.16f\n",i,cr[i]/(L/2*N));
            }
            fclose(fpp);
        }
    /*
        FILE *ffffp=fopen(name7, "a");
        if(ffffp==NULL)
            printf("XX\n");
       
        else{
            
            for (int i=1; i<=L/2; i++ ) {
                if(lncr[i]!=0)
                    fprintf(ffffp,"%d %.10f\n",i,lncr[i]);
            }
            fclose(ffffp);
        }
	*/
}


void writeDrr(){
    
        FILE *fp=fopen(name5, "a");
        if(fp==NULL)
            printf("XX\n");
        
        else{
            for (int i=2; i<=L/2; i++ )
                fprintf(fp,"%d\t%.10f\t%.10f\t%.10f\n",i,dr3[i]/(L*N),dr[i]/(L*N),dr2[i]/(L*N));
             fclose(fp);
        }
    
    
        FILE *ffp=fopen(name4, "a");
        if(ffp==NULL)
            printf("XX\n");
        
        else{
            for (int i=2; i<=L/2; i++ )
                fprintf(ffp,"%d\t%.10f\t%.10f\t%.10f\n",i,drsta3[i]/(L*N),drsta[i]/(L*N),drsta2[i]/(L*N));
            fclose(ffp);
        }
    
   /*
        FILE *f4p=fopen(name9, "a");
        if(f4p==NULL)
            printf("XX\n");
        
        else{
            for (int i=2; i<=L/2; i++ )
                fprintf(f4p,"%d\t%.10f\n",i,dr2[i]/(L*N));
            fclose(f4p);
        }
    
        FILE *f2p=fopen(name10, "a");
        if(f2p==NULL)
            printf("XX\n");
    
        else{
            for (int i=2; i<=L/2; i++ )
                fprintf(f2p,"%d\t%.10f\n",i,drsta2[i]/(L*N));
            fclose(f2p);
        }
   */
}
/*
void Writenum(){
 
    for (int i=1; i<=L/2; i++ ) {
        
        FILE *fffffpp=fopen(name9, "a");
        if(fffffpp==NULL)
            printf("XX\n");
        
        else{
            for (int j=0; j<Nc[i]; j++)
                fprintf(fffffpp,"%d\n",i);
            fclose(fffffpp);
        }
    }
    
    for (int i=1; i<=L/2; i++ ) {
        
        FILE *fffffppp=fopen(name10, "a");
        if(fffffppp==NULL)
            printf("XX\n");
        
        else{
            for (int j=0; j<Nd[i]; j++)
                fprintf(fffffppp,"%d\n",i);
            fclose(fffffppp);
        }
    }
}
 */

void WriteE(){
    
    FILE *fffp=fopen(name6, "a");
    if(fffp==NULL)
        printf("XX\n");
    
    else{
        fprintf(fffp,"%d\t%.10f\t%.10f\n",livpoint,Esum,Emax);
        fclose(fffp);
    }
}

void writedata(){
    
    printf("===========Q1Bond===========\n");
    
    for(int i=0;i<L;i++){
        if(bondj[i].jj!=Jmin)
            printf("%d Jbond:(%d,%d) J:%lf\n",i,bondj[i].Jleft,bondj[i].Jright,bondj[i].jj);
    }
    printf("===========Q2Bond===========\n");
    for(int i=0;i<L;i++){
        if(bondq2[i].qq!=Jmin)
            printf("%d Qbond:(%d,%d,%d,%d) Q:%lf\n",i,bondq2[i].Q2k,bondq2[i].Q2l,bondq2[i].Q2m,bondq2[i].Q2n,bondq2[i].qq);
    }
    
    printf("===========Q3Bond===========\n");
    for(int i=0;i<L;i++){
        if(bondq3[i].qq!=Jmin)
            printf("%d Qbond:(%d,%d,%d,%d,%d,%d) Q:%lf\n",i,bondq3[i].Qk,bondq3[i].Ql,bondq3[i].Qm,bondq3[i].Qn,bondq3[i].Qo,bondq3[i].Qp,bondq3[i].qq);
    }
    printf("============================\n");
    
    for (int i=0; i<L; i++) {
        if(act[i]!=0)
            printf("remain:%d\n",i);
    }
    
    printf("actpoint:%d\n",livpoint);
}



int main(int argc, char** argv){
    
    long int iseed = -1;
    //int  myid, nprocs;
    //MPI_Init(&argc,&argv);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    //MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    for(int i = 1; i<argc; i++ ){
        if (strstr(argv[i], "-L" )) L        = atoi(argv[i]+2);
        if (strstr(argv[i], "-N" )) N        = atoi(argv[i]+2);
    	if (strstr(argv[i], "-D" )) D        = atof(argv[i]+2);
        if (strstr(argv[i], "-I" )) iseed    = atoi(argv[i]+2);
        if (strstr(argv[i], "-C" )) Check    = atoi(argv[i]+2);
        if (strstr(argv[i], "-W" )) Write    = atoi(argv[i]+2);
        if (strstr(argv[i], "-G" )) Gap      = atoi(argv[i]+2);
        if (strstr(argv[i], "-B" )) Bins     = atoi(argv[i]+2);
        if (strstr(argv[i], "-E" )) EE       = atoi(argv[i]+2);
        if (strstr(argv[i], "-M" )) Num      = atoi(argv[i]+2);
        if (strstr(argv[i], "-K" )) CCheck   = atoi(argv[i]+2);
    }
/*
    sprintf(name,"localJl%dD%lg_%d",L,D,myid);
    sprintf(name2,"localQl%dD%lg_%d",L,D,myid);
    sprintf(name3,"crl%dD%lg_%d",L,D,myid);
    sprintf(name4,"drstal%dD%lg_%d",L,D,myid);
    sprintf(name5,"drl%dD%lg_%d",L,D,myid);
    sprintf(name6,"EEl%dD%lg_%d",L,D,myid);
    sprintf(name7,"lncrl%dD%lg_%d",L,D,myid);
    sprintf(name8,"lndrstal%dD%lg_%d",L,D,myid);
  */
    sprintf(name,"localJl%dD%lg",L,D);
    sprintf(name2,"localQl%dD%lg",L,D);
    sprintf(name3,"crl%dD%lg",L,D);
    sprintf(name4,"drstal%dD%lg",L,D);
    sprintf(name5,"drl%dD%lg",L,D);
    sprintf(name6,"EEl%dD%lg",L,D);
    sprintf(name7,"lncrl%dD%lg",L,D);
    sprintf(name8,"lndrstal%dD%lg",L,D);

    allocateArray();
    time_t now;
    time(&now);
    if(iseed<0) time ((time_t*)&iseed);    //     iseed>0 -> no random
    r250_srandom(iseed);        // assign seed to random number generator

    
    for (int kk=0; kk<Bins; kk++) {
        resetArray();
        int trigger=0;
        
        for (int z=0; z<N; z++) {
            resetVariable();
            initialize();
            if (Check) writedata();
        
            while(true){
                
                Emax=max();
               // if (Check) printf("%lf(%d),%d\n",max(),type,bondindex);
                
                
                if (EE) {
                    AddE();
                    WriteE();
                }

                if (type==0)     //rg when jmax
                    rgq1();
                
                else if(type==1)
                    rgq2();

                else if(type ==2)             //rg when qmax
                    rgq3();
                
                //if (Check) writedata();
                
                maxx=max();
               
                if(maxx!=Bondmin && maxx!=Jmin){
                    localJ=maxx;
                    maxxtype=type;
                }
                
                //printf("localJ=%lf,localQ=%lf\n",localJ,localQ);

                
                //if(maxx==Bondmin||livpoint<2){
                    if(livpoint<2){
    
                    //printf("localJ=%lf,localQ=%lf\n",localJ,localQ);
                    break;
                    }
                
                if (livpoint<0) {
                    printf("localJ=%lf,localQ=%lf\n",localJ,localQ);
                    printf("iseed=%d\n",iseed);
                    writedata();
                    trigger=1;
                    exit(1);
                    }
            }
            if(CCheck){
               // writepair();
                WriteYR();
               // writesinglet();
            }
            
            transposeloop();
            //Dimer();
            Dimer2();
            Dimerstar();
            Cc();
            
            if(Gap){
               // if (maxxtype==0)
                writefile();
            }
        }
        
        //printf("HI\n");
        //if (Num)
            //Writenum();
            
        if (Write){
            //printf("HI\n");

            writeCrr();
            writeDrr();
        }
    }
    //MPI_Finalize();
    freeArray();
}

