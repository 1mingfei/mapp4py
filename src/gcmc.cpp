/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "elements.h"
#include <mpi.h>
#include "gcmc.h"
#include "memory.h"
#include "random.h"
#include "neighbor.h"
#include "ff_md.h"
#include "MAPP.h"
#include "atoms_md.h"
#include "comm.h"
#include "dynamic_md.h"
using namespace MAPP_NS;
#define GCMCDEBUG
/*--------------------------------------------
 constructor
 --------------------------------------------*/
GCMC::GCMC(AtomsMD*& __atoms,ForceFieldMD*&__ff,DynamicMD*& __dynamic,elem_type __gas_type0,type0 __mu0,elem_type __gas_type1,type0 __mu1,type0 __T,int seed):
atoms(__atoms),
ff(__ff),
world(__atoms->comm.world),
dynamic(__dynamic),
gas_type0(__gas_type0),
gas_type1(__gas_type1),
mu0(__mu0),
mu1(__mu1),
T(__T),
cut_sq(__ff->cut_sq),
s_hi(__atoms->comm.s_hi),
s_lo(__atoms->comm.s_lo),
natms_lcl(__atoms->natms_lcl),
natms_ph(__atoms->natms_ph),
ielem(gas_type0)
{
        ngas_lclArr = new int [nGasType];
        muArr = new type0 [nGasType];// nullptr; //-mingfei change lambda to a vector
        muArr[0] = mu0;
        muArr[1] = mu1;
        gas_typeArr = new elem_type [nGasType]; //-mingfei change lambda to a vector
        gas_typeArr[0] = gas_type0;
        gas_typeArr[1] = gas_type1;
        gas_massArr = new type0 [nGasType]; //-mingfei change lambda to a vector
        lambdaArr   = new type0 [nGasType]; //-mingfei change lambda to a vector
        sigmaArr    = new type0 [nGasType]; //-mingfei change sigma to a vector
        z_facArr    = new type0 [nGasType]; //-mingfei change z_facArr to a vector
        ngasArr     = new int [nGasType];
 

    random=new Random(seed);
    s_trials=new type0*[__dim__];
    *s_trials=NULL;
    del_ids=NULL;
    del_ids_sz=del_ids_cpcty=0;
    vars=lcl_vars=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GCMC::~GCMC()
{
    delete [] del_ids;
    delete [] s_trials;
    delete random;
    // -mingfei release simgaArr
    delete [] ngas_lclArr;
    delete [] muArr;
    delete [] gas_typeArr;
    delete [] gas_massArr;
    delete [] lambdaArr;
    delete [] sigmaArr;
    delete [] z_facArr;
    delete [] ngasArr;


}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::add_del_id(int* new_ids,int no)
{
    if(del_ids_sz+no>del_ids_cpcty)
    {
        int* del_ids_=new int[del_ids_sz+no];
        memcpy(del_ids_,del_ids,del_ids_sz*sizeof(int));
        del_ids_cpcty=del_ids_sz+no;
        delete [] del_ids;
        del_ids=del_ids_;
    }
    memcpy(del_ids+del_ids_sz,new_ids,sizeof(int)*no);
    del_ids_sz+=no;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int GCMC::get_new_id()
{
    if(del_ids_sz)
    {
        del_ids_sz--;
        return del_ids[del_ids_sz];
    }
    else
    {
        max_id++;
        return max_id;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::init()
{
    cut=ff->cut[ielem][0];
    for(size_t i=1;i<atoms->elements.nelems;i++)
        cut=MAX(cut,ff->cut[ielem][i]);
    
    gas_mass0=atoms->elements.masses[gas_type0];
    gas_mass1=atoms->elements.masses[gas_type1];
    //-mingfei
    gas_massArr[0]=atoms->elements.masses[gas_type0];
    gas_massArr[1]=atoms->elements.masses[gas_type1];
    //-mingfei 

    kbT=atoms->kB*T;
    beta=1.0/kbT;
     
    lambda=atoms->hP/sqrt(2.0*M_PI*kbT*gas_mass0);
    sigma=sqrt(kbT/gas_mass0);
    z_fac=1.0;
    for(int i=0;i<__dim__;i++) z_fac/=lambda;
    z_fac*=exp(beta*mu0);
    for (int i = 0;i < nGasType;i++)
    {
        lambdaArr[i] = atoms->hP/sqrt(2.0*M_PI*kbT*gas_massArr[i]);
        sigmaArr[i] = sqrt(kbT/gas_massArr[i]);
        z_facArr[i] = 1.0;
        for(int j=0;j<__dim__;j++) z_facArr[i] /= lambdaArr[i];
        z_facArr[i] *= exp(beta*muArr[i]);
        ngas_lclArr[i] = 0;
        ngasArr[i] = 0;
#ifdef GCMCDEBUG
        FILE* fp_debug=NULL;
        if(atoms->comm_rank==0)
        {
            fp_debug=fopen("gcmc_debug","a");
            fprintf(fp_debug,"i\t%d\n",i);
            fprintf(fp_debug,"sigma\t%e\n",sigmaArr[i]);
            fprintf(fp_debug,"gas_mass0\t%e\n",gas_mass0);
            fprintf(fp_debug,"gas_mass1\t%e\n",gas_mass1);
            fprintf(fp_debug,"gas_massArr\t%e\n",gas_massArr[i]);
            fprintf(fp_debug,"mu\t%e\n",muArr[i]);
        }
#endif

    }

    vol=1.0;
    for(int i=0;i<__dim__;i++)vol*=atoms->H[i][i];
    
    unsigned int max_id_=0;
    unsigned int* id=atoms->id->begin();
    for(int i=0;i<natms_lcl;i++)
        max_id_=MAX(id[i],max_id_);
    MPI_Allreduce(&max_id_,&max_id,1,MPI_UNSIGNED,MPI_MAX,world);
    for(int i=0;i<del_ids_sz;i++)
        max_id=MAX(max_id,del_ids[i]);
        
    ngas_lcl=0;
    elem_type* elem=atoms->elem->begin();
    for(int i=0;i<natms_lcl;i++) 
    {
        if(elem[i]==gas_type0) ngas_lcl++; 
        for (int j=0; j< nGasType;j++) // -mingfei 
        {
            if(elem[i]==gas_typeArr[j]) ngas_lclArr[j]++; //-mingfei
        }
    }
    MPI_Allreduce(&ngas_lcl,&ngas,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(ngas_lclArr,ngasArr,nGasType,MPI_INT,MPI_SUM,world);//-mingfei
#ifdef GCMCDEBUG
        FILE* fp_debug=NULL;
        if(atoms->comm_rank==0)
        {
            fp_debug=fopen("gcmc_debug","a");
            fprintf(fp_debug,"ngasArr\t%d\t%d\n",ngasArr[0],ngasArr[1]);
        }
#endif


}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::fin()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_setup()
{
    int sz=0;
    max_ntrial_atms=1;
    for(int i=0;i<__dim__;i++)
    {
        type0 tmp=0.0;
        for(int j=i;j<__dim__;j++)
            tmp+=atoms->B[j][i]*atoms->B[j][i];
        cut_s[i]=sqrt(tmp)*cut;
        

        
        s_lo_ph[i]=s_lo[i]-cut_s[i];
        s_hi_ph[i]=s_hi[i]+cut_s[i];
        nimages_per_dim[i][0]=static_cast<int>(floor(s_hi_ph[i]));
        nimages_per_dim[i][1]=-static_cast<int>(floor(s_lo_ph[i]));
        max_ntrial_atms*=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];
        sz+=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];

    }
    
    *s_trials=new type0[sz];
    for(int i=1;i<__dim__;i++)
        s_trials[i]=s_trials[i-1]+1+nimages_per_dim[i-1][0]+nimages_per_dim[i-1][1];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_dismantle()
{
    delete [] *s_trials;
    *s_trials=NULL;
}



