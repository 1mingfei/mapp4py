/*
 * author: 1mingfei
 * date: 01/30/2019
 * goal: semi-carnonical ensemble implementation
 */
#include "atoms_md.h"
#include "pgcmc_2.h"
#include "psgmc_2.h"
#include "memory.h"
#include "random.h"
#include "elements.h"
#include "neighbor_md.h"
#include "ff_md.h"
#include "xmath.h"
#include "MAPP.h"
#include "comm.h"
#include "dynamic_md.h"
#include "print.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
PSGMC_2::PSGMC_2(AtomsMD*& __atoms, ForceFieldMD*& __ff,\
               DynamicMD*& __dynamic, int __m, elem_type __gas_type0, elem_type __gas_type1, type0 __mu0,\
               type0 __mu1, int seed, type0 __T):
    PGCMC_2(__atoms, __ff, __dynamic, __m, __gas_type0, __mu0, __T, seed, __gas_type1, __mu1, 0.5),
{
    n_comm=0;
    comm_id=NULL;
    success=NULL;
    ntrial_atms=NULL;
    gcmc_mode=NULL;
    gcmc_world=NULL;
    curr_comms=NULL;
    vars_comm=NULL;
    lcl_vars_comm=NULL;
    int_buff=NULL;
    roots=NULL;
    comms=NULL;
    
    s_buff=NULL;
    s_x_buff=NULL;
    cell_coord_buff=NULL;
    
    n_cells=0;
    head_atm=NULL;
    
    max_ntrial_atms=0;
    cell_coord_buff=NULL;
    s_x_buff=NULL;
    s_buff=NULL;

    lcl_random=new Random(seed+atoms->comm_rank);
    
    /*--------------------------------------------------
     find the relative neighbor list for cells here we 
     figure out the number of neighboring cells and also 
     allocate the memory for rel_neigh_lst and 
     rel_neigh_lst_coord. values for rel_neigh_lst_coord
     is assigned here,but finding the values for 
     rel_neigh_lst requires knowledge of the box and 
     domain, it will be differed to create()
     --------------------------------------------------*/
    int countr[__dim__];
    for(int i=0;i<__dim__;i++)
        countr[i]=-m;
    int max_no_neighs=1;
    for(int i=0;i<__dim__;i++)
        max_no_neighs*=2*m+1;
    nneighs=0;
    rel_neigh_lst_coord=new int[max_no_neighs*__dim__];
    int* rel_neigh_lst_coord_=rel_neigh_lst_coord;
    int sum;
    int rc_sq=m*m;
    for(int i=0;i<max_no_neighs;i++)
    {
        sum=0;
        for(int j=0;j<__dim__;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*std::abs(countr[j]);
        }
        
        if(sum<rc_sq)
        {
            for(int j=0;j<__dim__;j++)
                rel_neigh_lst_coord_[j]=countr[j];
            rel_neigh_lst_coord_+=__dim__;
            nneighs++;
        }
        
        countr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(countr[j]==m+1)
            {
                countr[j]=-m;
                countr[j+1]++;
            }
    }
    rel_neigh_lst_coord_=new int[nneighs*__dim__];
    memcpy(rel_neigh_lst_coord_,rel_neigh_lst_coord,nneighs*__dim__*sizeof(int));
    delete [] rel_neigh_lst_coord;
    rel_neigh_lst_coord=rel_neigh_lst_coord_;
}

/*--------------------------------------------
 destructor
 --------------------------------------------*/
PSGMC_2::~PSGMC_2()
{
    delete [] rel_neigh_lst_coord;
    delete lcl_random;   
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PSGMC_2::attmpt()
{
    create_comm_pattern();
    if(im_root)
    {
        if(ngas_lclArr[0] || ngas_lclArr[1]) //if existing any lcl gas 
        {
            gcmc_mode[0]=EX_MODE;
            comm_buff[0]=3;
            //use total lcl gas element instead
            igas=static_cast<int>((ngas_lclArr[0]+ngas_lclArr[1])*lcl_random->uniform());

            if (ngas_lclArr[0])
            {
                if (igas < ngas_lclArr[0])
                {
                    n = igas;
                    igasType = 0;
                }
                else
                {
                    n = (igas - ngas_lclArr[0]);
                    igasType = 1;
                }
            }
            else
            {
                n = igas;
                igasType = 1;
            }
            
            int icount=-1;
            elem_type* elem=atoms->elem->begin();
            // I don't know why sina do this, but I will keep it (start)
            del_idx=0;
            for(;icount!=n;del_idx++)
                if(elem[del_idx]==gas_typeArr[igasType]) icount++;
            del_idx--;
            gas_id=atoms->id->begin()[del_idx];
            // I don't know why sina do this, but I will keep it (end)

            //this can be optimized
            memcpy(s_buff[0],s_vec_p->begin()+del_idx*__dim__,__dim__*sizeof(type0));
            memcpy(comm_buff+1,s_buff[0],__dim__*sizeof(type0));
        }
        else
        {
            gcmc_mode[0]=NOEX_MODE;
            comm_buff[0]=0;
        }

        
        MPI_Bcast(comm_buff,comm_buff_size,MPI_BYTE,roots[0],*curr_comms[0]);

    }
    else
        for(int i=0;i<n_curr_comms;i++)
        {
            MPI_Bcast(comm_buff,comm_buff_size,MPI_BYTE,roots[i],*curr_comms[i]);
            memcpy(s_buff[i],comm_buff+1,__dim__*sizeof(type0));
            if(comm_buff[0]==3)
                gcmc_mode[i]=EX_MODE;
            else
                gcmc_mode[i]=NOEX_MODE;
        }
    

    prep_s_x_buff();

    if(tag_vec_p) reset_tag();

    ff->pre_xchng_energy_timer(this);

    delta_u=ff->xchng_energy_timer(this);

    root_succ=false;
    if(im_root) decide();

    for(int i=0;i<n_curr_comms;i++)
        MPI_Bcast(&success[i],1,MPI_INT,roots[i],*curr_comms[i]);

    if(tag_vec_p) success2tag();

    ff->post_xchng_energy_timer(this);

    finalize();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PSGMC_2::decide()
{
    if(prev_p!=-1)
    {
        MPI_Recv(&int_buff_sz,1,MPI_INT,prev_p,0,*gcmc_world,MPI_STATUS_IGNORE);
        MPI_Recv(int_buff,int_buff_sz,MPI_INT,prev_p,1,*gcmc_world,MPI_STATUS_IGNORE);
    }
    else
    {
        int_buff_sz=2;
        int_buff[0]=0;
        int_buff[1]=del_ids_sz;
    }
    
    success[0]=-1;
    type0 fac;

    if(gcmc_mode[0]==EX_MODE)
    {
        //igasType : the gas type TO BE changed
        //if 0 (false), then we are flipping 0 -> 1;
        //if 1 (true), then we are flipping 1 -> 0;
        type0 deltaMu = igasType ? (__mu0 - __mu1) : (__mu1 - __mu0);
        fac=static_cast<type0>(ngas_lclArr[igasType])*exp(beta*(delta_u-deltaMu))/(z_facArr[igasType]*vol_lcl);
        if(lcl_random->uniform()<fac)
        {
//#ifdef GCMCDEBUG
//            tot_delta_u_lcl-=delta_u;
//#endif
            root_succ=true;
            success[0]=0;
            /* int_buff[0]--;
            int_buff[int_buff_sz]=gas_id;
            int_buff_sz++;
            */
        }
    }
    
    if(next_p!=-1)
    {
        MPI_Send(&int_buff_sz,1,MPI_INT,next_p,0,*gcmc_world);
        MPI_Send(int_buff,int_buff_sz,MPI_INT,next_p,1,*gcmc_world);
    }
    
    if(root_succ)
    {
        if(gcmc_mode[0]==EX_MODE) ex_succ();
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PSGMC_2::ex_succ()
{
    atoms->elem->begin()[del_idx]=gas_typeArr[igasType];

    if(tag_vec_p) tag_vec_p->begin()[natms_lcl-1]=-1;
    if(!dof_empty)
    {
        bool* dof=atoms->x_dof->begin()+(natms_lcl-1)*__dim__;
        for(int i=0;i<__dim__;i++) dof[i]=true;
    }
    
    
    Algebra::DyadicV<__dim__>(gas_massArr[igasType],atoms->x_d->begin()+(natms_lcl-1)*__dim__,mvv_lcl);
    memcpy(s_vec_p->begin()+(natms_lcl-1)*__dim__,s_buff[0],__dim__*sizeof(type0));
    
}

