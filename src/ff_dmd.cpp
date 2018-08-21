#include "ff_dmd.h"
#include "atoms_dmd.h"
#include "neighbor_dmd.h"
#include "xmath.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::ForceFieldDMD(AtomsDMD* __atoms):
ForceField(__atoms),
atoms(__atoms),
ave_mu(NULL),
rsq_crd(NULL),
r_crd(NULL)
{
    Memory::alloc(ave_mu,nelems);
    Memory::alloc(cut_sk,nelems,nelems);
    Memory::alloc(rsq_crd,nelems);
    Memory::alloc(r_crd,nelems);
    neighbor=new NeighborDMD(__atoms,cut_sk,rsq_crd);
    f=new Vec<type0>(atoms,__dim__,"f");
    f_alpha=new DMDVec<type0>(atoms,0.0,"f_alpha");
    c_d=new DMDVec<type0>(atoms,0.0,"c_d");
    mu=new DMDVec<type0>(atoms,0.0,"mu");
}
/*--------------------------------------------
 
 --------------------------------------------*/
ForceFieldDMD::~ForceFieldDMD()
{
    delete mu;
    delete c_d;
    delete f_alpha;
    delete f;
    delete neighbor;
    Memory::dealloc(r_crd);
    Memory::dealloc(rsq_crd);
    Memory::dealloc(cut_sk);
    Memory::dealloc(ave_mu);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::calc_ndof()
{
    nx_dof=atoms->natms*__dim__;
    if(!atoms->x_dof->is_empty())
    {
        const int n=atoms->natms_lcl*__dim__;
        int nx_dof_lcl=n;
        bool* x_dof=atoms->x_dof->begin();
        for(int i=0;i<n;i++) if(!x_dof[i]) nx_dof_lcl--;
        MPI_Allreduce(&nx_dof_lcl,&nx_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
    
    nalpha_dof=atoms->natms*atoms->c_dim;
    if(!atoms->alpha_dof->is_empty())
    {
        const int n=atoms->natms_lcl*atoms->c_dim;
        int nalpha_dof_lcl=n;
        bool* alpha_dof=atoms->alpha_dof->begin();
        for(int i=0;i<n;i++) if(!alpha_dof[i]) nalpha_dof_lcl--;
        MPI_Allreduce(&nalpha_dof_lcl,&nalpha_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
    
    
    nc_dof=atoms->natms*atoms->c_dim;
    if(!atoms->c_dof->is_empty())
    {
        const int n=atoms->natms_lcl*atoms->c_dim;
        int nc_dof_lcl=n;
        bool* c_dof=atoms->c_dof->begin();
        for(int i=0;i<n;i++) if(!c_dof[i]) nc_dof_lcl--;
        MPI_Allreduce(&nc_dof_lcl,&nc_dof,1,Vec<int>::MPI_T,MPI_SUM,world);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::reset()
{
    type0* __f=f->begin();
    const int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++) __f[i]=0.0;
    __f=f_alpha->begin();
    const int m=atoms->natms_lcl*c_dim;
    for(int i=0;i<m;i++) __f[i]=0.0;
    for(int i=0;i<__nvoigt__+2;i++) __vec_lcl[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::impose_dof(type0* __f,type0* __f_alpha)
{
    if(!dof_empty)
    {
        bool* dof=atoms->x_dof->begin();
        const int n=atoms->natms_lcl*__dim__;
        for(int i=0;i<n;i++) __f[i]=dof[i] ? __f[i]:0.0;
    }
    
    if(!dof_alpha_empty)
    {
        bool* dof=atoms->alpha_dof->begin();
        const int n=atoms->natms_lcl*c_dim;
        for(int i=0;i<n;i++) __f_alpha[i]=dof[i] ? __f_alpha[i]:0.0;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::norm_sq(type0* __f,type0* __f_alpha)
{
    type0 err_lcl=0.0,err;
    int n=atoms->natms_lcl*__dim__;
    for(int i=0;i<n;i++) err_lcl+=__f[i]*__f[i];
    n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) err_lcl+=__f_alpha[i]*__f_alpha[i];
    
    MPI_Allreduce(&err_lcl,&err,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return err;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::reset_c_d()
{
    type0* __c_d=c_d->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++) __c_d[i]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::pre_init()
{
    dof_empty=atoms->x_dof->is_empty();
    dof_alpha_empty=atoms->alpha_dof->is_empty();
    dof_c_empty=atoms->c_dof->is_empty();
    
    type0 tmp;
    max_cut=0.0;
    for(size_t i=0;i<nelems;i++)
    {
        for(size_t j=0;j<i+1;j++)
        {
            tmp=cut[i][j];
            max_cut=MAX(max_cut,tmp);
            cut_sq[i][j]=cut_sq[j][i]=tmp*tmp;
            tmp+=atoms->comm.skin;
            cut_sk[i][j]=cut_sk[j][i]=tmp;
            cut_sk_sq[i][j]=cut_sk_sq[j][i]=tmp*tmp;
        }
        rsq_crd[i]=r_crd[i]*r_crd[i];
    }
    c_dim=atoms->c->dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::post_fin()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::value()
{
    //timer->start(FORCE_TIME_mode);
    __vec_lcl[0]=__vec_lcl[1+__nvoigt__]=0.0;
    __energy_calc();
    type0 en;
    MPI_Allreduce(&__vec_lcl[0],&en,1,Vec<type0>::MPI_T,MPI_SUM,world);
    //timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::derivative()
{
    reset();
    __force_calc();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];

    impose_dof(f->begin(),f_alpha->begin());
    return __vec;
}
#include "dynamic_dmd.h"
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::calc_thermo()
{
    DynamicDMD __dynamic(atoms,this,false,{},{},{});
    __dynamic.init();
    calc_ndof();
    derivative();
    
    
    
    type0* mu_sum_lcl=NULL;
    Memory::alloc(mu_sum_lcl,nelems);
    for(size_t i=0;i<nelems;i++) mu_sum_lcl[i]=ave_mu[i]=0.0;
    
    type0* sum_lcl=NULL;
    type0* sum=NULL;
    Memory::alloc(sum_lcl,nelems);
    Memory::alloc(sum,nelems);
    for(size_t i=0;i<nelems;i++) sum_lcl[i]=sum[i]=0.0;
    
    
    type0* mu_vec=mu->begin();
    type0 const* c=atoms->c->begin();
    elem_type const* elem_vec=atoms->elem->begin();
    const int n=atoms->natms_lcl*c_dim;
    for(int i=0;i<n;i++)
    if(c[i]!=-1.0)
    {
        mu_sum_lcl[elem_vec[i]]+=mu_vec[i];
        ++sum_lcl[elem_vec[i]];
    }
    MPI_Allreduce(mu_sum_lcl,ave_mu,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    MPI_Allreduce(sum_lcl,sum,static_cast<int>(nelems),Vec<type0>::MPI_T,MPI_SUM,world);
    for(size_t i=0;i<nelems;i++) if(sum[i]!=0.0) ave_mu[i]/=sum[i];
    
    Memory::dealloc(sum_lcl);
    Memory::dealloc(sum);
    Memory::dealloc(mu_sum_lcl);
    
    __dynamic.fin();


}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldDMD::derivative_timer(type0(*&S)[__dim__])
//{
//    reset();
//    force_calc();
//    type0* fvec=f->begin();
//    type0* xvec=atoms->x->begin();
//    if(dof_empty)
//    {
//        const int natms_lcl=atoms->natms_lcl;
//        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
//            Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//    }
//    else
//    {
//        bool* dof=atoms->dof->begin();
//        const int natms_lcl=atoms->natms_lcl;
//        for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__,dof+=__dim__)
//        {
//            Algebra::Do<__dim__>::func([&dof,&fvec](int i){fvec[i]=dof[i] ? fvec[i]:0.0;});
//            Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//        }
//    }
//    
//    
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    
//    if(!dof_alpha_empty)
//    {
//        type0* fvec=f_alpha->begin();
//        bool* dof=atoms->dof_alpha->begin();
//        const int n=atoms->natms_lcl*c_dim;
//        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
//    }
//    
//    atoms->fe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
//    atoms->s=__vec[1+__nvoigt__];
//}
/*--------------------------------------------
 this does not sound right hs to be check later
 --------------------------------------------*/
//void ForceFieldDMD::derivative_timer(bool affine,type0(*&S)[__dim__])
//{
//    reset();
//    force_calc();
//    if(!affine)
//    {
//        type0* fvec=f->begin();
//        type0* xvec=atoms->x->begin();
//        if(dof_empty)
//        {
//            const int natms_lcl=atoms->natms_lcl;
//            for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__)
//                Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//        }
//        else
//        {
//            bool* dof=atoms->dof->begin();
//            const int natms_lcl=atoms->natms_lcl;
//            for(int i=0;i<natms_lcl;i++,fvec+=__dim__,xvec+=__dim__,dof+=__dim__)
//            {
//                Algebra::Do<__dim__>::func([&dof,&fvec](int i){fvec[i]=dof[i] ? fvec[i]:0.0;});
//                Algebra::DyadicV<__dim__>(xvec,fvec,__vec_lcl+1);
//            }
//        }
//    }
//    
//    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
//    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
//    Algebra::NONAME_DyadicV_mul_MLT(__vec+1,atoms->B,S);
//    const type0 vol=atoms->vol;
//    Algebra::Do<__nvoigt__>::func([this,&vol](int i){__vec[i+1]/=-vol;});
//    
//    if(!dof_alpha_empty)
//    {
//        type0* fvec=f_alpha->begin();
//        bool* dof=atoms->dof_alpha->begin();
//        const int n=atoms->natms_lcl*c_dim;
//        for(int i=0;i<n;i++) fvec[i]=dof[i] ? fvec[i]:0.0;
//    }
//    
//    atoms->fe=__vec[0];
//    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
//    atoms->s=__vec[1+__nvoigt__];
//}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::force_calc_static()
{
    reset();
    __force_calc_static();
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    
    impose_dof(f->begin(),f_alpha->begin());
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::c_d_calc()
{
    reset_c_d();
    __c_d_calc();
    
    impose_dof(f->begin(),f_alpha->begin());
    err=sqrt(norm_sq(f->begin(),f_alpha->begin()));
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::c_dd_norm()
{
    reset_c_d();
    __c_d_calc();
    __J(c_d,f_alpha);
    int n=atoms->natms_lcl*c_dim;
    type0* __f_alpha=f_alpha->begin();
    type0 norm_lcl=0.0;
    for(int i=0;i<n;i++)
        norm_lcl+=__f_alpha[i]*__f_alpha[i];
    type0 norm;
    MPI_Allreduce(&norm_lcl,&norm,1,Vec<type0>::MPI_T,MPI_SUM,world);
    return sqrt(norm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::J(Vec<type0>* x,Vec<type0>* Jx)
{
    __J(x,Jx);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0* ForceFieldDMD::J(Vec<type0>* x_ptr,Vec<type0>* alpha_ptr,Vec<type0>* Jx_ptr,Vec<type0>* Jalpha_ptr)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    __J(x_ptr,alpha_ptr,Jx_ptr,Jalpha_ptr);
    impose_dof(Jx_ptr->begin(),Jalpha_ptr->begin());
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__,Vec<type0>::MPI_T,MPI_SUM,world);
    return __vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::prepJ_n_res(Vec<type0>* x_ptr,Vec<type0>* alpha_ptr)
{
    Algebra::zero<__nvoigt__+2>(__vec_lcl);
    __prepJ_n_res(x_ptr,alpha_ptr);
    MPI_Allreduce(__vec_lcl,__vec,__nvoigt__+2,Vec<type0>::MPI_T,MPI_SUM,world);
    Algebra::Do<__nvoigt__>::func([this](int i){__vec[i+1]*=-1.0;});
    Algebra::DyadicV_2_MSY(__vec+1,F_H);
    atoms->fe=__vec[0];
    type0 vol_neg=-atoms->vol;
    Algebra::Do<__nvoigt__>::func([this,&vol_neg](int i){__vec[i+1]/=vol_neg;});
    Algebra::DyadicV_2_MSY(__vec+1,atoms->S_fe);
    atoms->s=__vec[1+__nvoigt__];
    
    impose_dof(x_ptr->begin(),alpha_ptr->begin());
    type0 err_sq=norm_sq(x_ptr->begin(),alpha_ptr->begin());
    err=sqrt(err_sq);
    return err_sq;
}





