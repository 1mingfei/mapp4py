#ifndef __MAPP__min_cg_fit__
#define __MAPP__min_cg_fit__
#include "min.h"
#include "min_vec.h"
#include "ff_md.h"
#include "atoms_md.h"
#include "dynamic_md.h"
#include "MAPP.h"
namespace MAPP_NS
{
    class MinCGFit:public Min
    {
    private:
    protected:
        
        void prep();
        
        VecTens<type0,1> h;
        VecTens<type0,1> x;
        VecTens<type0,1> x0;
        VecTens<type0,1> x_d;
        VecTens<type0,1> f;
        VecTens<type0,1> f0;
        vec* ext_vec_0;
        vec* ext_vec_1;
        type0 S_tmp[__dim__][__dim__];
        
        
        
    public:
        MinCGFit(type0,bool(&)[__dim__][__dim__],bool,type0,class LineSearch*,vec*,vec*);
        ~MinCGFit();
        virtual void run(int);
        template<class C>
        void run(C*,int);
        virtual void init();
        virtual void fin();
        
        void force_calc();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();

        
        
        class AtomsMD* atoms;
        class ForceFieldMD* ff;
        class ExportMD* xprt;
        class DynamicMD* dynamic;
        
        typedef struct
        {
            PyObject_HEAD
            MinCGFit* min;
            LineSearch::Object* ls;
            ExportMD::Object* xprt;
        }Object;
        
        static PyTypeObject TypeObject;
        static PyObject* __new__(PyTypeObject*,PyObject*, PyObject*);
        static int __init__(PyObject*, PyObject*,PyObject*);
        static PyObject* __alloc__(PyTypeObject*,Py_ssize_t);
        static void __dealloc__(PyObject*);
        
        static PyMethodDef methods[];
        static void setup_tp_methods();
        static void ml_run(PyMethodDef&);
        
        static PyGetSetDef getset[];
        static void setup_tp_getset();
        static void getset_export(PyGetSetDef&);
        
        static int setup_tp();
        
        
        
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class C>
void MinCGFit::run(C* ls,int nsteps)
{
    int step=atoms->step;
    
    force_calc();
    
    int nevery_xprt=xprt==NULL ? 0:xprt->nevery;
    if(nevery_xprt) xprt->write(step);
    

    
    type0 e_prev,e_curr=atoms->pe;
    type0 f0_f0,f_f,f_f0;
    type0 ratio,alpha;
    int err=nsteps==0? MIN_F_MAX_ITER:LS_S;
    h=f;
    f0_f0=f*f;
    int istep=0;
    for(;err==LS_S;istep++)
    {
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }
        
        x0=x;
        f0=f;
        
        e_prev=e_curr;
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            f_h=f0_f0;
        }
        prep();
        err=ls->min(this,e_curr,alpha,1);
        
        if(err!=LS_S)
        {
            // this was a bullshit step so we have to decrease the setup once to compensate
            // the last step was the previous one
            istep--;
            continue;
        }
        
        force_calc();
        
        if(nevery_xprt && (istep+1)%nevery_xprt==0) xprt->write(step+istep+1);
        //this was a successfull step but the last one
        if(e_prev-e_curr<e_tol) err=MIN_S_TOLERANCE;
        if(istep+1==nsteps) err=MIN_F_MAX_ITER;
        if(err) continue;
        
        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h=ratio*h+f;
        f0_f0=f_f;
    }
    

    
    
    if(nevery_xprt && istep%nevery_xprt) xprt->write(step+istep);
    

    
    
    atoms->step+=istep;
}
#endif
