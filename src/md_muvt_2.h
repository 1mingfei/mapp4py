#ifndef __MAPP__md_muvt_2__
#define __MAPP__md_muvt_2__
#include "md_nvt.h"
namespace MAPP_NS
{
    class MDMuVT_2:public MDNVT
    {
    private:
        int seed;
        type0 mu0;
        type0 mu1;
        std::string gas_elem_name0;
        std::string gas_elem_name1;
        type0 ratio;
        int nevery;
        int nattempts;

    protected:
        void update_x_d__x(type0);
        void update_x_d__x_w_dof(type0);
        void update_x_d();
        void update_x_d_w_dof();
        void pre_run_chk(AtomsMD*,ForceFieldMD*);
        void pre_init();
    public:
        MDMuVT_2(type0,type0,type0,std::string,int,std::string,type0,type0);
        ~MDMuVT_2();
        void init();
        void run(int);
        void fin();
        
        typedef struct
        {
            PyObject_HEAD
            MDMuVT_2* md;
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
        static void getset_nevery(PyGetSetDef&);
        static void getset_nattempts(PyGetSetDef&);
        static void getset_seed(PyGetSetDef&);
        static void getset_gas_element0(PyGetSetDef&);
        static void getset_gas_element1(PyGetSetDef&);
        
        
        static int setup_tp();
        
    };
}
#endif
