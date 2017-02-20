#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL ARRAY_API
#include <numpy/arrayobject.h>
#include "MAPP.h"
#include <mpi.h>
#include "comm.h"
#include "example.h"
#include "atoms_styles.h"
#include "md_styles.h"
#include "min_styles.h"
#include "read_styles.h"
#define GET_FILE(file_name) reinterpret_cast<PyFileObject*>(PySys_GetObject((char*)#file_name))->f_fp
using namespace MAPP_NS;
/*--------------------------------------------*/
PyMethodDef MAPP::methods[]={[0 ... 2]={NULL}};
/*--------------------------------------------*/
void MAPP::setup_methods()
{
    methods[0].ml_name="pause_slave_out";
    methods[0].ml_meth=(PyCFunction)pause_out;
    methods[0].ml_flags=METH_NOARGS;
    methods[0].ml_doc="pauses stdout & stderr of non-root processes";
    
    methods[1].ml_name="resume_slave_out";
    methods[1].ml_meth=(PyCFunction)resume_out;
    methods[1].ml_flags=METH_NOARGS;
    methods[1].ml_doc="resumes stdout & stderr of non-root processes";
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::pause_out(PyObject* self)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(!glbl_rank) Py_RETURN_NONE;
    if(!__devnull__) __devnull__=fopen(devnull_path,"w");

    GET_FILE(stdout)=__devnull__;
    GET_FILE(stderr)=__devnull__;
    mapp_out=__devnull__;
    mapp_err=__devnull__;
    Py_RETURN_NONE;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::resume_out(PyObject*)
{
    MPI_Barrier(MPI_COMM_WORLD);
    GET_FILE(stdout)=__stdout__;
    GET_FILE(stderr)=__stderr__;
    mapp_out=__stdout__;
    mapp_err=__stderr__;
    if(__devnull__)
    {
        fclose(__devnull__);
        __devnull__=NULL;
    }
    Py_RETURN_NONE;
}
/*--------------------------------------------*/
const char* MAPP::devnull_path;
int MAPP::glbl_rank=0;
FILE* MAPP_NS::MAPP::__devnull__(NULL);
FILE* MAPP_NS::MAPP::__stdout__(NULL);
FILE* MAPP_NS::MAPP::__stderr__(NULL);
FILE* MAPP_NS::MAPP::mapp_out(NULL);
FILE* MAPP_NS::MAPP::mapp_err(NULL);
/*--------------------------------------------*/
void MAPP::init_module(void)
{
    PyObject* posixpath=PyImport_ImportModule("posixpath");
    PyObject* devnull_path_op=PyObject_GetAttrString(posixpath,"devnull");
    devnull_path=PyString_AsString(devnull_path_op);
    Py_DECREF(devnull_path_op);
    Py_DECREF(posixpath);
    glbl_rank=Communication::get_rank();
    
    mapp_out=__stdout__=GET_FILE(stdout);
    mapp_err=__stderr__=GET_FILE(stderr);
    
    import_array();
    setup_methods();
    PyObject* module=Py_InitModule3("mapp",methods,"MIT Atomistic Parallel Package");
    if(module==NULL) return;
    
    MAPP_MPI::setup_tp();
    if(PyType_Ready(&MAPP_MPI::TypeObject)<0) return;
    Py_INCREF(&MAPP_MPI::TypeObject);
    PyModule_AddObject(module,"mpi",reinterpret_cast<PyObject*>(&MAPP_MPI::TypeObject));
    
    ExamplePython::setup_tp();
    if(PyType_Ready(&ExamplePython::TypeObject)<0) return;
    Py_INCREF(&ExamplePython::TypeObject);
    PyModule_AddObject(module,"xmpl",reinterpret_cast<PyObject*>(&ExamplePython::TypeObject));
    
    PyObject* md=MAPP::MD::init_module();
    if(md==NULL) return;
    PyModule_AddObject(module,"md",md);
    
    PyObject* dmd=MAPP::DMD::init_module();
    if(dmd==NULL) return;
    PyModule_AddObject(module,"dmd",dmd);
}
/*--------------------------------------------*/
PyMethodDef MAPP::MD::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void MAPP::MD::setup_methods()
{
    ReadCFGMD::ml_cfg(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::MD::init_module(void)
{
    setup_methods();
    PyObject* module=Py_InitModule3("md",methods,"Molecular Dynamics (MD) module");
    if(module==NULL) return NULL;
    
    AtomsMD::setup_tp();
    if(PyType_Ready(&AtomsMD::TypeObject)<0) return NULL;
    Py_INCREF(&AtomsMD::TypeObject);
    PyModule_AddObject(module,"atoms",reinterpret_cast<PyObject*>(&AtomsMD::TypeObject));
    
    MDNVT::setup_tp();
    if(PyType_Ready(&MDNVT::TypeObject)<0) return NULL;
    Py_INCREF(&MDNVT::TypeObject);
    PyModule_AddObject(module,"nvt",reinterpret_cast<PyObject*>(&MDNVT::TypeObject));
    
    
    MDNST::setup_tp();
    if(PyType_Ready(&MDNST::TypeObject)<0) return NULL;
    Py_INCREF(&MDNST::TypeObject);
    PyModule_AddObject(module,"nst",reinterpret_cast<PyObject*>(&MDNST::TypeObject));
    
    
    MDMuVT::setup_tp();
    if(PyType_Ready(&MDMuVT::TypeObject)<0) return NULL;
    Py_INCREF(&MDMuVT::TypeObject);
    PyModule_AddObject(module,"muvt",reinterpret_cast<PyObject*>(&MDMuVT::TypeObject));
    
    
    MinCG::setup_tp();
    if(PyType_Ready(&MinCG::TypeObject)<0) return NULL;
    Py_INCREF(&MinCG::TypeObject);
    PyModule_AddObject(module,"min_cg",reinterpret_cast<PyObject*>(&MinCG::TypeObject));
    
    
    MinLBFGS::setup_tp();
    if(PyType_Ready(&MinLBFGS::TypeObject)<0) return NULL;
    Py_INCREF(&MinLBFGS::TypeObject);
    PyModule_AddObject(module,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGS::TypeObject));
    
    return module;
}
/*--------------------------------------------*/
PyMethodDef MAPP::DMD::methods[]={[0 ... 1]={NULL}};
/*--------------------------------------------*/
void MAPP::DMD::setup_methods()
{
    ReadCFGDMD::ml_cfg(methods[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* MAPP::DMD::init_module(void)
{
    setup_methods();
    PyObject* module=Py_InitModule3("dmd",methods,"Diffusive Molecular Dynamics (DMD) module");
    if(module==NULL) return NULL;
    
    AtomsDMD::setup_tp();
    if(PyType_Ready(&AtomsDMD::TypeObject)<0) return NULL;
    Py_INCREF(&AtomsDMD::TypeObject);
    PyModule_AddObject(module,"atoms",reinterpret_cast<PyObject*>(&AtomsDMD::TypeObject));
    
    MinCGDMD::setup_tp();
    if(PyType_Ready(&MinCGDMD::TypeObject)<0) return NULL;
    Py_INCREF(&MinCGDMD::TypeObject);
    PyModule_AddObject(module,"min_cg",reinterpret_cast<PyObject*>(&MinCGDMD::TypeObject));
    
    MinLBFGSDMD::setup_tp();
    if(PyType_Ready(&MinLBFGSDMD::TypeObject)<0) return NULL;
    Py_INCREF(&MinLBFGSDMD::TypeObject);
    PyModule_AddObject(module,"min_lbfgs",reinterpret_cast<PyObject*>(&MinLBFGSDMD::TypeObject));
    
    return module;
}
