#include "atoms_md.h"
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsMD::AtomsMD(MPI_Comm& world):
Atoms(world),
x_d(NULL),
elem(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsMD::AtomsMD(Communication& comm):
Atoms(comm),
x_d(NULL),
elem(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
AtomsMD::~AtomsMD()
{
    delete x_d;
    delete elem;
}
/*--------------------------------------------
 
 --------------------------------------------*/
#include "random.h"
#include "elements.h"
#include "xmath.h"
void AtomsMD::create_T(type0 T,int seed)
{
    if(!x_d) x_d=new Vec<type0>(this,__dim__);
    Random rand(seed+comm_rank);
    type0* __x_d=x_d->begin();
    type0* m=elements->masses;
    elem_type* __elem=elem->begin();
    type0 fac;
    if(dof)
    {
        bool* __dof=dof->begin();
        for(int i=0;i<natms;i++)
        {
            fac=sqrt(kB*T/m[*__elem]);
            Algebra::Do<__dim__>::func([&fac,&rand,&__x_d,__dof,this](const int i){if(__dof[i]) __x_d[i]=rand.gaussian()*fac;});
            __x_d+=__dim__;
            __dof+=__dim__;
            ++__elem;
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            fac=sqrt(kB*T/m[*__elem]);
            Algebra::Do<__dim__>::func([&fac,&rand,&__x_d,this](const int i){__x_d[i]=rand.gaussian()*fac;});
            __x_d+=__dim__;
            ++__elem;
        }
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
#include "ff_styles.h"
PyObject* AtomsMD::__new__(PyTypeObject* type,PyObject* args,PyObject* kwds)
{
    Object* __self=reinterpret_cast<Object*>(type->tp_alloc(type,0));
    PyObject* self=reinterpret_cast<PyObject*>(__self);
    return self;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AtomsMD::__init__(PyObject* self,PyObject* args,PyObject* kwds)
{
    FuncAPI<> func("__init__");
    if(func(args,kwds)==-1) return -1;
    Object* __self=reinterpret_cast<Object*>(self);
    __self->atoms=NULL;
    __self->ff=NULL;
    return 0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
PyObject* AtomsMD::__alloc__(PyTypeObject* type,Py_ssize_t)
{
    Object* __self=new Object;
    __self->ob_type=type;
    __self->ob_refcnt=1;
    __self->atoms=NULL;
    __self->ff=NULL;
    return reinterpret_cast<PyObject*>(__self);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsMD::__dealloc__(PyObject* self)
{
    Object* __self=reinterpret_cast<Object*>(self);
    delete __self->ff;
    delete __self->atoms;
    delete __self;
}
/*--------------------------------------------*/
PyTypeObject AtomsMD::TypeObject ={PyObject_HEAD_INIT(NULL)};
/*--------------------------------------------*/
void AtomsMD::setup_tp()
{
    TypeObject.tp_name="atoms_md";
    TypeObject.tp_doc="I will add doc here";
    
    TypeObject.tp_flags=Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    TypeObject.tp_basicsize=sizeof(Object);
    
    TypeObject.tp_new=__new__;
    TypeObject.tp_init=__init__;
    TypeObject.tp_alloc=__alloc__;
    TypeObject.tp_dealloc=__dealloc__;
    
    setup_tp_getset();
    TypeObject.tp_getset=getset;
    setup_tp_methods();
    TypeObject.tp_methods=methods;
    
    TypeObject.tp_base=&Atoms::TypeObject;
}
/*--------------------------------------------*/
PyGetSetDef AtomsMD::getset[]={[0 ... 0]={NULL,NULL,NULL,NULL,NULL}};
/*--------------------------------------------*/
void AtomsMD::setup_tp_getset()
{
    
}
/*--------------------------------------------*/
PyMethodDef AtomsMD::methods[]={[0 ... 6]={NULL,NULL,0,NULL}};
/*--------------------------------------------*/
void AtomsMD::setup_tp_methods()
{
    ml_create_T(methods[0]);
    ForceFieldLJ::ml_new(methods[1]);
    ForceFieldEAM::ml_new(methods[2],methods[3],methods[4]);
    ForceFieldFS::ml_new(methods[5]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AtomsMD::ml_create_T(PyMethodDef& tp_method)
{
    tp_method.ml_flags=METH_VARARGS | METH_KEYWORDS;
    tp_method.ml_name="create_T";
    tp_method.ml_doc="I will add doc here";
    
    tp_method.ml_meth=(PyCFunction)(PyCFunctionWithKeywords)
    [](PyObject* self,PyObject* args,PyObject* kwds)->PyObject*
    {
        FuncAPI<type0,int> f("create_T",{"T","seed"});
        f.logics<0>()[0]=VLogics("gt",0.0);
        f.logics<1>()[0]=VLogics("gt",0);
        if(f(args,kwds)) return NULL;
        
        AtomsMD::Object* __self=reinterpret_cast<AtomsMD::Object*>(self);
        __self->atoms->create_T(f.val<0>(),f.val<1>());
        Py_RETURN_NONE;
    };
}
