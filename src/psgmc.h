#ifndef __MAPP__PSGMC_2__
#define __MAPP__PSGMC_2__
#include "pgcmc_2.h"
namespace MAPP_NS
{
    class PSGMC_2:public PGCMC_2
    {
    public:
        PSGMC_2(class AtomsMD*&, class ForceFieldMD*&,class DynamicMD*&,int,elem_type,elem_type,type0,type0,int,type0);
        ~PSGMC_2();
    private:
        void attmpt();
        void decide();
        void ex_succ();
    };




}
#endif
