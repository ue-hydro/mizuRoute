#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"
/**
 * Below is the implementation of the C interface for SUMMA. When Summa calls a function 
 * the functions below are the ones that are invoked first. 
 * The openWQ object is then passed from Fortran to these functions so that the OpenWQ object
 * can be called. The openWQ object methods are defined above.
 */
// Interface functions to create Object
CLASSWQ_OPENWQ* create_openwq() {
    return new ClassWQ_OpenWQ();
}

void delete_openwq(CLASSWQ_OPENWQ* openWQ) {
    delete openWQ;
}

int openwq_decl(
    ClassWQ_OpenWQ *openWQ,
    int nRch
    //int hruCount,              // num HRU
    //int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
    //int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
    //int nSoil_2openwq,        // num layers of snoil (variable)
    //int nRunoff_2openwq,      // num layers of runoff (fixed to 1)
    //int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    //int nYdirec_2openwq       // num of layers in y-dir (set to 1 because not used in summa)
    ){            

    return openWQ->decl(
        nRch
        //hruCount, 
        //nCanopy_2openwq, 
        //nSnow_2openwq, 
        //nSoil_2openwq, 
        //nRunoff_2openwq,
        //nAquifer_2openwq, 
        //nYdirec_2openwq
        );

}

