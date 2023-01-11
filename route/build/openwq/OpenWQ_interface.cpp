#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"
/**
 * Below is the implementation of the C interface for SUMMA. When Summa calls a function 
 * the functions below are the ones that are invoked first. 
 * The openWQ object is then passed from Fortran to these functions so that the OpenWQ object
 * can be called. The openWQ object methods are defined above.
 */
// Interface functions to create Object
CLASSWQ_openwq* create_openwq() {
    return new CLASSWQ_openwq();
}

void delete_openwq(CLASSWQ_openwq* openWQ) {
    delete openWQ;
}

int openwq_decl(
    CLASSWQ_openwq *openWQ,
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

int openwq_run_time_start(
    CLASSWQ_openwq *openWQ,
//    bool last_hru_flag, 
//    int hru_index, 
//    int nSnow_2openwq, 
//    int nSoil_2openwq, 
      int simtime_mizuroute[] 
//    double soilMoist_depVar_summa_frac[],                  
//    double soilTemp_depVar_summa_K[],
//    double airTemp_depVar_summa_K,
//    double sweWatVol_stateVar_summa_m3[],
//    double canopyWatVol_stateVar_summa_m3,
//    double soilWatVol_stateVar_summa_m3[],
//    double aquiferWatVol_stateVar_summa_m3
) {
    
    return openWQ->openwq_run_time_start(
//        last_hru_flag,
//        hru_index, 
//        nSnow_2openwq, 
//        nSoil_2openwq, 
          simtime_mizuroute 
//        soilMoist_depVar_summa_frac,                  
//        soilTemp_depVar_summa_K,
//        airTemp_depVar_summa_K,
//        sweWatVol_stateVar_summa_m3,
//        canopyWatVol_stateVar_summa_m3,
//        soilWatVol_stateVar_summa_m3,
//        aquiferWatVol_stateVar_summa_m3
    );   
}

int openwq_run_time_end(
    CLASSWQ_openwq *openWQ, 
    int simtime_mizuroute[]) {

    return openWQ->openwq_run_time_end(
        simtime_mizuroute);
}
