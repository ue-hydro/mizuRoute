// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"


// Constructor
// initalize numHRUs value
ClassWQ_OpenWQ::ClassWQ_OpenWQ() {}

// Deconstructor
ClassWQ_OpenWQ::~ClassWQ_OpenWQ() {}

//time_t ClassWQ_OpenWQ::convert_time(
//    int year, 
//    int month, 
//    int day, 
//    int hour, 
//    int minute) {

//    std::time_t sim_time;
//    std::tm tm{};
//    tm.tm_year = year - 1900; // -1900 is needed to get the conversion to produce the correct output
//    tm.tm_mon = month - 1;
//    tm.tm_hour = hour;
//    tm.tm_mday = day;
//    tm.tm_min = minute;
//    sim_time = timegm(&tm);

//    return sim_time;
//}

int ClassWQ_OpenWQ::decl(
    int nRch
    //int num_HRU,                // num HRU
    //int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
    //int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
    //int nSoil_2openwq,        // num layers of snoil (variable)
    //int nRunoff_2openwq,      // num layers in the runoff of SUMMA
    //int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    //int nYdirec_2openwq       // num of layers in y-dir (set to 1 because not used in summa)
    ){           
    
    OpenWQ_hostModelconfig_ref = new OpenWQ_hostModelconfig();
    OpenWQ_couplercalls_ref = new OpenWQ_couplercalls();
    OpenWQ_json_ref = new OpenWQ_json();
    OpenWQ_wqconfig_ref = new OpenWQ_wqconfig();
    OpenWQ_units_ref = new OpenWQ_units();
    OpenWQ_utils_ref = new OpenWQ_utils();
    OpenWQ_readjson_ref = new OpenWQ_readjson();
    OpenWQ_initiate_ref = new OpenWQ_initiate();
    OpenWQ_watertransp_ref = new OpenWQ_watertransp();
    OpenWQ_chem_ref = new OpenWQ_chem();
    OpenWQ_extwatflux_ss_ref = new OpenWQ_extwatflux_ss();
    OpenWQ_output_ref = new OpenWQ_output();
    
    //this->num_HRU = num_HRU;
    this->nRch = nRch;

    if (OpenWQ_hostModelconfig_ref->HydroComp.size()==0) {

        // Compartment names
        // Make sure to use capital letters for compartment names
        OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(rivernetwork_nRch_openwq,"RIVER_NETWORK_REACHES", nRch, 1, 1));      // Canopy
    //    OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(snow_index_openwq,"ILAYERVOLFRACWAT_SNOW", num_HRU, nYdirec_2openwq, max_snow_layers));  // snow (layerd)
    //    OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(runoff_index_openwq,"RUNOFF", num_HRU, nYdirec_2openwq, nRunoff_2openwq));               // Runoff
    //    OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(soil_index_openwq,"ILAYERVOLFRACWAT_SOIL", num_HRU, nYdirec_2openwq, nSoil_2openwq));    // Soil (layerd)
    //    OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(aquifer_index_openwq,"SCALARAQUIFER", num_HRU, nYdirec_2openwq, nAquifer_2openwq));      // Aquifer


    //    OpenWQ_vars_ref = new OpenWQ_vars(OpenWQ_hostModelconfig_ref->HydroComp.size());

        // External fluxes
        // Make sure to use capital letters for external fluxes
    //    OpenWQ_hostModelconfig_ref->HydroExtFlux.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"PRECIP",   num_HRU,nYdirec_2openwq,1));

        // Dependencies
        // to expand BGC modelling options
    //    OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"SM",        num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq));
    //    OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(1,"Tair_K",    num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq));
    //    OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(2,"Tsoil_K",   num_HRU,nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq));

        // Master Json
    //    OpenWQ_wqconfig_ref->OpenWQ_masterjson = "openWQ_master.json";


    //    OpenWQ_couplercalls_ref->InitialConfig(
    //        *OpenWQ_hostModelconfig_ref,
    //        *OpenWQ_json_ref,                // create OpenWQ_json object
    //        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
    //        *OpenWQ_units_ref,               // functions for unit conversion
    //        *OpenWQ_utils_ref,                // utility methods/functions
    //        *OpenWQ_readjson_ref,            // read json files
    //        *OpenWQ_vars_ref,
    //        *OpenWQ_initiate_ref,            // initiate modules
    //        *OpenWQ_watertransp_ref,         // transport modules
    //        *OpenWQ_chem_ref,                // biochemistry modules
    //        *OpenWQ_extwatflux_ss_ref,       // sink and source modules)
    //        *OpenWQ_output_ref);
            
    }
    return 0;
}


