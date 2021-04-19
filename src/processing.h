#ifndef __PROCESSING_H__
#define __PROCESSING_H__

#include "simulation.h"
#include "simulation_file.h"
#include "results.h"


// Run simulation with 'sim' configuration and 'num_photons_processed_'
// photons already processed, save received photons in 'file' and write
// results in 'results'. If 'exit_request' is true in any time, simulation is
// safely aborted. 'num_photons_processed_' and 'file' can be NULL
void run_simulation(Simulation * sim, Results * results, SimulationFile * file,
                    long * num_photons_processed_, bool * exit_request);

// Return results obtained from apply every configuration from array 'sim' of
// length 'sweep_length' on photons from simulation file 'file'
Results * get_results_from_file(Simulation * sim, SimulationFile * file,
                                unsigned sweep_length);


#endif