#include "processing.h"

#include <stdbool.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

#include "photon_block.h"
#include "simulation.h"
#include "status.h"
#include "urand.h"



void run_simulation(Simulation * sim, Results * results, SimulationFile * file,
                    long * num_photons_processed_, bool * exit_request)
{
    long * num_photons_processed; // Numbers of photons already processed
    long num_photons_remaining;   // Number of photons to be processed
    long num_photons;  // Photons to be processes in every task
    bool local_malloc; // If 'num_photons_processed_' is allocated locally

    // Set seeds of random module (one per thread)
    urand_init();
    local_malloc = false;

    //Init simulation n and boundary_pos array
    sim->med_boundary_pos = (float*)malloc((sim->med_layers-1)*sizeof(float));
    sim->med_n_water_variables = (float*)malloc((sim->med_layers)*sizeof(float));
    if(strcmp(sim->phase_json,"") != 0)
        initParametrosPantallaFaseFromJson(sim->phase_json, sim);
    else
        sim->phase_layers = 0;
    // If num_photons_processed_ == NULL, then there are no photons processed
    if(num_photons_processed_ != NULL)
        num_photons_processed = num_photons_processed_;
    else
    {
        num_photons_processed = malloc(sizeof(*num_photons_processed));
        local_malloc = true;
        *num_photons_processed = 0;
    }

    num_photons_remaining = sim->photons - *num_photons_processed;

    start_processing_time_tracking();
    // Parallel region starts
    #pragma omp parallel
    {
        #pragma omp single
        {
            while(num_photons_remaining > 0  &&  !*exit_request)
            {
                // Set the number of photons to be processed in this iteration
                if(num_photons_remaining > sim->block_size)
                    num_photons = sim->block_size;
                else
                    num_photons = num_photons_remaining;
                num_photons_remaining -= num_photons;

                // Create a task to process 'num_photons' photons
                #pragma omp task firstprivate(num_photons)
                {
                    PhotonBlock * block; // Block of photons to process
                    Results * lresults;  // Local results obtained from 'block'

                    // If there is an exit request, task cannot be executed
                    if(!*exit_request)
                    {
                        // Create the block
                        block = photon_block_new(num_photons);
                        init_water_n_and_boundarys(sim);
                        // Process block (emission, tranmission and reception)
                        photon_block_process(block, sim);
                        // Obtain results (statistics, etc) from processed block
                        lresults = results_new_from_photon_block(block);

                        #pragma omp critical(results)
                        {
                            // Add local results to shared results
                            results_add(results, lresults);
                            *num_photons_processed += num_photons;
                        }
                        #pragma omp critical(write)
                        {
                            // Write photons to simulation file
                            simulation_file_write_photons(file, block);
                            if(!*exit_request)
                                update_status(block->length,
                                              results->photons_received);
                        }

                        results_destroy(lresults);
                        photon_block_destroy(block);
                    }
                }
            }
        }
    }

    set_processing_time(results);
    free(sim->med_n_water_variables);
    free(sim->med_boundary_pos);
    free(sim->phase_layer_pos);
    for(int i = 0; i < sim->phase_layers; i++){
        for(int j = 0; j < sim->phase_max_x; j++){
            free(sim->phase_derivadasX[i][j]);
            free(sim->phase_derivadasY[i][j]);
        }
    }
    for(int i = 0; i < sim->phase_layers; i++){
        free(sim->phase_derivadasX[i]);
        free(sim->phase_derivadasY[i]);
    }
    free(sim->phase_derivadasX);
    free(sim->phase_derivadasY);

    if(local_malloc)
        free(num_photons_processed);
}


Results * get_results_from_file(Simulation * sim, SimulationFile * file,
                                 unsigned sweep_length)
{
    PhotonBlock * block; // Block of photons read
    Results * results;   // Results obtained from photons processed
    long num_photons;    // Number of photons in the simulation file
    unsigned k;          // Loop index

    // Initialize variables
    results = results_new_array(sweep_length);
    block = photon_block_new(sim->block_size);
    num_photons = simulation_file_get_num_photons(file);

    start_status_tracking(num_photons);
    do
    {
        // Read photons from simulation file
        simulation_file_read_photons(file, block);
        for(k = 0; k < sweep_length; ++k)
        {
            // Mark photons detected by the corresponding configuration
            photon_block_mark_received(block, &sim[k]);
            // Add results obtained only from photons marked as received
            results_add_photons_received(&results[k], block);
        }
        update_status(block->length, 0);
    }
    while(block->length == sim->block_size);
    end_status_tracking();

    photon_block_destroy(block);

    return results;
}