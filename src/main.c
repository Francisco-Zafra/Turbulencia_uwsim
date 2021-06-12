#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <error.h>

#include "configuration.h"
#include "utilities.h"
#include "status.h"
#include "processing.h"
#include "settings.h"
#include "simulation.h"
#include "simulation_file.h"

/* --- Signals management --------------------------------------------------- */

static bool exit_request = false;

void sigint_handler(int signal)
{
    exit_request = true;
    print_stop_warning(stderr);
}


/* --- Auxiliar functions --------------------------------------------------- */

void check_validity(Settings * settings)
{
    SourceType source_type;
    unsigned opts;

    if(!settings_is_fixed(settings, "surface", "depth"))
    {
        //unsigned opts = settings_get_number_of_options_fixed(settings);
        if(settings_is_fixed(settings, "surface", "sigma"))
            error(1, 0, "You cannot set standard deviation of the surface inclination if there isn't surface\n");
    }

    if(settings_is_fixed(settings, 0, "continue"))
    {
        // In this case 'photon' option is mandatory and 'quiet' is
        // optional
        opts = settings_get_number_of_options_fixed(settings);
        if(!settings_is_fixed(settings, 0, "photons"))
            error(1, 0,"'photons' option must be set with 'continue' option\n");
        if(opts > 3  ||
           (opts > 2  &&  !settings_is_fixed(settings, 0, "quiet")))
            error(1, 0, "Only 'photons' and 'quiet' options can be set with 'continue' option\n");
        if(settings_is_sweep_active(settings))
            error(1, 0, "Sweep is not allowed in 'continue' mode\n");
    }

    if(settings_is_fixed(settings, 0, "input file"))
    {
        // Persistent options (they can not be changed after simulation)
        if(settings_is_fixed(settings, 0, "photons")  ||
           settings_is_fixed(settings, "receptor", "x")  ||
           settings_is_fixed(settings, "receptor", "y")  ||
           settings_is_fixed(settings, "receptor", "z")  ||
           settings_is_fixed(settings, "receptor", "theta")  ||
           settings_is_fixed(settings, "receptor", "sigma")  ||
           settings_is_fixed(settings, "receptor", "phi")  ||
           settings_is_fixed(settings, 0, "rouletting")  ||
           settings_is_fixed(settings, "medium", "attenuation")  ||
           settings_is_fixed(settings, "medium", "albedo")  ||
           settings_is_fixed(settings, "source", "type")  ||
           settings_is_fixed(settings, "source", "divergence")  ||
           settings_is_fixed(settings, "source", "beam waist")  ||
           settings_is_fixed(settings, "source", "theta")  ||
           settings_is_fixed(settings, "source", "sigma")  ||
           settings_is_fixed(settings, "source", "phi")  ||
           settings_is_fixed(settings, "surface", "depth")  ||
           settings_is_fixed(settings, "medium", "index")  ||
           settings_is_fixed(settings, "receptor", "index")  ||
           settings_is_fixed(settings, "surface", "index")  ||
           settings_is_fixed(settings, "surface", "sigma"))
            error(1, 0, "Incompatible options\n");
    }
    else
    {
        if(settings_is_fixed(settings, 0, "time"))
           error(1, 0, "Time can not be set without input file\n");
            // TODO: make it possible (slower speed)
    }

    source_type = (SourceType)settings_get_index(settings, "source", "type");

    if(source_type == ideal)
    {
        if(settings_is_fixed(settings, "source", "divergence")  ||
        settings_is_fixed(settings, "source", "beam waist")  ||
        settings_is_fixed(settings, "source", "full angle"))
            error(1, 0, "An ideal source doesn't need 'divergence', 'beam waist' or 'full angle' to be set\n");
    }
    else if(source_type == lambertian)
    {
        if(settings_is_fixed(settings, "source", "divergence")  ||
        settings_is_fixed(settings, "source", "beam waist")  ||
        !settings_is_fixed(settings, "source", "full angle"))
            error(1, 0, "A lambertian source only needs 'full angle' to be set\n");
    }
    else if(source_type == gaussian)
    {
        if(!settings_is_fixed(settings, "source", "divergence")  ||
        !settings_is_fixed(settings, "source", "beam waist")  ||
        settings_is_fixed(settings, "source", "full angle"))
            error(1, 0, "A gaussian source only needs divergence and beam waist to be set\n");
    }
}


void check_compatibility(Settings * settings, Simulation * sim,
                         SimulationFile * file)
{
    Settings * father_settings;
    Simulation * father_sim;
    unsigned sweep_length;
    unsigned k;

    father_settings = settings_new_from(settings);
    simulation_file_read_all_settings(file, father_settings);
    father_sim = simulation_new_from_settings_start(father_settings);
    sweep_length = settings_get_sweep_length(settings);

    for(k = 0; k < sweep_length; ++k)
        if(!simulation_is_compatible(&sim[k], father_sim))
            error(1, 0, "Incompatible values\n");

    settings_destroy(father_settings);
    simulation_destroy(father_sim);
}


/* --- Main ----------------------------------------------------------------- */

int main(int argc , char* argv[])
{
    Settings * settings;
    Simulation * sim;
    Results * results;
    SimulationFile * sim_file;
    FILE * json_file;
    char * filename;
    char * input_filename;
    long num_photons_processed;
    unsigned sweep_length;
    unsigned k;
    bool use_storage;
    bool readable_mode;

    settings = settings_new();
    num_photons_processed = 0;
    sim = NULL;
    results = NULL;
    filename = NULL;

    // Default values
    char * SUP_p[] = {"hg", "tthg", 0};
    char * SUP_s[] = {"ideal", "isotropic", "gaussian", "lambertian", 0};
    char * SUP_M[] = {"all", "los", "nlos", 0};

    settings_select_family(settings, "medium");
    settings_add_str(settings, "phase function", 'P', DOC_P, ID_P, 0, "hg",
                     SUP_p);

    settings_add_dbl(settings, "layers", 'L', DOC_L, ID_L, 0, 1.0);
    settings_add_dbl(settings, "varZ", 705, DOC_705, ID_705, 0, 1);

    settings_add_dbl(settings, "var_n_water", 706, DOC_706, ID_706, 0, 0.05);
    settings_add_dbl(settings, "boundary_max_theta", 707, DOC_707, ID_707, 0, 0.5);

    settings_add_dbl(settings, "g", 'g', DOC_g, ID_g, 0, 0.924);
    settings_add_dbl(settings, "attenuation", 'c', DOC_c, ID_c, 0, 2.19);
    settings_add_dbl(settings, "albedo", 'w', DOC_w, ID_w, 0, 0.83);
    settings_add_dbl(settings, "index", 700, DOC_700, ID_700,
                     FILENAME_HIDDEN, 1.33);
    
    settings_add_str(settings, "phase_json", 710, DOC_710, ID_710, 0, "", NULL);

    settings_select_family(settings, "surface");
    settings_add_dbl(settings, "depth", 'D', DOC_D, ID_D, 0, INFINITY);
    settings_add_dbl(settings, "sigma", 'W', DOC_W, ID_W, 0, 0.0);
    settings_add_dbl(settings, "index", 702, DOC_702, ID_702,
                     FILENAME_HIDDEN, 1.0);
    
    settings_select_family(settings, "floor");
    settings_add_dbl(settings, "depth", 708, DOC_708, ID_708, 0, INFINITY);
    settings_add_dbl(settings, "index", 709, DOC_709, ID_709,
                     FILENAME_HIDDEN, 1.5);

    settings_select_family(settings, "source");
    settings_add_str(settings, "type", 'T', DOC_T, ID_T, 0, "ideal", SUP_s);
    settings_add_dbl(settings, "full angle", 'F', DOC_F, ID_F, FILENAME_HIDDEN,
                     INFINITY);
    settings_add_dbl(settings, "divergence", 'd', DOC_d, ID_d, FILENAME_HIDDEN,
                     INFINITY);
    settings_add_dbl(settings, "beam waist", 'b', DOC_b, ID_b, FILENAME_HIDDEN,
                     INFINITY);
    settings_add_dbl(settings, "theta", 'E', DOC_E, ID_E, 0, 0.0);
    settings_add_dbl(settings, "sigma", 'S', DOC_S, ID_S, 0, 0.0);
    settings_add_dbl(settings, "phi", 'H', DOC_H, ID_H, 0, 0.0);

    settings_select_family(settings, "receptor");
    settings_add_dbl(settings, "x", 'x', DOC_x, ID_x, 0, 0.0);
    settings_add_dbl(settings, "y", 'y', DOC_y, ID_y, 0, 0.0);
    settings_add_dbl(settings, "z", 'z', DOC_z, ID_z, 0, 1.0);
    settings_add_dbl(settings, "x offset", 'X', DOC_X, ID_X, 0, 0.0);
    settings_add_dbl(settings, "y offset", 'Y', DOC_Y, ID_Y, 0, 0.0);
    settings_add_dbl(settings, "aperture", 'a', DOC_a, ID_a, 0, INF);
    settings_add_dbl(settings, "fov", 'f', DOC_f, ID_f, 0, 180.0);
    settings_add_dbl(settings, "index", 703, DOC_703, ID_703,
                     FILENAME_HIDDEN, 1.0);
    settings_add_dbl(settings, "theta", 'e', DOC_e, ID_e, 0, 0.0);
    settings_add_dbl(settings, "sigma", 's', DOC_s, ID_s, 0, 0.0);
    settings_add_dbl(settings, "phi", 'h', DOC_h, ID_h, 0, 0.0);

    settings_select_family(settings, 0);
    settings_add_dbl(settings, "time", 't', DOC_t, 0, 0, INFINITY);
    settings_add_dbl(settings, "photons", 'n', DOC_n, ID_n, 0, 0.0);
    settings_add_dbl(settings, "block size", 704, DOC_704, 0, FILENAME_HIDDEN,
                     1e5);
    settings_add_dbl(settings, "length factor", 'l', DOC_l, ID_l,
                     FILENAME_HIDDEN, 5);
    settings_add_dbl(settings, "max events", 'm', DOC_m, ID_m, FILENAME_HIDDEN,
                     255);
    settings_add_str(settings, "mode", 'M', DOC_M, ID_M, 0, "all", SUP_M);
    settings_add_dbl(settings, "rouletting", 701, DOC_701, ID_701,
                     FILENAME_HIDDEN, 0);
    settings_add_str(settings, "input file", 'i', DOC_i, 0,
                     FILENAME_HIDDEN, "", 0);
    settings_add_basic(settings, "readable", 'r', DOC_r);
    settings_add_basic(settings, "quiet", 'q', DOC_q);
    settings_add_str(settings, "continue", 'C', DOC_C, 0,
                     FILENAME_HIDDEN | JSON_HIDDEN, "", 0);
    settings_add_str(settings, "output", 'o', DOC_o, 0,
                     FILENAME_HIDDEN | JSON_HIDDEN, "", 0);

    settings_parse_input(settings, argc, argv);
    check_validity(settings);
    sweep_length = settings_get_sweep_length(settings);

    if(settings_is_fixed(settings, 0, "quiet"))
        stdout = freopen("/dev/null", "a", stdout);

    if(settings_is_fixed(settings, 0, "output"))
        filename = settings_get_str(settings, 0, "output");

    if(settings_is_fixed(settings, 0, "input file"))
    {
        input_filename = settings_get_str(settings, 0, "input file");
        if(!file_exists(input_filename))
            error(1, 0, "%s doesn't exists\n", input_filename);
        sim_file = simulation_file_new_binary(input_filename);
        simulation_file_read_not_fixed_settings(sim_file, settings);
        update_time_option_if_required(settings, sim_file);
        sim = simulation_new_array_from_settings(settings);
        // We have to ensure that, for example, FOV used in the simulation file is greater than the required one
        check_compatibility(settings, sim, sim_file);
        print_settings(stdout, settings);

        results = get_results_from_file(sim, sim_file, sweep_length);
    }
    else if(*settings_get_dbl(settings, 0, "photons") > 0.0)
    {
        // if sim->photons == 0 there is nothing to simulate
        results = results_new_array(sweep_length);
        use_storage = sweep_length == 1;
        readable_mode = settings_is_fixed(settings, 0, "readable");

        if(use_storage)
        {
            if(settings_is_fixed(settings, 0, "continue"))
            {
                signal(SIGINT, sigint_handler);
                filename = settings_get_str(settings, 0, "continue");
                if(!file_exists(filename))
                    error(1, 0, "%s doesn't exists\n", filename);
                sim_file = simulation_file_new_binary(filename);
                simulation_file_read_not_fixed_settings(sim_file, settings);
                num_photons_processed = (long)simulation_file_get_dbl(sim_file,
                                                                      ID_n);
            }
            else if(readable_mode)
            {
                if(filename == NULL)
                    filename = get_filename(settings, RD_EXTENSION);
                sim_file = simulation_file_new_readable(filename);
            }
            else
            {
                signal(SIGINT, sigint_handler);
                if(filename == NULL)
                    filename = get_filename(settings, BN_EXTENSION);
                sim_file = simulation_file_new_binary(filename);
                simulation_file_write_settings(sim_file, settings);
            }
        }
        else
            sim_file = NULL;

        sim = simulation_new_array_from_settings(settings);
        print_settings(stdout, settings);

        start_status_tracking(sim->photons * sweep_length -
                              num_photons_processed);
        run_simulation(sim, results, sim_file, &num_photons_processed,
                               &exit_request);
        for(k = 1; k < sweep_length; ++k)
            run_simulation(&sim[k], &results[k], sim_file, NULL, &exit_request);
        end_status_tracking();

        // Overwrite real number of photons processed
        if(use_storage  &&  !readable_mode)
        {
            settings_set_dbl(settings, 0, "photons",
                             (double)num_photons_processed);
            simulation_file_overwrite_setting(sim_file, settings, ID_n);
            simulation_file_destroy(sim_file);
        }
    }
    else
    {
        print_only_settings(stdout, settings);
    }

    if(sim != NULL)
    {
        for(k = 0; k < sweep_length; ++k)
            results_set_events_threshold(&results[k], sim[k].events_threshold);
        simulation_destroy_array(sim);
    }

    if(results != NULL)
    {
        if(settings_is_sweep_active(settings))
        {
            if(filename == NULL)
                filename = get_filename(settings, SW_EXTENSION);
            json_file = fopen(filename, "w");
            print_settings_and_results(json_file, settings, results);
            fclose(json_file);
        }

        print_results(stdout, results, sweep_length);

        results_destroy_array(results, sweep_length);
    }

    settings_destroy(settings);

    return 0;
}