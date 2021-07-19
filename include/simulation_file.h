#ifndef __SIMULATION_FILE_H__
#define __SIMULATION_FILE_H__

#include <stdio.h>
#include <stdbool.h>

#include "settings.h"
#include "photon_block.h"

typedef enum
{
    closed_mode,
    read_mode,
    write_mode,
    rdwr_mode,
    append_mode,
    reading_photons_mode
}
file_mode_t;

typedef struct
{
    char * filename;
    FILE * stream;
    bool binary_format;
    mode_t mode;
}
SimulationFile;





SimulationFile * simulation_file_new_binary(char * filename);
SimulationFile * simulation_file_new_readable(char * filename);
void simulation_file_write_photons(SimulationFile * file, PhotonBlock * block);
void simulation_file_read_photons(SimulationFile * file, PhotonBlock * block);
void simulation_file_write_settings(SimulationFile * file, Settings * set);
void simulation_file_overwrite_setting(SimulationFile * file, Settings  *set,
                                       unsigned id);
void simulation_file_read_all_settings(SimulationFile * file, Settings * set);
void simulation_file_read_not_fixed_settings(SimulationFile * file,
                                             Settings * set);
long simulation_file_get_num_photons(SimulationFile * file);
float simulation_file_get_minimun_time(SimulationFile * file,
                                       reception_mode_t mode);
double simulation_file_get_dbl(SimulationFile * file, unsigned id);
void simulation_file_destroy(SimulationFile * sim_file);


#endif