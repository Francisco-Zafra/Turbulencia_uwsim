#ifndef __DATA_MANAGEMENT_H__
#define __DATA_MANAGEMENT_H__


#include <stdbool.h>
#include <stdio.h>

#include "configuration.h"
#include "settings.h"
#include "results.h"
#include "simulation.h"
#include "simulation_file.h"


bool file_exists(char * filename);
char * get_filename(Settings * settings, char * extension);
char * get_filename_extension(char *filename);

void print_settings_and_results(FILE * stream, Settings * settings,
                              Results * results);
void print_only_settings(FILE * stream, Settings * settings);
void print_settings(FILE * stream, Settings * settings);
void print_results(FILE * stream, Results * results, unsigned sweep_length);
void print_stop_warning(FILE * file);

void update_time_option_if_required(Settings * settings,SimulationFile * file);


#endif