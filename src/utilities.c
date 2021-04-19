#include "utilities.h"

#include <stdbool.h>
#include <math.h>
#include <string.h>	// strchr()
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>	// isprint()
#include "json-glib/json-glib.h"


#include "settings.h"
#include "simulation_file.h"



char * get_filename_extension(char *filename)
{
    char * dot = strrchr(filename, '.');
    if(dot == NULL || dot == filename)
        return "";
    return dot + 1;
}


bool file_exists(char * filename)
{
    return access(filename, F_OK) == 0;
}


char * get_filename(Settings * settings, char * extension)
{
    char * filename;
    unsigned index;

    filename = malloc(255);

    index = 0;

    do
    {
        sprintf(filename, "%s_%02u%s", settings_get_filename(settings),
                index, extension);
        index++;
    }
    while(file_exists(filename));

    return filename;
}


char * get_string_from_json(JsonObject * jobject)
{
    JsonNode * jnode;

    jnode = json_node_new(JSON_NODE_OBJECT);
    jnode = json_node_init_object(jnode, jobject);

    return json_to_string(jnode, true);
}


void print_settings_and_results(FILE * stream, Settings * settings,
                              Results * results)
{
    JsonObject * jall;
    JsonObject * jsettings;
    JsonObject * jresults;
    char * string;
    unsigned sweep_length;

    sweep_length = settings_get_sweep_length(settings);
    jsettings = settings_get_json_object(settings);
    jresults = results_get_json_object(results, sweep_length);
    jall = json_object_new();
    json_object_set_object_member(jall, "settings", jsettings);
    json_object_set_object_member(jall, "results", jresults);
    string = get_string_from_json(jall);
    fprintf(stream, "%s\n", string);

    free(string);
}


void print_only_settings(FILE * stream, Settings * settings)
{
    JsonObject * jall;
    JsonObject * jsettings;
    char * string;

    jsettings = settings_get_json_object(settings);
    jall = json_object_new();
    json_object_set_object_member(jall, "settings", jsettings);
    string = get_string_from_json(jall);
    fprintf(stream, "%s\n", string);

    free(string);
}


void print_settings(FILE * stream, Settings * settings)
{
    JsonObject * jall;
    JsonObject * jsettings;
    char * string;

    jsettings = settings_get_json_object(settings);
    jall = json_object_new();
    json_object_set_object_member(jall, "settings", jsettings);
    string = get_string_from_json(jall);
    string[strlen(string) - 1] = 0;
    string[strlen(string) - 1] = ',';
    fprintf(stream, "%s\n", string);

    free(string);
}


void print_results(FILE * stream, Results * results, unsigned sweep_length)
{
    JsonObject * jall;
    JsonObject * jresults;
    char * string;
    char * cut_string;

    jresults = results_get_json_object(results, sweep_length);
    jall = json_object_new();
    json_object_set_object_member(jall, "results", jresults);
    string = get_string_from_json(jall);
    cut_string = string + 2;
    fprintf(stream, "%s\n", cut_string);

    free(string);
}


void print_stop_warning(FILE * file)
{
    fprintf(file, "\033[1;35m");   // Set output color to green
    fprintf(file, "\nStopping simulation safely\n");
    fprintf(file, "\033[0m");   // Reset output color to default
}


void update_time_option_if_required(Settings * settings, SimulationFile * file)
{
    reception_mode_t mode;
    float ballistic_tof;
    float time;

    if(settings_is_fixed(settings, 0, "time"))
    {
        mode = (reception_mode_t)settings_get_index(settings, 0, "mode");
        ballistic_tof = simulation_file_get_minimun_time(file, mode);
        if(settings_is_sweep_option(settings, 0, "time"))
        {
            settings_increase_sweep_start_value(settings, ballistic_tof);
        }
        else
        {
            time = *settings_get_dbl(settings, 0, "time");
            time += ballistic_tof;
            settings_set_dbl(settings, 0, "time", time);
        }
    }
}