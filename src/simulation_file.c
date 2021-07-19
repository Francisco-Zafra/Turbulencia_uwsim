#include "simulation_file.h"

#include <errno.h>
#include <error.h>
#include <stdint.h>
#include <math.h>   // INFINITY

#include "photon_block.h"


/*
    File structure: header - settings - body
*/


#define SIMULATION_FILE_VERSION 2
#define SIMULATION_FILE_ID 0xF603C984


typedef uint16_t option_id_t;

typedef struct
{
    option_id_t id;
    Value value;
}
StoredOption;

#pragma pack(push) // Save the current alignment configuration
#pragma pack(1) // To avoid padding
typedef struct
{
    float weight;
    float tof;
    float distance;
    float x, y;
    float theta;
    uint8_t events;
    bool reflected;
}
StoredPhoton;
#pragma pack(pop) // Restore the alignment configuration


typedef struct
{
    uint32_t id;
    uint16_t version;
    settings_length_t length;
    uint64_t reserved;  // For future changes
}
Header;


// This function acts exactly like fread, only that if the bytes are not all
// read, the function exits reporting the error
void freadn(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    unsigned num_read;

    num_read = fread(ptr, size, nmemb, stream);
    if (num_read < nmemb)
        error(1, errno, "Unable to read from file %p\n", stream);
}


void verify_binary_format(SimulationFile * file)
{
    if(!file->binary_format)
        error(1, 0, "Trying to read or write in binary mode from text file\n");
}


// stream must be en read mode
settings_length_t read_length(FILE * stream)
{
    Header  header;

    fseek(stream, 0, SEEK_SET);
    freadn(&header, sizeof(header), 1, stream);

    if(header.id != SIMULATION_FILE_ID)
        error(1, 0, "Unrecognized file\n");
    if(header.version != SIMULATION_FILE_VERSION)
        error(1, 0, "Version %u not supported\n", header.version);

    return header.length;
}


void write_header(FILE * stream, settings_length_t length)
{
    Header header;
    header.id = SIMULATION_FILE_ID;
    header.version = SIMULATION_FILE_VERSION;
    header.length = length;
    fwrite(&header, sizeof(header), 1, stream);
}


void move_to_settings_start(FILE * stream)
{
    fseek(stream, sizeof(Header), SEEK_SET);
}


void move_to_photons_start(FILE * stream)
{
    settings_length_t length;

    length = read_length(stream);
    move_to_settings_start(stream);
    fseek(stream, length * sizeof(StoredOption), SEEK_CUR);
}


void ensure_mode(SimulationFile * file, file_mode_t mode)
{
    if(file->mode != mode)
    {
        if(file->mode != closed_mode)
            if(file->stream != NULL)
                fclose(file->stream);
        file->mode = mode;
        switch(mode)
        {
            case closed_mode:
                break;
            case read_mode:
                file->stream = fopen(file->filename, "r");
                break;
            case write_mode:
                file->stream = fopen(file->filename, "w");
                break;
            case rdwr_mode:
                file->stream = fopen(file->filename, "r+");
                break;
            case append_mode:
                file->stream = fopen(file->filename, "a");
                break;
            case reading_photons_mode:
                file->stream = fopen(file->filename, "r");
                move_to_photons_start(file->stream);
                break;
            default:
                error(1, 0, "file mode not recognized\n");
        }
        if(file == NULL)
            error(1, errno, "Unable to open '%s'\n", file->filename);
    }
}


// This function assumes that stream is in read mode
void move_to_value_position(FILE * stream, unsigned id)
{
    StoredOption stored_opt;
    settings_length_t length;
    settings_length_t k;
    bool found;

    k = 0;
    found = false;

    length = read_length(stream);
    move_to_settings_start(stream);

    while(!found  &&  k < length)
    {
        freadn(&stored_opt, sizeof(stored_opt), 1, stream);
        if(stored_opt.id == id)
            found = true;
        ++k;
    }

    if(!found)
        error(1, 0, "Option %u not found\n", id);

    fseek(stream, -sizeof(Value), SEEK_CUR);
}


void read_conditionally(SimulationFile * file, Settings * set, bool read_fixed)
{
    settings_length_t length;
    settings_length_t k;
    Option * option;
    StoredOption opt;

    verify_binary_format(file);
    ensure_mode(file, read_mode);
    length = read_length(file->stream);
    move_to_settings_start(file->stream);
    for(k = 0; k < length; ++k)
    {
        freadn(&opt, sizeof(opt), 1, file->stream);
        option = settings_get_option_by_id(set, opt.id);
        if(option != NULL)
        {
            if(read_fixed  ||  (!read_fixed  &&  !option->fixed))
                memcpy(&option->value, &opt.value, sizeof(option->value));
        }
        else
        {
            error(1, 0, "option in file %u not allowed\n", opt.id);
        }
    }
}


void simulation_file_write_settings(SimulationFile * file, Settings * set)
{
    StoredOption stored_opt;
    settings_length_t k;

    verify_binary_format(file);
    ensure_mode(file, write_mode);
    write_header(file->stream, settings_get_number_of_options_storable(set));
    for(k = 0; k < set->length; ++k)
    {
        if(settings_is_storable(set, set->options[k].id))
        {
            stored_opt.id = set->options[k].id;
            memcpy(&stored_opt.value, &set->options[k].value,
                sizeof(stored_opt.value));
            fwrite(&stored_opt, sizeof(stored_opt), 1, file->stream);
        }
    }
}


void simulation_file_overwrite_setting(SimulationFile * file, Settings  *set,
                                       unsigned id)
{
    Option * option;

    if(file == NULL)
        return;

    option = settings_get_option_by_id(set, id);

    verify_binary_format(file);
    ensure_mode(file, rdwr_mode);
    move_to_value_position(file->stream, id);
    fwrite(&option->value, sizeof(option->value), 1, file->stream);
}


void simulation_file_read_all_settings(SimulationFile * file, Settings * set)
{
    read_conditionally(file, set, true);
}


void simulation_file_read_not_fixed_settings(SimulationFile * file,
                                             Settings * set)
{
    read_conditionally(file, set, false);
}


SimulationFile * simulation_file_new(char * filename, bool is_binary)
{
    SimulationFile * file;

    file = malloc(sizeof(SimulationFile));

    file->binary_format = is_binary;
    file->filename = strdup(filename);
    file->stream = NULL;
    file->mode = closed_mode;

    return file;
}


SimulationFile * simulation_file_new_binary(char * filename)
{
    return simulation_file_new(filename, true);
}


SimulationFile * simulation_file_new_readable(char * filename)
{
    return simulation_file_new(filename, false);
}


void simulation_file_write_photons(SimulationFile * file,
                                   PhotonBlock * photon_block)
{
    Photon * photons;
    long k;

    if(file == NULL)
        return;

    photons = photon_block->photons;
    ensure_mode(file, append_mode);

    if(file->binary_format)
    {
        for(k = 0; k < photon_block->length; ++k)
            if(photons[k].received)
                fwrite(&photons[k], sizeof(StoredPhoton), 1, file->stream);
    }
    else
    {
        for(k = 0; k < photon_block->length; ++k)
            if(photons[k].received)
                fprintf(file->stream, "%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",
                        photons[k].tof,
                        photons[k].weight, photons[k].x, photons[k].y,
                        photons[k].theta, photons[k].distance,
                        photons[k].events, photons[k].reflected);
    }
}


// In success, this function return true    // TODO: use this function
bool read_photon(SimulationFile * file, Photon * photon)
{
    return (bool)fread(photon, sizeof(StoredPhoton), 1, file->stream);
}


void simulation_file_read_photons(SimulationFile * file, PhotonBlock * block)
{
    size_t num_read;
    long k;

    k = 0;
    ensure_mode(file, reading_photons_mode);

    do
    {
        num_read = fread(&block->photons[k], sizeof(StoredPhoton), 1,
                         file->stream);
        ++k;
    }
    while(k < block->length  &&  num_read != 0);

    if(k < block->length)
        block->length = k - 1;
}


void simulation_file_destroy(SimulationFile * sim_file)
{
    fclose(sim_file->stream);
    free(sim_file->filename);
    free(sim_file);
}


double simulation_file_get_dbl(SimulationFile * file, unsigned id)
{
    double value;

    move_to_value_position(file->stream, id);
    freadn(&value, sizeof(double), 1, file->stream);

    return value;
}


long simulation_file_get_num_photons(SimulationFile * file)
{
    long num_photons;
    long start, end;

    num_photons = 0;
    ensure_mode(file, read_mode);

    move_to_photons_start(file->stream);
    start = ftell(file->stream);
    fseek(file->stream, 0, SEEK_END);   // Move to end of file
    end = ftell(file->stream);
    num_photons = (end - start) / sizeof(StoredPhoton);

    return num_photons;
}


float simulation_file_get_minimun_time(SimulationFile * file,
                                       reception_mode_t mode)
{
    StoredPhoton photon;
    size_t num_read;
    float minimun_time;

    minimun_time = INFINITY;
    ensure_mode(file, reading_photons_mode);

    // Find the minimun travel_time value
    num_read = fread(&photon, sizeof(StoredPhoton), 1, file->stream);
    while(num_read != 0)
    {
        if(photon.tof < minimun_time  &&  (mode == all  ||
           (mode == los  &&  !photon.reflected)  ||
           (mode == nlos  &&  photon.reflected)))
            minimun_time = photon.tof;
        num_read = fread(&photon, sizeof(StoredPhoton), 1, file->stream);
    }

    ensure_mode(file, closed_mode);

    return minimun_time;
}