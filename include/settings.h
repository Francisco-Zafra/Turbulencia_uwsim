#ifndef __SETTINGS_H__
#define __SETTINGS_H__


#include <stdbool.h>
#include <stdint.h>
#include "json-glib/json-glib.h"


// Max length for name, family, and value.str fields
#define MAX_STRING_LENGTH   256u

// Options flags
#define JSON_HIDDEN         0x01
#define FILENAME_HIDDEN     0x02
#define NOT_PERSISTENT      0x08

// If this definition is changed, OPTION_NOT_FOUND value in settings.c must be
// reconsidered
typedef uint16_t settings_length_t;



typedef enum {tp_basic, tp_string, tp_double} OptionType;

typedef union
{
    double dbl;
    char str[MAX_STRING_LENGTH];
}
Value;


typedef struct
{
    OptionType type;
    unsigned id;
    char * name;
    char * family;
    int key;
    char * description;
    Value value;
    char * arg;
    char ** supported_args;
    unsigned index;
    bool fixed;
    bool json_printable;
    bool filename_included;
}
Option;


typedef struct
{
    bool active;
    Option * option;
    double start_value;
    double resolution;
    unsigned length;
}
Sweep;


typedef struct
{
    settings_length_t length;
    Option * options;
    Sweep sweep;
    char * current_family;
}
Settings;



Settings * settings_new();
Settings * settings_new_from(Settings * settings);
void settings_destroy(Settings * settings);

void settings_parse_input(Settings * settings, int argc, char * argv[]);

char * settings_get_filename(Settings * settings);

void settings_select_family(Settings * settings, char * family);

void settings_add_basic(Settings * settings,  char * name, int key,
                        char * description);
void settings_add_str(Settings * settings, char * name, int key,
                      char * description, unsigned id, unsigned flags,
                      char * value, char ** supported_values);
void settings_add_dbl(Settings * settings, char *name, int key,
                      char * description, unsigned id, unsigned flags,
                      double value);

void settings_set_str(Settings * settings, char * family, char * name,
                      char * value);
void settings_set_dbl(Settings * settings, char * family, char * name,
                      double value);

char * settings_get_str(Settings * settings, char * family, char * name);
unsigned settings_get_index(Settings * settings, char * family, char * name);
double * settings_get_dbl(Settings * settings, char * family, char * name);

bool settings_option_exists(Settings * settings, char * family, char * name);
Option * settings_get_option_by_id(Settings * settings, unsigned id);

bool settings_is_fixed(Settings * settings, char * family, char * name);
bool settings_is_storable(Settings * settings, unsigned id);
settings_length_t settings_get_number_of_options_fixed(Settings * settings);
settings_length_t settings_get_number_of_options_storable(Settings * settings);

bool settings_is_sweep_active(Settings * settings);
unsigned settings_get_sweep_length(Settings * settings);
double settings_get_sweep_resolution(Settings * settings);
void settings_apply_sweep(Settings * settings);
void settings_init_sweep(Settings * settings);
void settings_increase_sweep_start_value(Settings * settings,
                                         double start_value);
bool settings_is_sweep_option(Settings * settings, char * family, char * name);

JsonObject * settings_get_json_object(Settings * settings);


#endif