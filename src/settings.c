#include "settings.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <argp.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>

#include "configuration.h"


/* --- Local defines -------------------------------------------------------- */

#define GLOBAL_FAMILY       ""   // String representing Global family setting
#define MAX_HELP_LENGTH     1024


const char *argp_program_version = PROGRAM_VERSION;
const char *argp_program_bug_address = PROGRAM_BUG_ADDRESS;


/* --- Parse auxiliar functions --------------------------------------------- */

bool is_sweep_string(char * sweep_str)
{
    char * temp;
    temp = strchr(sweep_str, ':');
    if(temp == NULL)
        return false;
    else
    {
        temp++;
        temp = strchr(temp, ':');
    }
    if(temp == NULL)
        return false;
    else
    {
        temp++;
        temp = strchr(temp, ':');
    }
    if(temp == NULL)
        return true;
    else
        return false;
}


void get_sweep_parameters(char * argument, double * start_value,
                          double * resolution, unsigned * length)
{
    char * start_value_str;
    char * resolution_str;
    char * end_value_str;
    double end_value;

    start_value_str = strdup(argument);
    resolution_str = strchr(argument, ':');
    *resolution_str = '\0';
    resolution_str++;
    end_value_str = strchr(resolution_str, ':');
    *end_value_str = '\0';
    end_value_str++;

    *start_value = atof(start_value_str);
    *resolution = atof(resolution_str);
    end_value = atof(end_value_str);
    *length = (unsigned)floor((end_value - *start_value) / *resolution) + 1;
}


// index is set to the string position within de list if is found
bool string_is_in_list(char * string, char * list[], unsigned * index)
{
    bool found;
    unsigned k;

    found = false;
    k = 0;

    if(list == NULL)
        found = true;
    else
        while(!found  &&  list[k] != NULL)
        {
            found = (strcmp(string, list[k]) == 0);
            ++k;
        }

    if(list != NULL && found)
        *index = k - 1;

    return found;
}


Option * settings_get_option_by_key(Settings * settings, int key)
{
    settings_length_t k;
    bool found;

    k = 0;
    found = false;

    while(!found  &&  k < settings->length)
    {
        if(settings->options[k].key == key)
            found = true;
        else
            ++k;
    }

    if(found == false)
        return NULL;
    else
        return &(settings->options[k]);
}


static int parse_opt(int key, char *arg, struct argp_state *state)
{
    Settings * settings;
    Option * option;
    Sweep * sweep;

    settings = state->input;
    option = settings_get_option_by_key(settings, key);
    if(option == NULL)
        return ARGP_ERR_UNKNOWN;

    option->fixed = true;
    if(arg != NULL)
        option->arg = strdup(arg);
    if(option->type == tp_string)
    {
        if(strlen(arg) > MAX_STRING_LENGTH + 1)
            error(1, 0, "Argument too long: %s", arg);
        if(string_is_in_list(arg, option->supported_args, &option->index))
            strcpy(option->value.str, arg);
        else
            error(1, 0, "Invalid argument '%s'\n", arg);

    }
    else if(option->type == tp_double)
    {
        option->value.dbl = atof(arg);

        if(is_sweep_string(arg))
        {
            sweep = &(settings->sweep);
            if(sweep->active)
                error(1, 0, "Only one sweep is allowed\n");

            sweep->active = true;
            sweep->option = option;
            get_sweep_parameters(arg, &sweep->start_value, &sweep->resolution,
                                 &sweep->length);
        }
    }

    return 0;
}


char * help_filter (int key, const char *text, void *input)
{
    Settings * settings;
    Option * option;
    char * new_text;
    unsigned k;

    settings = (Settings *)input;
    option = settings_get_option_by_key(settings, key);
    k = 0;

    if(option != NULL  &&  option->type == tp_string  &&
       option->supported_args != NULL)
    {
        new_text = malloc(MAX_HELP_LENGTH);

        sprintf(new_text, "%s (Supported values are:", text);
        // TODO: use asprintf

        while(option->supported_args[k] != 0)
        {
            if(k == 0)
                sprintf(new_text, "%s '%s'", new_text,
                        option->supported_args[k]);
            else
                sprintf(new_text, "%s, '%s'", new_text,
                        option->supported_args[k]);
            ++k;
        }

        sprintf(new_text, "%s)", new_text);

        return new_text;
    }
    else
    {
        return (char *)text;
    }
}


bool is_only_name(Settings * settings, char * name)
{
    unsigned k;
    unsigned count;

    count = 0u;

    for(k = 0; k < settings->length; ++k)
    {
        if(strcmp(name,settings->options[k].name) == 0)
            ++count;
    }

    if(count == 0)
        error(1, 0, "name %s not found\n", name);

    return count == 1;
}


char * option_get_printable_name(Option * option, Settings * settings)
{
    char * whitespace;
    char * family;
    char * name;
    char * printable_name;

    family = strdup(option->family);
    name = strdup(option->name);


    // White spaces are replaced by '-'
    while((whitespace = strchr(name, ' ')) != NULL)
        *whitespace = '-';
    if(strcmp(family, GLOBAL_FAMILY) != 0  &&
       !is_only_name(settings, option->name))
    {
        while((whitespace = strchr(family, ' ')) != NULL)
            *whitespace = '-';
        printable_name = malloc(strlen(family) + strlen(name) + 1);
        sprintf(printable_name, "%s-%s", family, name);
    }
    else
    {
        printable_name = strdup(name);
    }

    free(family);
    free(name);

    return printable_name;
}


/* --- Option auxiliar functions -------------------------------------------- */

Option * settings_get_option_or_null(Settings * settings, char * family,
                                     char * name)
{
    unsigned k;
    bool found;

    found = false;
    k = 0u;

    if(family == NULL)
        family = strdup(GLOBAL_FAMILY);

    while(k < settings->length  &&  !found)
    {
        if(strcmp(family, settings->options[k].family) == 0  &&
            strcmp(name, settings->options[k].name) == 0)
            found = true;
        ++k;
    }

    if(!found)
        return NULL;
    else
        return &(settings->options[k - 1]);
}


Option * settings_get_option(Settings * settings, char * family, char * name)
{
    Option * option;

    option = settings_get_option_or_null(settings, family, name);

    if(option == NULL)
        error(1, 0, "Option '%s' - '%s' doesn's exists\n", family, name);

    return option;
}


Option * settings_get_option_by_id(Settings * settings, unsigned id)
{
    Option * option;
    settings_length_t k;
    bool found;

    option = NULL;
    k = 0;
    found = false;

    while(!found  &&  k < settings->length)
    {
        if(settings->options[k].id == id)
        {
            option = &settings->options[k];
            found = true;
        }

        ++k;
    }

    return option;
}


bool settings_option_exists(Settings * settings, char * family, char * name)
{
    return settings_get_option_or_null(settings, family, name) != NULL;
}


Option * settings_add_option(Settings * settings, char *name, int key,
                             char * description, unsigned id, unsigned flags,
                             OptionType type)
{
    Option * option;
    unsigned index;

    if(settings_option_exists(settings, settings->current_family, name))
        error(1, 0, "Option with family '%s' and name '%s' already exists",
              settings->current_family, name);

    index = settings->length;
    settings->length++;
    settings->options = realloc(settings->options,
                                settings->length * sizeof(Option));
    option = &(settings->options[index]);
    option->type = type;
    option->id = id;
    option->name = strdup(name);
    option->family = strdup(settings->current_family);
    option->key = key;
    option->description = strdup(description);
    option->arg = NULL;
    option->index = 0;
    option->fixed = false;

    option->json_printable = true;
    option->filename_included = true;
    if((flags & JSON_HIDDEN) != 0)
        option->json_printable = false;
    if((flags & FILENAME_HIDDEN) != 0)
        option->filename_included = false;

    if(!isprint(key))
        option->filename_included = false;

    return option;
}


/* --- Public functions ----------------------------------------------------- */

Settings * settings_new()
{
    Settings * settings;

    settings = malloc(sizeof(Settings));
    settings->length = 0;
    settings->options = NULL;
    settings->sweep.active = false;
    settings->sweep.length = 1;
    settings->sweep.resolution = 0.0;
    settings->current_family = NULL;

    return settings;
}


Settings * settings_new_from(Settings * settings)
{
    Settings * settings_new;
    settings_length_t k;

    settings_new = malloc(sizeof(Settings));

    settings_new->length = settings->length;
    settings_new->current_family = NULL;
    settings_new->sweep = settings->sweep;
    settings_new->options = malloc(settings->length * sizeof(Option));
    for(k = 0; k < settings->length; ++k)
        settings_new->options[k] = settings->options[k];

    return settings_new;
}


void settings_destroy(Settings * settings)
{
    if(settings->options != NULL)
        free(settings->options);
    if(settings->current_family != NULL)
        free(settings->current_family);
    free(settings);
}


void settings_parse_input(Settings * settings, int argc, char * argv[])
{
    struct argp_option * options;
    settings_length_t k;

    options = malloc(sizeof(struct argp_option) * (settings->length + 1));

    for(k = 0; k < settings->length; ++k)
    {
        options[k].name = option_get_printable_name(&settings->options[k],
                                                    settings);
        options[k].key = settings->options[k].key;
        options[k].doc = strdup(settings->options[k].description);
        options[k].flags = 0;
        options[k].group = 0;
        if(settings->options[k].type == tp_string)
            options[k].arg = strdup("<arg>");
        else if(settings->options[k].type == tp_double)
            options[k].arg = strdup("<value>");
        else if(settings->options[k].type == tp_basic)
            options[k].arg = 0;
        else
            error(1, 0, "Option %i not implemented", settings->options[k].type);

    }
    memset(&options[k], 0, sizeof(struct argp_option)); // To mark the array end

    struct argp argp = {options, parse_opt, 0, 0, 0, help_filter};
    argp_parse(&argp, argc, argv, 0, 0, settings);
}


bool option_is_filename_printable(Option * option)
{
    return (option->type == tp_string  ||  option->type == tp_double)  &&
           isprint(option->key);
}


char * settings_get_filename(Settings * settings)
{
    Option * option;
    char * filename;
    char * filename_dup;
    char * value;
    settings_length_t k;
    bool first;

    filename = malloc(MAX_STRING_LENGTH);
    value = malloc(MAX_STRING_LENGTH);
    *filename = '\0';
    first = true;

    for(k = 0; k < settings->length; ++k)
    {
        option = &settings->options[k];
        if(!option->filename_included  ||  !isprint(option->key))
            continue;

        if(option->type == tp_double)
        {
            if(settings->sweep.active  &&  settings->sweep.option == option)
                strcpy(value, option->arg);
            else
                sprintf(value, "%.3g", option->value.dbl);
        }
        else if(option->type == tp_string)
        {
            strcpy(value, option->value.str);
        }

        if(first)
        {
           sprintf(filename, "%c%s", option->key, value);
           first = false;
        }
        else
        {
            filename_dup = strdup(filename);
            sprintf(filename, "%s_%c%s", filename_dup, option->key, value);
            free(filename_dup);
        }
    }

    free(value);

    return filename;
}


void settings_select_family(Settings * settings, char * family)
{
    if(settings->current_family != NULL)
        free(settings->current_family);

    if(family == NULL)
        settings->current_family = strdup(GLOBAL_FAMILY);
    else
        settings->current_family = strdup(family);
}


void settings_add_basic(Settings * settings,  char * name, int key,
                        char * description)
{
    settings_add_option(settings, name, key, description, 0,
                        JSON_HIDDEN | FILENAME_HIDDEN, tp_basic);
}

// The value parameter has to be an string. If supported_values is a null
// pointer, then all arguments are supported. If not, it must be a strings
// array, last element being 0 (NULL)
void settings_add_str(Settings * settings, char * name, int key,
                      char * description, unsigned id, unsigned flags,
                      char * value, char ** supported_values)
{
    Option * option;

    option = settings_add_option(settings, name, key, description, id, flags,
                        tp_string);
    strcpy(option->value.str, value);
    option->supported_args = supported_values;
}


void settings_add_dbl(Settings * settings, char *name, int key,
                      char * description, unsigned id, unsigned flags,
                      double value)
{
    Option * option;

    option = settings_add_option(settings, name, key, description, id, flags,
                        tp_double);
    option->value.dbl = value;
}


void settings_set_str(Settings * settings, char * family, char *name,
                      char * value)
{
    strcpy(settings_get_option(settings, family, name)->value.str, value);
}


void settings_set_dbl(Settings * settings, char * family, char *name,
                      double value)
{
    settings_get_option(settings, family, name)->value.dbl = value;
}


char * settings_get_str(Settings * settings, char * family, char * name)
{
    return settings_get_option(settings, family, name)->value.str;
}


unsigned settings_get_index(Settings * settings, char * family, char * name)
{
    Option * option;
    option =  settings_get_option(settings, family, name);

    if(!string_is_in_list(option->value.str, option->supported_args,
                         &option->index))
        error(1, 0, "Option has not index support\n");

    return option->index;
}


double * settings_get_dbl(Settings * settings, char * family, char * name)
{
    return &(settings_get_option(settings, family, name)->value.dbl);
}


bool settings_is_fixed(Settings * settings, char * family, char * name)
{
    return settings_get_option(settings, family, name)->fixed;
}


bool settings_is_storable(Settings * settings, unsigned id)
{
    return settings_get_option_by_id(settings, id)->id != 0;
}


settings_length_t settings_get_number_of_options_fixed(Settings * settings)
{
    settings_length_t k;
    settings_length_t nopts;

    nopts = 0;

    for(k = 0; k < settings->length; ++k)
    {
        if(settings->options[k].fixed)
            ++nopts;
    }

    return nopts;
}


settings_length_t settings_get_number_of_options_storable(Settings * settings)
{
    settings_length_t length;
    settings_length_t k;

    length = 0;

    for(k = 0; k < settings->length; ++k)
        if(settings_is_storable(settings, settings->options[k].id))
            ++length;

    return length;
}


bool settings_is_sweep_active(Settings * settings)
{
    return settings->sweep.active;
}


// If sweep is not active, length is 1
unsigned settings_get_sweep_length(Settings * settings)
{
    if(settings->sweep.active)
        return settings->sweep.length;
    else
        return 1u;
}


double settings_get_sweep_resolution(Settings * settings)
{
    return settings->sweep.resolution;
}


void settings_apply_sweep(Settings * settings)
{
    if(settings->sweep.active)
    {
        settings->sweep.option->value.dbl += settings->sweep.resolution;
    }
}


void settings_init_sweep(Settings * settings)
{
    if(settings->sweep.active)
        settings->sweep.option->value.dbl = settings->sweep.start_value;
}


void settings_increase_sweep_start_value(Settings * settings,
                                         double start_value)
{
    settings->sweep.start_value += start_value;
}


bool settings_is_sweep_option(Settings * settings, char * family, char * name)
{
    return settings->sweep.option ==
           settings_get_option(settings, family, name);
}


bool option_is_json_printable(Option * option)
{
    return (option->type == tp_string  ||  option->type == tp_double);
}


void settings_add_option_to_json_object(Settings * settings,
                                        settings_length_t index,
                                        JsonObject * jobject)
{
    Option * option;

    option = &(settings->options[index]);

    if(settings->sweep.active  &&  option == settings->sweep.option)
    {
        json_object_set_string_member(jobject, option->name, "sweep");
    }
    else
    {
        if(option->type == tp_double)
        {
            if(option->value.dbl == INFINITY)
                json_object_set_null_member(jobject, option->name);
            else
                json_object_set_double_member(jobject, option->name,
                                                option->value.dbl);
        }
        else if(option->type == tp_string)
        {
            if(strcmp(option->value.str, "") == 0)
                json_object_set_null_member(jobject, option->name);
            else
                json_object_set_string_member(jobject, option->name,
                                            option->value.str);
        }
        else
        {
            error(1, 0, "Option '%s' - '%s' is not json printable\n",
                option->family, option->name);
        }
    }
}


JsonObject * settings_get_json_object(Settings * settings)
{
    JsonObject * jsettings;
    JsonObject * jfamily;
    JsonArray * jsweep;
    Option * option;
    char * family;
    settings_length_t k;
    unsigned sweep_length;
    double * sweep_value;

    jsettings = json_object_new();

    for(k = 0; k < settings->length; ++k)
    {
        option = &settings->options[k];
        family = option->family;
        if(option->json_printable  &&  strcmp(family, GLOBAL_FAMILY) == 0)
            settings_add_option_to_json_object(settings, k, jsettings);
    }

    for(k = 0; k < settings->length; ++k)
    {
        option = &settings->options[k];
        family = option->family;
        if(option->json_printable  &&  strcmp(family, GLOBAL_FAMILY) != 0)
        {
            if(json_object_has_member(jsettings, family))
            {
                jfamily = json_object_get_object_member(jsettings, family);
            }
            else
            {
                jfamily = json_object_new();
                json_object_set_object_member(jsettings, family, jfamily);
            }
            settings_add_option_to_json_object(settings, k, jfamily);
        }
    }

    if(settings->sweep.active)
    {
        jsweep = json_array_new();
        json_object_set_array_member(jsettings, "sweep", jsweep);
        sweep_value = &(settings->sweep.option->value.dbl);
        sweep_length = settings_get_sweep_length(settings);
        settings_init_sweep(settings);
        for(k = 0; k < sweep_length; ++k)
        {
            json_array_add_double_element(jsweep, *sweep_value);
            settings_apply_sweep(settings);
        }
    }

    return jsettings;
}