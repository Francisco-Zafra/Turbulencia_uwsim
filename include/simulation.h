#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <stdbool.h>

#include "settings.h"


typedef enum {ideal, isotropic, gaussian, lambertian} SourceType;
typedef enum {hg, tthg} SPFType;    // Scattering Phase Function Type
typedef enum {all, los, nlos} reception_mode_t;


typedef struct
{
    float g;
    float g_plus_one;
    float g_by_two;
    float one_divided_by_g_by_two;
    float g_square;
    float g_square_plus_one;
    float one_minus_g_square;
}
HGParameters;


typedef struct
{
    float m11, m12, m13;
    float m21, m22, m23;
    float m31, m32, m33;
}
RotationMatrix;

typedef struct
{
    float x;
    float y;
    float z;
}
Vector;

typedef struct
{
    // Global values
    long photons;
    long block_size;
    float max_time;
    unsigned events_threshold;
    reception_mode_t mode;

    // Source values
    SourceType sou_type;
    float sou_lamb_exponent;
    float sou_divergence;
    float sou_beam_waist;
    float sou_sigma;
    //float sou_theta;
    //float sou_phi;
    RotationMatrix sou_mat;
    Vector sou_vx;
    Vector sou_vy;
    Vector sou_vz;
    bool source_deflected;
    bool source_vibrating;

    // Reception values
    float rec_x;
    float rec_y;
    float rec_z;
    float rec_x_offset;
    float rec_y_offset;
    float rec_sigma;
    float rec_radius;
    float rec_radius_square;
    float rec_half_fov;
    float rec_cos_half_fov;
    float rec_n_air;
    float rec_n_quotient;
    float rec_cos_critical_angle;
    RotationMatrix rec_mat;
    Vector rec_vx;
    Vector rec_vy;
    Vector rec_nz;
    bool receptor_deflected;
    bool receptor_vibrating;

    // Medium values
    float med_albedo;
    float med_minus_inv_c;
    unsigned rouletting_threshold;
    float med_n_water;
    float med_speed_m_ns;
    uint8_t med_layers;
    float med_boundary_var_z;
    float* med_n_water_variables;
    float* med_boundary_pos;
    float med_var_n_water;
    float med_boundary_max_theta;

    // Surface
    float sur_sigma;
    float sur_n_air;
    float sur_n_quotient;
    float sur_cos_critical_angle;
    float receptor_depth;
    bool there_is_surface;
    bool surface_vibrating;

    // Scattering Phase Function
    SPFType spf_type;
    float spf_alpha;
    HGParameters spf_g;
    HGParameters spf_minus_h;

    bool rotation_required;

    //Floor
    float floor_depth;
    bool there_is_floor;
    float floor_n;

    //Phase screens
    double*** phase_derivadasX;
    double*** phase_derivadasY;
    int phase_max_x;
    int phase_max_y;
    float phase_resolution; //Multiplico 1e3 por metros para pasarlo a mm
    uint8_t phase_layers;
    float* phase_layer_pos;
    char* phase_json;
}
Simulation;


float dot_product(float x1, float y1, float z1, float x2, float y2, float z2);
void combine_vectors(Vector * sum, Vector * a, float a_scale, Vector * b,
                     float b_scale);
void get_rotation_matrix(RotationMatrix * mat, float theta, Vector * vector);
void translate_from_system(RotationMatrix * mat, float * x, float * y,
                           float * z);
void translate_to_system(RotationMatrix * mat, float * x, float * y, float * z);

Simulation * simulation_new_from_settings_start(Settings * settings);
Simulation * simulation_new_array_from_settings(Settings * settings);
void simulation_destroy(Simulation * sim);
void simulation_destroy_array(Simulation * sim);
bool simulation_is_compatible(Simulation * sim, Simulation * father_sim);
float get_gaussian(float standard_deviation);
float inv_cdf_spf(float q, Simulation * sim); // TODO: mover a otro lado
void init_water_n_and_boundarys(Simulation* sim);
void initParametrosPantallaFaseFromJson(char* json, Simulation* sim);

#endif