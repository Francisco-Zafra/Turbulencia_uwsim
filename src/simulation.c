#include "simulation.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>

#include "configuration.h"
#include "urand.h"


//#define COS_THETA_LUT_SIZE


typedef struct
{
    double * photons_dbl;
    double * block_size_dbl;
    reception_mode_t mode;
    double * rouletting_threshold;
    double * max_events;
    double * length_factor;
    double * max_time;

    double * med_albedo;
    double * med_attenuation;
    double * med_n_water;
    double * med_layers;
    double * med_boundary_var_z;

    double * sur_depth;
    double * sur_sigma_deg;
    double * sur_n_air;

    SourceType sou_type;
    double * sou_full_angle_deg;
    double * sou_divergence_mrad;
    double * sou_beam_waist_mm;
    double * sou_theta_deg;
    double * sou_phi_deg;
    double * sou_sigma_deg;

    double * rec_x;
    double * rec_y;
    double * rec_z;
    double * rec_x_offset;
    double * rec_y_offset;
    double * rec_aperture_cm;
    double * rec_fov_deg;
    double * rec_theta_deg;
    double * rec_sigma_deg;
    double * rec_phi_deg;
    double * rec_n_air;

    SPFType spf_type;
    double * spf_g;
}
Fundamentals;


float degrees_to_radians(float degrees)
{
    return degrees * M_PI / 180.0f;	// radians
}


void combine_vectors(Vector * sum, Vector * a, float a_scale, Vector * b,
                     float b_scale)
{
    sum->x = a->x * a_scale + b->x * b_scale;
    sum->y = a->y * a_scale + b->y * b_scale;
    sum->z = a->z * a_scale + b->z * b_scale;
}

float dot_product(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return x1 * x2 + y1 * y2 + z1 * z2;
}


void translate_from_system(RotationMatrix * mat, float * x, float * y,
                           float * z)
{
    float nx, ny;

    nx = mat->m11 * *x + mat->m12 * *y + mat->m13 * *z;
    ny = mat->m21 * *x + mat->m22 * *y + mat->m23 * *z;
    *z = mat->m31 * *x + mat->m32 * *y + mat->m33 * *z;
    *x = nx;
    *y = ny;
}


void translate_to_system(RotationMatrix * mat, float * x, float * y, float * z)
{
    float nx, ny;

    nx = mat->m11 * *x + mat->m21 * *y + mat->m31 * *z;
    ny = mat->m12 * *x + mat->m22 * *y + mat->m32 * *z;
    *z = mat->m13 * *x + mat->m23 * *y + mat->m33 * *z;
    *x = nx;
    *y = ny;
}

void get_rotation_matrix(RotationMatrix * mat, float theta, Vector * vector)
{
    float cos_theta;
    float sin_theta;
    float ux_sin_theta;
    float uy_sin_theta;
    float uz_sin_theta;
    float one_m_cos_theta;
    float ux_uy_one_m_cos_theta;
    float ux_uz_one_m_cos_theta;
    float uy_uz_one_m_cos_theta;
    float ux_sq, uy_sq, uz_sq;

    cos_theta = cosf(theta);
    sin_theta = sinf(theta);
    ux_sin_theta = vector->x * sin_theta;
    uy_sin_theta = vector->y * sin_theta;
    uz_sin_theta = vector->z * sin_theta;
    one_m_cos_theta = 1.0f - cos_theta;
    ux_uy_one_m_cos_theta = vector->x * vector->y * one_m_cos_theta;
    ux_uz_one_m_cos_theta = vector->x * vector->z * one_m_cos_theta;
    uy_uz_one_m_cos_theta = vector->y * vector->z * one_m_cos_theta;
    ux_sq = vector->x * vector->x;
    uy_sq = vector->y * vector->y;
    uz_sq = vector->z * vector->z;

    mat->m11 = ux_sq + (uy_sq + uz_sq) * cos_theta;
    mat->m21 = ux_uy_one_m_cos_theta + uz_sin_theta;
    mat->m31 = ux_uz_one_m_cos_theta - uy_sin_theta;

    mat->m12 = ux_uy_one_m_cos_theta - uz_sin_theta;
    mat->m22 = uy_sq + (ux_sq + uz_sq) * cos_theta;
    mat->m32 = uy_uz_one_m_cos_theta + ux_sin_theta;

    mat->m13 = ux_uz_one_m_cos_theta + uy_sin_theta;
    mat->m23 = uy_uz_one_m_cos_theta - ux_sin_theta;
    mat->m33 = uz_sq + (ux_sq + uy_sq) * cos_theta;
}


// Return numbers gaussian distribution with mean 0 and standard deviation
// 'standard_deviation'
float get_gaussian(float standard_deviation)
{
    static bool ready = false;
    static float  gstore;
    float gauss;
    float v1,v2;
    float r;
    float fac;

    if (ready)
    {
        ready = false;
        gauss = gstore;
    }
    else
    {
        do
        {
            v1 = 2 * urand();
            v2 = 2 * urand();
            r = v1 * v1 + v2 * v2;
        }
        while (r > 1.0f);

        fac = sqrt(-2 * log(r) / r);
        gstore = v1 * fac;
        gauss = v2 * fac;
        ready = true;
    }

    return gauss * standard_deviation;
}


float inv_cdf_hg(float q, HGParameters * hgp)
{
    if(hgp->g != 0.0f)
    {
        float frac = hgp->one_minus_g_square /
                     (hgp->g_plus_one - hgp->g_by_two * q);
        return hgp->one_divided_by_g_by_two *
            (hgp->g_square_plus_one - frac * frac);
    }

    return 1.0f - 2.0f * urand();
}


float inv_cdf_spf(float q, Simulation * sim)
{
    if(sim->spf_type == hg)
    {
        return inv_cdf_hg(q, &sim->spf_g);
    }
    else	// TTHG
    {
        if(urand() < sim->spf_alpha)
            return inv_cdf_hg(q, &sim->spf_g);
        else
            return inv_cdf_hg(q, &sim->spf_minus_h);
    }
}


void set_hg_parameters(HGParameters * hgp, float g)
{
    hgp->g = g;
    hgp->g_plus_one = hgp->g + 1.0f;
    hgp->g_by_two = hgp->g * 2.0f;
    hgp->one_divided_by_g_by_two = 1.0f / hgp->g_by_two;
    hgp->g_square = hgp->g * hgp->g;
    hgp->g_square_plus_one = hgp->g_square + 1.0f;
    hgp->one_minus_g_square = 1.0f - hgp->g_square;
}


void get_rotation_vector(Vector * vector, float phi)
{
    vector->x = -cosf(phi);
    vector->y = sinf(phi);
    vector->z = 0.0f;
}


void get_fundamentals_from_settings(Fundamentals * fund, Settings * set)
{
    fund->photons_dbl = settings_get_dbl(set, 0, "photons");
    fund->block_size_dbl = settings_get_dbl(set, 0, "block size");
    fund->mode = (reception_mode_t)settings_get_index(set, 0, "mode");
    fund->rouletting_threshold = settings_get_dbl(set, 0, "rouletting");
    fund->max_events = settings_get_dbl(set, 0, "max events");
    fund->length_factor = settings_get_dbl(set, 0, "length factor");
    fund->max_time = settings_get_dbl(set, 0, "time");

    fund->med_albedo = settings_get_dbl(set, "medium", "albedo");
    fund->med_attenuation = settings_get_dbl(set, "medium", "attenuation");
    fund->med_n_water = settings_get_dbl(set, "medium", "index");
    fund->med_layers = settings_get_dbl(set, "medium", "layers");  // se extrae el número de capas
    fund->med_boundary_var_z = settings_get_dbl(set, "medium", "varZ"); //varZ de boundary

    fund->sur_depth = settings_get_dbl(set, "surface", "depth");
    fund->sur_sigma_deg = settings_get_dbl(set, "surface", "sigma");
    fund->sur_n_air = settings_get_dbl(set, "surface", "index");

    fund->sou_type = (SourceType)settings_get_index(set, "source", "type");
    fund->sou_full_angle_deg = settings_get_dbl(set, "source", "full angle");
    fund->sou_divergence_mrad = settings_get_dbl(set, "source", "divergence");
    fund->sou_beam_waist_mm = settings_get_dbl(set, "source", "beam waist");
    fund->sou_theta_deg = settings_get_dbl(set, "source", "theta");
    fund->sou_phi_deg = settings_get_dbl(set, "source", "phi");
    fund->sou_sigma_deg = settings_get_dbl(set, "source", "sigma");

    fund->rec_x = settings_get_dbl(set, "receptor", "x");
    fund->rec_y = settings_get_dbl(set, "receptor", "y");
    fund->rec_z = settings_get_dbl(set, "receptor", "z");
    fund->rec_x_offset = settings_get_dbl(set, "receptor", "x offset");
    fund->rec_y_offset = settings_get_dbl(set, "receptor", "y offset");
    fund->rec_aperture_cm = settings_get_dbl(set, "receptor", "aperture");
    fund->rec_fov_deg = settings_get_dbl(set, "receptor", "fov");
    fund->rec_theta_deg = settings_get_dbl(set, "receptor", "theta");
    fund->rec_sigma_deg = settings_get_dbl(set, "receptor", "sigma");
    fund->rec_phi_deg = settings_get_dbl(set, "receptor", "phi");
    fund->rec_n_air = settings_get_dbl(set, "receptor", "index");

    fund->spf_type = (SPFType)settings_get_index(set, "medium",
                                                 "phase function");
    fund->spf_g = settings_get_dbl(set, "medium", "g");
}


void get_simulation_from_fundamentals(Simulation * sim, Fundamentals * fund,
                                 float time_resolution)
{
    Vector rvector;
    float half_full_angle;
    float n;
    float scattering;
    float length_factor;
    float h;
    float g;
    float theta, phi;
    float med_attenuation;
    float link_length;


    sim->photons = (long)*fund->photons_dbl;
    sim->block_size = (long)*fund->block_size_dbl;
    sim->mode = fund->mode;
    sim->max_time = (float)*fund->max_time;

    sim->sou_type = fund->sou_type;
    half_full_angle = (float)*fund->sou_full_angle_deg * M_PI / 180.0f / 2.0f;
    n = -logf(2.0f) / logf(cosf(half_full_angle));
    sim->sou_lamb_exponent = (1.0f / (1.0f + n));
    sim->sou_divergence = (float)*fund->sou_divergence_mrad * 0.001f;
    sim->sou_beam_waist = (float)*fund->sou_beam_waist_mm * 0.001f;
    sim->sou_sigma = degrees_to_radians((float)*fund->sou_sigma_deg);
    if(sim->sou_sigma == 0.0f)
        sim->source_vibrating = false;
    else
        sim->source_vibrating = true;

    sim->sou_vx.x = 1.0f;
    sim->sou_vx.y = 0.0f;
    sim->sou_vx.z = 0.0f;

    sim->sou_vy.x = 0.0f;
    sim->sou_vy.y = 1.0f;
    sim->sou_vy.z = 0.0f;

    sim->sou_vz.x = 0.0f;
    sim->sou_vz.y = 0.0f;
    sim->sou_vz.z = 1.0f;

    theta = degrees_to_radians((float)*fund->sou_theta_deg);
    if(theta != 0.0f)
    {
        sim->source_deflected = true;
        phi = degrees_to_radians((float)*fund->sou_phi_deg);
        get_rotation_vector(&rvector, phi);
        get_rotation_matrix(&sim->sou_mat, theta, &rvector);
        translate_from_system(&sim->sou_mat, &sim->sou_vx.x, &sim->sou_vx.y,
                              &sim->sou_vx.z);
        translate_from_system(&sim->sou_mat, &sim->sou_vy.x, &sim->sou_vy.y,
                              &sim->sou_vy.z);
        translate_from_system(&sim->sou_mat, &sim->sou_vz.x, &sim->sou_vz.y,
                              &sim->sou_vz.z);
    }
    else
        sim->source_deflected = false;

    sim->rec_x = (float)*fund->rec_x;
    sim->rec_y = (float)*fund->rec_y;
    sim->rec_z = (float)*fund->rec_z;
    sim->rec_x_offset = (float)*fund->rec_x_offset;
    sim->rec_y_offset = (float)*fund->rec_y_offset;
    sim->rec_radius = (float)(*fund->rec_aperture_cm / (2.0 * 100.0));
    sim->rec_radius_square = sim->rec_radius * sim->rec_radius;
    sim->rec_half_fov = (float)(*fund->rec_fov_deg / 180.0 / 2.0 * M_PI);
    sim->rec_cos_half_fov = cosf(sim->rec_half_fov);


    sim->rec_vx.x = 1.0f;
    sim->rec_vx.y = 0.0f;
    sim->rec_vx.z = 0.0f;
    sim->rec_vy.x = 0.0f;
    sim->rec_vy.y = 1.0f;
    sim->rec_vy.z = 0.0f;
    sim->rec_nz.x = 0.0f;
    sim->rec_nz.y = 0.0f;
    sim->rec_nz.z = 1.0f;

    theta = degrees_to_radians((float)*fund->rec_theta_deg);
    if(theta != 0.0f)
    {
        sim->receptor_deflected = true;
        phi = degrees_to_radians((float)*fund->rec_phi_deg);
        get_rotation_vector(&rvector, phi);
        get_rotation_matrix(&sim->rec_mat, theta, &rvector);

        translate_from_system(&sim->rec_mat, &sim->rec_vx.x, &sim->rec_vx.y,
                              &sim->rec_vx.z);
        translate_from_system(&sim->rec_mat, &sim->rec_vy.x, &sim->rec_vy.y,
                              &sim->rec_vy.z);
        translate_from_system(&sim->rec_mat, &sim->rec_nz.x,
                              &sim->rec_nz.y, &sim->rec_nz.z);
    }
    else
    {
        sim->receptor_deflected = false;
    }

    sim->rec_sigma = degrees_to_radians((float)*fund->rec_sigma_deg);
    if(sim->rec_sigma == 0.0)
        sim->receptor_vibrating = false;
    else
        sim->receptor_vibrating = true;

    sim->med_layers = (uint8_t)*fund->med_layers;  //número de capas
    sim->med_boundary_var_z = (float)*fund->med_boundary_var_z;
    sim->med_albedo = (float)*fund->med_albedo;
    med_attenuation = (float)*fund->med_attenuation;
    sim->med_minus_inv_c = -1.0f / med_attenuation;
    if(sim->med_albedo == 0.0f)
        sim->rouletting_threshold = 0;
    else if(*fund->rouletting_threshold == 0.0)
        sim->rouletting_threshold = 65535;  // MAX unit16_t, just in case
    else
        sim->rouletting_threshold = ceilf(logf(*fund->rouletting_threshold)/
                                              logf(sim->med_albedo));

    length_factor = (float)*fund->length_factor;
    scattering = med_attenuation * sim->med_albedo;
    link_length = sqrtf(powf(sim->rec_z, 2.0) + powf(sim->rec_y, 2.0) +
                  powf(sim->rec_x, 2.0));
    sim->events_threshold = ceilf(scattering * link_length * length_factor);
    if(sim->events_threshold < 10)
        sim->events_threshold = 10;
    if((double)sim->events_threshold > *fund->max_events)
        sim->events_threshold = (unsigned)*fund->max_events;
    if(sim->events_threshold > 255)
        sim->events_threshold = 255;

    sim->med_n_water = (float)*fund->med_n_water;
    sim->rec_n_air = (float)*fund->rec_n_air;
    sim->sur_n_air = (float)*fund->sur_n_air;
    sim->rec_n_quotient = sim->med_n_water / sim->rec_n_air;
    sim->rec_cos_critical_angle = sqrtf(1.0f - powf(sim->rec_n_air /
                                                    sim->med_n_water, 2.0f));

    sim->sur_n_quotient = sim->med_n_water / sim->sur_n_air;
    sim->sur_cos_critical_angle = sqrtf(1.0f - powf(sim->sur_n_air /
                                                    sim->med_n_water, 2.0f));
    sim->med_speed_m_ns = SPEED_VACUUM / sim->med_n_water / 1e9f;

    float source_surface_dist = (float)*fund->sur_depth;
    if(source_surface_dist != INFINITY)
    {
        sim->there_is_surface = true;
        sim->receptor_depth = source_surface_dist - sim->rec_y;
    }
    else
        sim->there_is_surface = false;
    sim->sur_sigma = degrees_to_radians((float)*fund->sur_sigma_deg);
    if(sim->sur_sigma == 0.0)
        sim->surface_vibrating = false;
    else
        sim->surface_vibrating = true;

    sim->spf_type = fund->spf_type;
    g = (float)*fund->spf_g;
    set_hg_parameters(&sim->spf_g, g);
    if(sim->spf_type == tthg)
    {
        h = -0.3061446 + 1.000568 * g - 0.01826332 * g*g + 0.03643748 * g*g*g;
        set_hg_parameters(&sim->spf_minus_h, -h);
        sim->spf_alpha = (h * (1 + h)) / ((g + h) * (1 + h - g));
    }

    sim->rotation_required = sim->receptor_deflected  ||
                             sim->receptor_vibrating;
}


Simulation * simulation_new_conditionally(Settings * set, bool sweep_array)
{
    Simulation * sim;
    Fundamentals fund;
    unsigned sweep_length;
    float time_resolution;
    unsigned k;

    time_resolution = INFINITY;

    if(sweep_array)
    {
        sweep_length = settings_get_sweep_length(set);
        if(settings_is_sweep_option(set, 0, "time"))
            time_resolution = (float)settings_get_sweep_resolution(set);
        settings_init_sweep(set);
    }
    else
        sweep_length = 1;

    sim = malloc(sweep_length * sizeof(Simulation));
    get_fundamentals_from_settings(&fund, set);

    settings_init_sweep(set);
    for(k = 0; k < sweep_length; ++k)
    {
        get_simulation_from_fundamentals(&sim[k], &fund, time_resolution);
        if(sweep_array)
            settings_apply_sweep(set);
    }
    // Reset initial
    settings_init_sweep(set);

    return sim;
}


// set must contain input data
Simulation * simulation_new_from_settings_start(Settings * set)
{
    return simulation_new_conditionally(set, false);
}


// set must contain input data
Simulation * simulation_new_array_from_settings(Settings * set)
{
    return simulation_new_conditionally(set, true);
}


void simulation_destroy(Simulation * sim)
{
    free(sim);
}


void simulation_destroy_array(Simulation * sim)
{
    free(sim);
}


bool simulation_is_compatible(Simulation * sim, Simulation * father_sim)
{
    float radius_margin;
    float x_dist, y_dist;
    float distance;

    if(father_sim->events_threshold < sim->events_threshold)
        return false;

    if((father_sim->mode != sim->mode)  &&  father_sim->mode != all)
        return false;

    if(sim->rec_half_fov > father_sim->rec_half_fov)
        return false;

    radius_margin = father_sim->rec_radius - sim->rec_radius;
    if(radius_margin < 0.0)
        return false;

    x_dist = sim->rec_x_offset - father_sim->rec_x_offset;
    y_dist = sim->rec_y_offset - father_sim->rec_y_offset;
    distance = sqrt(x_dist * x_dist + y_dist * y_dist);
    if(distance > radius_margin)
        return false;

    return true;
}