#include "photon_block.h"

#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#include "urand.h"

#include <stdio.h> //añadido para el printf

#define UZ_THRESHOLD 0.99999f



/* --- Auxiliar functions --------------------------------------------------- */

// Translate vector of photons to system defined by 'mat'
void translate_photon_to_system(Photon * photon, RotationMatrix * mat)
{
    translate_to_system(mat, &photon->x, &photon->y, &photon->z);
    translate_to_system(mat, &photon->ux, &photon->uy, &photon->uz);
}


// Translate vector of photons from system defined by 'mat'
void translate_photon_from_system(Photon * photon, RotationMatrix * mat)
{
    translate_from_system(mat, &photon->x, &photon->y, &photon->z);
    translate_from_system(mat, &photon->ux, &photon->uy, &photon->uz);
}


// Return true if photon is intercepted by receptor plane with normal vector
// (nx,ny,nz)
bool is_intercepted(Photon * photon, float nx, float ny,
                    float nz)
{
    return dot_product(photon->x, photon->y, photon->z, nx, ny ,nz) > 0;
}

// Reflect photon to the other side of the surface with direction reflected
void reflect_photon(Photon * photon, Simulation * sim)
{
    RotationMatrix surface_mat; // Rotation matrix for vibration modeling
    Vector rvector;             // Rotation vector for 'surface_mat'
    float distance_to_surface;  // Distance traveled since surface reflection
    float rp, rs;
    float R;
    float cos_exit;
    float theta, phi;

    photon->reflected = true;
    // Calculate distance traveled since surface reflection
    distance_to_surface = (photon->y - sim->receptor_depth) /
                            photon->uy;
    // If receptor depth is 4m, then y position of the photon to be in the
    // surface must be 4m
    photon->y = sim->receptor_depth;

    // If surface is not horizontal, translate coordinates
    if(sim->surface_vibrating)
    {
        // If surface is not horizontal, then x and y change after reflection, so they are brought back to the interception position
        photon->x -= distance_to_surface * photon->ux;
        photon->z -= distance_to_surface * photon->uz;

        // Generate random rotation vector
        phi = urand() * 2.0f * M_PI;
        theta = get_gaussian(sim->sur_sigma);
        rvector.x = cosf(phi);
        rvector.z = sinf(phi);
        rvector.y = 0.0f;
        // Get rotation matrix from vector
        get_rotation_matrix(&surface_mat, theta, &rvector);
        // Translate photon direction to rotated system
        translate_to_system(&surface_mat, &photon->ux,
                            &photon->uy, &photon->uz);
    }

    // If incident angle is lower than critical angle, apply Fresnel equations
    if(photon->uy > sim->sur_cos_critical_angle)
    {
        cos_exit = sqrtf( 1.0f - sim->sur_n_quotient * sim->sur_n_quotient *
                         (1.0f - photon->uy * photon->uy));
        rp = (sim->sur_n_air * photon->uy - sim->med_n_water * cos_exit) /
             (sim->sur_n_air * photon->uy + sim->med_n_water * cos_exit);
        rs = (sim->med_n_water * photon->uy - sim->sur_n_air * cos_exit) /
             (sim->med_n_water * photon->uy + sim->sur_n_air * cos_exit);
        R = (rp * rp + rs  *rs) / 2.0f;
        photon->weight *= R;
    }

    // Incident angle is equal to reflexed angle
    photon->uy = -photon->uy;

    // If surface is not horizontal, inverse translation
    if(sim->surface_vibrating)
    {
        translate_from_system(&surface_mat, &photon->ux,
                              &photon->uy, &photon->uz);

        // Photon must travel the remain distance to fulfill Beer-Lambert law
        photon->x += distance_to_surface * photon->ux;
        photon->z += distance_to_surface * photon->uz;
    }

    // Photon must travel the remain distance to fulfill Beer-Lambert law
    photon->y += distance_to_surface * photon->uy;
}


// Update photon direction with polar angle cosine 'cos_theta' and azimuthal
// angle 'phi'
void update_photon_direction(Photon * photon, float cos_theta, float phi)
{
    float ux_new, uy_new;
    float sin_theta;
    float cos_phi;
    float sin_phi;
    float sqrt_uz;

    sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
    cos_phi = cosf(phi);
    sin_phi = sinf(phi);

    if(photon->uz > UZ_THRESHOLD)
    {
        // uz very close to 1
        photon->ux = sin_theta * cos_phi;
        photon->uy = sin_theta * sin_phi;
        photon->uz = cos_theta;
    }
    else if(photon->uz < -UZ_THRESHOLD)
    {
        // uz very close to -1
        photon->ux = sin_theta * cos_phi;
        photon->uy = sin_theta * sin_phi;
        photon->uz = -cos_theta;
    }
    else
    {
        sqrt_uz = sqrtf(1.0f - photon->uz * photon->uz);

        ux_new = (photon->ux * cos_theta) + (sin_theta / sqrt_uz) *
                 (photon->ux * photon->uz * cos_phi - photon->uy * sin_phi);
        uy_new = (photon->uy * cos_theta) + (sin_theta/sqrt_uz) *
                 (photon->uy * photon->uz * cos_phi + photon->ux * sin_phi);
        photon->uz = ((-sin_theta) * cos_phi) * sqrt_uz + photon->uz *
                      cos_theta;
        photon->ux = ux_new;
        photon->uy = uy_new;
    }
}


// Generate initial 'photon' state with configuration 'sim'
void photon_emit(Photon * photon, Simulation * sim)
{
    RotationMatrix vib_mat; // Rotation matrix for source vibration modeling
    Vector rvector;         // Rotation vector for 'vib_mat'
    float norm_radius;      // 'radius' normalized by 'beam_waist'
    float radius;           // Distance to the center of the source
    float theta;
    float cos_theta;
    float sin_theta;
    float phi;
    float cos_phi;
    float sin_phi;

    if(sim->sou_type == gaussian)
    {
        // Gaussian source
        norm_radius = sqrtf(-logf(urand()));
        radius = sim->sou_beam_waist * norm_radius;
        theta = sim->sou_divergence * norm_radius;
        sin_theta = sinf(theta);

        phi = (2.0f * M_PI) * urand();
        cos_phi = cosf(phi);
        sin_phi = sinf(phi);

        photon->x = radius * cos_phi;
        photon->y = radius * sin_phi;
        photon->z = 0.0f;

        photon->uz = cosf(theta);
        photon->ux = sin_theta * cos_phi;
        photon->uy = sin_theta * sin_phi;

        if(sim->source_deflected)
            // Apply source deflection
            translate_photon_from_system(photon, &sim->sou_mat);
        if(sim->source_vibrating)
        {
            // Genrate random rotation vector
            phi = (2.0f * M_PI) * urand();
            theta = get_gaussian(sim->sou_sigma);
            combine_vectors(&rvector, &sim->sou_vx, cosf(phi), &sim->sou_vy,
                            sinf(phi));
            // Get rotation matrix from vector
            get_rotation_matrix(&vib_mat, theta, &rvector);
            // Apply vibration
            translate_photon_from_system(photon, &vib_mat);
        }

        // Move photon coordinates to receptor origin
        photon->x -= sim->rec_x;
        photon->y -= sim->rec_y;
        photon->z -= sim->rec_z;
    }
    else // Point sources
    {
        // Set initial position with respect to the receptor
        photon->x = -sim->rec_x;
        photon->y = -sim->rec_y;
        photon->z = -sim->rec_z;

        if(sim->sou_type == isotropic)
        {
            // Isotropic source (there is no need to model vibration)
            cos_theta = 1.0f - 2.0f * urand();
            sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
            phi = 2.0f * M_PI * urand();
            cos_phi = cosf(phi);
            sin_phi = sinf(phi);
            photon->ux = cos_phi * sin_theta;
            photon->uy = sin_phi * sin_theta;
            photon->uz = cos_theta;
        }
        else // Point sources that are inclination dependent
        {
            if(sim->sou_type == ideal)
            {
                // Ideal source
                photon->ux = sim->sou_vz.x;
                photon->uy = sim->sou_vz.y;
                photon->uz = sim->sou_vz.z;
            }
            else if(sim->sou_type == lambertian)
            {
                // Lambertian source
                cos_theta = powf(urand(), sim->sou_lamb_exponent);
                phi = (2.0f * M_PI) * urand();

                if(sim->source_deflected)
                {
                    photon->ux = sim->sou_vz.x;
                    photon->uy = sim->sou_vz.y;
                    photon->uz = sim->sou_vz.z;
                    update_photon_direction(photon, cos_theta, phi);
                }
                else
                {
                    cos_phi = cosf(phi);
                    sin_phi = sinf(phi);
                    sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
                    photon->uz = cos_theta;
                    photon->ux = sin_theta * cos_phi;
                    photon->uy = sin_theta * sin_phi;
                }
            }

            if(sim->source_vibrating)
            {
                // Translate only the direction (position is translated
                // directly)
                phi = (2.0f * M_PI) * urand();
                theta = get_gaussian(sim->sou_sigma);
                update_photon_direction(photon, cosf(theta), phi);
            }
        }
    }

    photon->weight = 1.0f;
    photon->distance = 0.0f;
    photon->events = 0u;
    photon->reflected = false;
    photon->intercepted = false;
    photon->received = false;
    photon->layer = 1u;    // se inicializa estando en la capa 1 por si se considera modelar el medio con múltiples capas y diferentes índices de refracción
}


// Move 'photon' from source plane to receptor plane, unless events threshold
// or rouletting threshold is exceded, with configuration 'sim'
void photon_move(Photon * photon, Simulation * sim)
{
    static unsigned rouletting = 0; // Index for rouletting sequence control
    Photon new_photon;              // Photon clone to make tests
    RotationMatrix mat;             // Rotation matrix to model vibration
    Vector rv;                      // Rotation vector for 'mat'
    float nx, ny, nz;               // Normal vector to receptor plane
    float r;
    float distance_to_surface;
    float distance_to_rec;
    float distance_to_boundary;
    float theta, phi;
    bool intercepted;
    bool out;
    bool boundary;
    //float n_water_variable[10] = {0, 1.33f, 1.23f, 1.3f, 1.29f, 1.32f, 1.15f, 1.33f, 1.28f, 1.3f};
    float n_quotient;
    float boundary_cos_critical_angle;
    float boundary_normal_x = 0, boundary_normal_y = 0, boundary_normal_z = 1; //Probamos con n = (0,0,1)
    float ux_exit, uy_exit, uz_exit;
    float beta;
    float rp,rs,R,T;

    // This initialization is important in the cases when there is surface but
    // not rotation is required
    nx = sim->rec_nz.x;
    ny = sim->rec_nz.y;
    nz = sim->rec_nz.z;

    if(sim->rotation_required)
    {
        if(sim->there_is_surface)
        {
            nx = sim->rec_nz.x;
            ny = sim->rec_nz.y;
            nz = sim->rec_nz.z;
            if(sim->receptor_vibrating)
            {
                phi = urand() * 2.0f * M_PI;
                theta = get_gaussian(sim->rec_sigma);
                combine_vectors(&rv, &sim->rec_vx, cosf(phi), &sim->rec_vy,
                                sinf(phi));
                get_rotation_matrix(&mat, theta, &rv);
                translate_from_system(&mat, &nx, &ny, &nz);
            }
        }
        else
        {
            // If there is no surface, photon vectors are translated to receptor
            // frame
            if(sim->receptor_vibrating)
            {
                // Generate rotation vector
                phi = urand() * 2.0f * M_PI;
                theta = get_gaussian(sim->rec_sigma);
                combine_vectors(&rv, &sim->rec_vx, cosf(phi), &sim->rec_vy,
                                sinf(phi));
                // Get rotation matrix from rotation vector
                get_rotation_matrix(&mat, theta, &rv);
                // Apply vibration
                translate_photon_to_system(photon, &mat);
            }
            if(sim->receptor_deflected)
                // Apply deflection
                translate_photon_to_system(photon, &sim->rec_mat);
        }

    }

    // If photon is intercepted before it has moved, it is discarded
    if(sim->there_is_surface)
    {
        if(is_intercepted(photon, nx, ny, nz))
            return;
    }
    else
    {
        // In this case, coordinates are related to the receptor frame
        if(photon->z >= 0.0f)
            return;
    }



    while(true)
    {
        // Update position
        r = sim->med_minus_inv_c * logf(urand());
        photon->distance += r;
        photon->x += r * photon->ux;
        photon->y += r * photon->uy;
        photon->z += r * photon->uz;

        // Tengo >2 layers? He pasado la proxima layer?
        boundary = (photon->layer < sim->med_layers) && 
                    dot_product(photon->x, photon->y, photon->z, boundary_normal_x, boundary_normal_y ,boundary_normal_z) > 
                    sim->med_boundary_pos[photon->layer - 1];
        
        while (boundary)    //chequeamos cambio de capa
        {
            //Calculate normal vector with polar and azimuthal angles
            theta = urand() * (sim->med_boundary_max_theta*M_PI/2.0f);
            //phi = get_gaussian(M_PI/2.0f);
            phi = urand() * (2.0f * M_PI);
            boundary_normal_x = cosf(phi)*sinf(theta);
            boundary_normal_y = sinf(phi)*sinf(theta);
            boundary_normal_z = cosf(theta);

            //Calcualte n quotient
            n_quotient=sim->med_n_water_variables[photon->layer]/sim->med_n_water_variables[photon->layer-1];
            boundary_cos_critical_angle= sqrtf(1.0f - powf(n_quotient, 2.0f));

            if((sim->med_n_water_variables[photon->layer-1]>sim->med_n_water_variables[photon->layer])&&(photon->uz < boundary_cos_critical_angle)){
                // si reflexión total en la transición de medios entonces se descarta el fotón
                break;
            }

            // Bring back photon to boundary position
            distance_to_boundary = (photon->z-sim->med_boundary_pos[photon->layer - 1]) / photon->uz;
            //printf("Distancia: %f\n", distance_to_boundary);
            photon->x -= distance_to_boundary * photon->ux;
            photon->y -= distance_to_boundary * photon->uy;
            photon->z -= distance_to_boundary * photon->uz;
            

            //Calculate new photon direction
            float b = 2*n_quotient*dot_product(photon->ux,photon->uy,photon->uz,boundary_normal_x,boundary_normal_y,boundary_normal_z);
            float c = 1 - powf(n_quotient, 2.0f);

            beta = (-b + sqrtf(powf(b,2) - 4*c)) / 2;        
            ux_exit = n_quotient * photon->ux + beta * boundary_normal_x;
            uy_exit = n_quotient * photon->uy + beta * boundary_normal_y;
            uz_exit = n_quotient * photon->uz + beta * boundary_normal_z;
            if ((dot_product(photon->ux,photon->uy,photon->uz,boundary_normal_x,boundary_normal_y,boundary_normal_z)>0 &&
                dot_product(ux_exit,uy_exit,uz_exit,boundary_normal_x,boundary_normal_y,boundary_normal_z)<0) ||
                (dot_product(photon->ux,photon->uy,photon->uz,boundary_normal_x,boundary_normal_y,boundary_normal_z)<0 &&
                dot_product(ux_exit,uy_exit,uz_exit,boundary_normal_x,boundary_normal_y,boundary_normal_z)>0))
            {
                beta = (-b - sqrtf(powf(b,2) - 4*c)) / 2;
                ux_exit = n_quotient * photon->ux + beta * boundary_normal_x;
                uy_exit = n_quotient * photon->uy + beta * boundary_normal_y;
                uz_exit = n_quotient * photon->uz + beta * boundary_normal_z;         
            }

            //printf("coordenada del fotón (%f,%f,%f)\n",photon->x, photon->y, photon->z);
            //printf("direccion del fotón (%f,%f,%f)\n",ux_exit, uy_exit, uz_exit);
            //printf("capa %d, pos %f\n",photon->layer, -(sim->rec_z-photon->layer*(sim->rec_z/sim->med_layers)));

            // Apply Fresnel equations
            rp = (sim->med_n_water_variables[photon->layer] * photon->uz -
                        sim->med_n_water_variables[photon->layer-1] * uz_exit) /
                        (sim->med_n_water_variables[photon->layer] * photon->uz +
                        sim->med_n_water_variables[photon->layer-1] * uz_exit);
            rs = (sim->med_n_water_variables[photon->layer-1] * photon->uz -
                        sim->med_n_water_variables[photon->layer] * uz_exit) /
                        (sim->med_n_water_variables[photon->layer-1] * photon->uz +
                        sim->med_n_water_variables[photon->layer] * uz_exit);
            R = (rp * rp + rs  *rs) / 2.0f;
            T = 1.0f - R;
            // Calculate final weight
            photon->weight *= T;
            photon->layer += 1;

            //Update photon trayectory 
            photon->ux = ux_exit;
            photon->uy = uy_exit;
            photon->uz = uz_exit;

            // Photon must travel the remain distance to fulfill Beer-Lambert law
            photon->x += distance_to_boundary * photon->ux;
            photon->y += distance_to_boundary * photon->uy;
            photon->z += distance_to_boundary * photon->uz;

            //printf("coordenada z del fotón %5.2f\n",photon->z);
            //printf("capa %d\n",photon->layer);

            //Por si me he saltado varias layers
            boundary = (photon->layer < sim->med_layers) && 
                        dot_product(photon->x, photon->y, photon->z, boundary_normal_x, boundary_normal_y ,boundary_normal_z) > 
                        sim->med_boundary_pos[photon->layer - 1];
        }
        
        // Check interceptions with receptor (and surface, if there is)
        if(sim->there_is_surface)
        {
            intercepted = is_intercepted(photon, nx, ny, nz);
            out = photon->y > sim->receptor_depth;

            if(intercepted && out)
            {
                // A photon clone is brought back to the reflection position
                new_photon = *photon;
                distance_to_surface = (new_photon.y - sim->receptor_depth) /
                                      new_photon.uy;
                new_photon.x -= distance_to_surface * new_photon.ux;
                new_photon.y = sim->receptor_depth;
                new_photon.z -= distance_to_surface * new_photon.uz;

                // If in this position, photon clone is intercepted, then
                // interception takes place before reflection
                intercepted = is_intercepted(&new_photon, nx, ny, nz);
                out = !intercepted;
            }

            if(out)
            {
                // Reflection
                if(photon->reflected)
                    // A double reflexion is avoided
                    break;
                reflect_photon(photon, sim);
                if(photon->y >= 0.0f)
                    // Photon out
                    break;

                // Before reflection, photon can be intercepted by receptor
                intercepted = is_intercepted(photon, nx, ny, nz);
            }

            if(intercepted)
            {
                // In this case, coordinates are related to the universal frame
                photon->intercepted = true;
                if(sim->receptor_vibrating)
                        translate_photon_to_system(photon, &mat);
                if(sim->receptor_deflected)
                    translate_photon_to_system(photon, &sim->rec_mat);
            }
        }
        else
        {
            // In this case, coordinates are related to the receptor frame
            photon->intercepted = photon->z >= 0.0f;
        }

        if(photon->intercepted)
        {
            // Bring back photon to interception position
            distance_to_rec = photon->z / photon->uz;
            photon->x -= distance_to_rec * photon->ux;
            photon->y -= distance_to_rec * photon->uy;
            photon->distance -= distance_to_rec;
            break;
        }

        // Check number of events
        if(photon->events == sim->events_threshold)
            break;
        // Update number of events
        photon->events++;

        // Check rouletting threshold
        if(photon->events >= sim->rouletting_threshold)
        {
            // Every 'ROULETTING_MAX_VALUE' photons, one survive
            if(rouletting < ROULETTING_MAX_VALUE)
            {
                rouletting++;
                break;
            }
            else
            {
                photon->weight *= ROULETTING_WEIGHT_INCREASE;
                rouletting = 0;
            }
        }

        // Generate scattering angles
        float cos_theta = inv_cdf_spf(urand(), sim);
        float phi = 2 * M_PI * urand();
        update_photon_direction(photon, cos_theta, phi);
    }
}


// Check reception and calculate final 'photon' state with configuration 'sim'
void photon_receive(Photon * photon, Simulation * sim)
{
    float uz_exit;         // Exit polar angle cosine
    float distance_square; // Distance square to the center of the receptor
    float x_dist, y_dist;
    float rp,rs,R,T;

    if(photon->intercepted)
    {
        // Check critical angle and mode
        if(photon->uz > sim->rec_cos_critical_angle  &&
            (sim->mode == all  ||
            (sim->mode == los && !photon->reflected)  ||
            (sim->mode == nlos && photon->reflected)))
        {
            // Calculate distance to the receptor center and exit angle
            x_dist = photon->x - sim->rec_x_offset;
            y_dist = photon->y - sim->rec_y_offset;
            distance_square = x_dist * x_dist + y_dist * y_dist;
            uz_exit = sqrtf(1.0f - sim->rec_n_quotient * sim->rec_n_quotient
                            * (1.0f - photon->uz * photon->uz));

            // Check entry
            if(distance_square <= sim->rec_radius_square  &&
                uz_exit >= sim->rec_cos_half_fov)
            {
                photon->received = true;

                // Apply Fresnel equations
                rp = (sim->rec_n_air * photon->uz -
                        sim->med_n_water * uz_exit) /
                        (sim->rec_n_air * photon->uz +
                        sim->med_n_water * uz_exit);
                rs = (sim->med_n_water * photon->uz -
                        sim->rec_n_air * uz_exit) /
                        (sim->med_n_water * photon->uz +
                        sim->rec_n_air * uz_exit);
                R = (rp * rp + rs  *rs) / 2.0f;
                T = 1.0f - R;
                // Calculate final weight
                photon->weight *= T;
                photon->weight *= powf(sim->med_albedo,
                                            photon->events);
                // Calculate time of flight
                photon->tof = photon->distance / sim->med_speed_m_ns;
                // Calculate polar angel of arrival
                photon->uz = uz_exit;
                if(photon->uz >= 1.0f)
                    photon->theta = 0.0f;
                else
                    photon->theta = acosf(photon->uz);
            }
        }
    }

}


/* --- Public functions ----------------------------------------------------- */

PhotonBlock * photon_block_new(long number_of_photons)
{
    PhotonBlock * block;

    block = malloc(sizeof(PhotonBlock));
    block->length = number_of_photons;
    block->photons = malloc(number_of_photons * sizeof(Photon));

    return block;
}


void photon_block_destroy(PhotonBlock * block)
{
    free(block->photons);
    free(block);
}


void photon_block_process(PhotonBlock * photon_block, Simulation * sim)
{
    unsigned k;

    for(k = 0; k < photon_block->length; ++k)
    {

        photon_emit(&photon_block->photons[k], sim);
        photon_move(&photon_block->photons[k], sim);
        photon_receive(&photon_block->photons[k], sim);
    }
}


void photon_block_mark_received(PhotonBlock * photon_block, Simulation * sim)
{
    Photon * photons;      // Array of photons to process
    float x_dist, y_dist;
    float distance_square;
    long k;

    photons = photon_block->photons;

    for(k = 0; k < photon_block->length; ++k)
    {
        photons[k].received = false;

        // Check conditions
        if((sim->mode == all ||  (sim->mode == nlos && photons[k].reflected)  ||
           (sim->mode == los && !photons[k].reflected))  &&
           photons[k].events <= sim->events_threshold  &&
           photons[k].tof < sim->max_time  &&
           photons[k].theta < sim->rec_half_fov)
        {
            // Calculate distance to the receptor center
            x_dist = photons[k].x - sim->rec_x_offset;
            y_dist = photons[k].y - sim->rec_y_offset;
            distance_square = x_dist * x_dist + y_dist * y_dist;
            // Check entry
            if(distance_square < sim->rec_radius_square)
                photons[k].received = true;
        }
    }
}