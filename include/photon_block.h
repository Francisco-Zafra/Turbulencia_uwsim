#ifndef __PHOTON_PROCESSING_H__
#define __PHOTON_PROCESSING_H__


#include "configuration.h"
#include "simulation.h"


// Photon data
typedef struct
{
    // Stored part
    float weight;
    float tof;
    float distance;
    float x, y;
    float theta;
    uint8_t events;
    bool reflected;

    // Temporal part
    float z;
    float ux, uy, uz;
    bool intercepted;
    bool received;
    uint8_t layer;       // añado layer para que se guarde en qué capa está del medio. Cada capa tendrá un índice de refracción
    uint8_t phase_layer;
}
Photon;

// Block of photons
typedef struct
{
    Photon * photons;
    long length;
}
PhotonBlock;


// Create new block of photons
PhotonBlock * photon_block_new(long number_of_photons);

// Destroy block of photons
void photon_block_destroy(PhotonBlock * block);

// Mark all received photons in 'photon_block' with configuration 'sim'
void photon_block_mark_received(PhotonBlock * photon_block, Simulation * sim);

// Process photons in 'photon_block' (emission, tranmission and interception)
void photon_block_process(PhotonBlock * photon_block, Simulation * sim);


#endif