#ifndef __STATISTICS_H__
#define __STATISTICS_H__


#include "photon_block.h"


typedef struct
{
    double sum;
    double square_sum;
}
Stat;

typedef struct
{
    double processing_time;
    long photons_received;
    long events_threshold;
    Stat * stats;
}
Results;



Results * results_new();
Results * results_new_from_photon_block(PhotonBlock * block);
Results * results_new_array(unsigned length);
void results_destroy(Results * results);
void results_destroy_array(Results * results, unsigned length);
void results_reset(Results * results);
void results_add(Results * results, Results * new_results);
void results_add_photons_received(Results * results, PhotonBlock * block);
void results_set_events_threshold(Results * results, long events_threshold);
JsonObject * results_get_json_object(Results * results, unsigned num_results);


#endif