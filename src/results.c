#include "results.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include "json-glib/json-glib.h"

#define NUMBER_OF_STATS		6

#define WEIGHT				0
#define THETA               1
#define TOF                 2
#define DISTANCE            3
#define RADIUS              4
#define EVENTS				5

static const char * stats_names[] =
{
    "weight",
    "theta",
    "tof",
    "distance",
    "radius",
    "events"
};


void init_results(Results * results)
{
    results->stats = malloc(NUMBER_OF_STATS * sizeof(Stat));
    results_reset(results);
}


Results * results_new()
{
    Results * results;

    results = malloc(sizeof(Results));
    init_results(results);

    return results;
}


Results * results_new_from_photon_block(PhotonBlock * block)
{
    Results * results;

    results = malloc(sizeof(Results));
    init_results(results);

    results_add_photons_received(results, block);

    return results;
}


Results * results_new_array(unsigned length)
{
    Results * results;
    unsigned k;

    results = malloc(length * sizeof(Results));
    for(k = 0; k < length; ++k)
        init_results(&results[k]);

    return results;
}


void results_reset(Results * results)
{
    results->photons_received = 0;
    results->processing_time = 0.0;
    memset(results->stats, 0, NUMBER_OF_STATS * sizeof(Stat));
}


void results_destroy(Results * results)
{
    free(results->stats);
    free(results);
}


void results_destroy_array(Results * results, unsigned length)
{
    unsigned k;

    for(k = 0; k < length; ++k)
        free(results[k].stats);
    free(results);
}


void results_add(Results * results, Results * new_results)
{
    unsigned k;

    results->photons_received += new_results->photons_received;
    for(k = 0; k < NUMBER_OF_STATS; ++k)
    {
        results->stats[k].sum += new_results->stats[k].sum;
        results->stats[k].square_sum += new_results->stats[k].square_sum;
    }
}


void results_add_photons_received(Results * res, PhotonBlock * block)
{
    Photon * photons;
    Stat * st;
    float weight;
    float weight_th;
    float weight_tof;
    float weight_dist;
    float weight_ev;
    float radius_sqr;
    int k;

    photons = block->photons;
    st = res->stats;

    for(k = 0; k < block->length; ++k)
    {
        if(photons[k].received)
        {
            res->photons_received++;

            weight = photons[k].weight;
            st[WEIGHT].sum += (double)weight;
            st[WEIGHT].square_sum += (double)(weight * weight);
            weight_th = weight * photons[k].theta;
            st[THETA].sum += (double)(weight_th);
            st[THETA].square_sum += (double)(weight_th * photons[k].theta);
            weight_tof = weight * photons[k].tof;
            st[TOF].sum += (double)(weight_tof);
            st[TOF].square_sum += (double)(weight_tof * photons[k].tof);
            weight_dist = weight * photons[k].distance;
            st[DISTANCE].sum += (double)(weight_dist);
            st[DISTANCE].square_sum += (double)(weight_dist *
                                       photons[k].distance);
            radius_sqr = photons[k].x * photons[k].x +
                         photons[k].y * photons[k].y;
            st[RADIUS].sum += (double)(weight * sqrtf(radius_sqr));
            st[RADIUS].square_sum += (double)(weight * radius_sqr);
            weight_ev = weight * photons[k].events;
            st[EVENTS].sum += (double)(weight_ev);
            st[EVENTS].square_sum += (double)(weight_ev * photons[k].events);
        }
    }
}


void results_set_events_threshold(Results * results, long events_threshold)
{
    results->events_threshold = events_threshold;
}


void add_processing_time_to_json_object(Results * results, unsigned num_results,
                                        JsonObject * jobject)
{
    JsonArray * jarray;
    unsigned k;

    if(num_results > 1)
    {
        jarray = json_array_new();
        json_object_set_array_member(jobject, "processing time", jarray);
        for(k = 0; k < num_results; ++k)
            json_array_add_double_element(jarray, results[k].processing_time);
    }
    else
    {
        json_object_set_double_member(jobject,"processing time",
                                      results->processing_time);
    }

}


void add_photons_received_to_json_object(Results * results,
                                         unsigned num_results,
                                         JsonObject * jobject)
{
    JsonArray * jarray;
    unsigned k;

    if(num_results > 1)
    {
        jarray = json_array_new();
        json_object_set_array_member(jobject, "photons received", jarray);
        for(k = 0; k < num_results; ++k)
            json_array_add_int_element(jarray, results[k].photons_received);
    }
    else
    {
        json_object_set_int_member(jobject,"photons received",
                                      results->photons_received);
    }

}


void add_events_threshold_to_json_object(Results * results,
                                         unsigned num_results,
                                         JsonObject * jobject)
{
    JsonArray * jarray;
    unsigned k;

    if(num_results > 1)
    {
        jarray = json_array_new();
        json_object_set_array_member(jobject, "events threshold", jarray);
        for(k = 0; k < num_results; ++k)
            json_array_add_int_element(jarray, results[k].events_threshold);
    }
    else
    {
        json_object_set_int_member(jobject,"events threshold",
                                      results->events_threshold);
    }

}


void add_stat_to_json_object(Results * results, unsigned num_results,
                             unsigned element, JsonObject * jobject)
{
    JsonObject * jstat;
    JsonArray * jsum;
    JsonArray * jsquares_sum;
    unsigned k;

    jstat = json_object_new();

    json_object_set_object_member(jobject, stats_names[element], jstat);

    if(num_results > 1)
    {
        jsum = json_array_new();
        jsquares_sum = json_array_new();
        json_object_set_array_member(jstat, "sum", jsum);
        json_object_set_array_member(jstat, "square sum", jsquares_sum);
        for(k = 0; k < num_results; ++k)
        {
            json_array_add_double_element(jsum, results[k].stats[element].sum);
            json_array_add_double_element(jsquares_sum,
                                          results[k].stats[element].square_sum);
        }
    }
    else
    {
        json_object_set_double_member(jstat, "sum",
                                      results->stats[element].sum);
        json_object_set_double_member(jstat, "square sum",
                                      results->stats[element].square_sum);
    }
}


JsonObject * results_get_json_object(Results * results, unsigned num_results)
{
    JsonObject * jresults;
    JsonObject * jstats;
    unsigned k;

    jresults = json_object_new();
    jstats = json_object_new();
/*
    if(stats_ballistic_time == INFINITY)
        json_object_set_null_member(jresults, "ballistic time");
    else
        json_object_set_double_member(jresults, "ballistic time",
                                      (double)stats_ballistic_time);
*/
    add_photons_received_to_json_object(results, num_results, jresults);
    add_events_threshold_to_json_object(results, num_results, jresults);
    json_object_set_object_member(jresults, "stats", jstats);
    for(k = 0; k < NUMBER_OF_STATS; ++k)
        add_stat_to_json_object(results, num_results, k, jstats);
    add_processing_time_to_json_object(results, num_results, jresults);

    return jresults;
}





