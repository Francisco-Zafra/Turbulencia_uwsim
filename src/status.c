#include "status.h"

#include <stdbool.h>
#include <time.h>
#include <stdio.h>


static double start_time_status;
static double photons_required;

static double processing_time_start;


void split_time(time_t total_seconds, unsigned * seconds,
                unsigned * minutes, unsigned * hours)
{
    *seconds = total_seconds % 60;
    *hours = total_seconds / 3600;
    *minutes = (total_seconds % 3600) / 60;
}


double get_monotonic_time()
{
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &current_time);
    return (double)current_time.tv_sec + 1e-9f * (double)current_time.tv_nsec;
}


void start_status_tracking(long photons_required_)
{
    start_time_status = get_monotonic_time();
    photons_required = (double)photons_required_;
}


void end_status_tracking()
{
    // unsigned k;
    // unsigned length;

    // length = 70;

    // fprintf(stderr, "\r");
    // for(k = 0; k < length; ++k)
    //     fprintf(stderr, " ");
    // fprintf(stderr, "\r");
    fprintf(stderr, "\n");
}


void update_status(long photons_processed, long photons_received)
{
    static double total_photons_processed = 0;
    unsigned seconds, minutes, hours;
    double photons_remaining;
    double time_spent;
    double time_spent_per_photon;
    double progress;
    time_t time_remaining;

    total_photons_processed += photons_processed;

    progress = (double)total_photons_processed / photons_required;

    photons_remaining = photons_required - total_photons_processed;

    time_spent = get_monotonic_time() - start_time_status;
    time_spent_per_photon = time_spent / total_photons_processed;
    time_remaining = (time_t)(time_spent_per_photon * photons_remaining);

    split_time(time_remaining, &seconds, &minutes, &hours);

    fprintf(stderr, "\r");   // Return to the first position in the line
    fprintf(stderr, "\033[1;36m");   // Set output color to green
    fprintf(stderr, "Progress: %.2f%%  -  ", progress * 100.0f);
    fprintf(stderr, "%02ih %02im %02is left", hours, minutes, seconds);
    if(photons_received != 0)
        fprintf(stderr, "  -  %.1e photons received", (float)photons_received);
    fprintf(stderr, "\033[0m");   // Reset output color to default
    fflush(stdout);
}


void start_processing_time_tracking()
{
    processing_time_start = get_monotonic_time();
}


void set_processing_time(Results * results)
{
    results->processing_time = get_monotonic_time() - processing_time_start;
}