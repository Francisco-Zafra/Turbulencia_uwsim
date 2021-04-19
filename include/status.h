#ifndef __STATUS_H__
#define __STATUS_H__


#include "results.h"


void start_status_tracking(long photons_required_);
void end_status_tracking();
void update_status(long photons_processed, long photons_received);

void start_processing_time_tracking();
void set_processing_time(Results * results);


#endif