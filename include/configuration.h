#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__


#define PROGRAM_VERSION                 "uwsim 3.0.2"
#define PROGRAM_BUG_ADDRESS             "pedroricopinazo@gmail.com"


#define INF                             INFINITY

#define DEFAULT_ROULETTING_THRESHOLD    (1e-6)
#define ROULETTING_SURVIVE_PROBABILITY  (1.0f/10.0f)
#define ROULETTING_WEIGHT_INCREASE      (1.0f/ROULETTING_SURVIVE_PROBABILITY)
#define ROULETTING_MAX_VALUE            ((unsigned)ROULETTING_WEIGHT_INCREASE)

#define DEFAULT_NUM_PHOTONS             (0)

//#define CDF_SIZE                      (100 * 10  + 10 * 170) // = 2700
//#define PHOTON_BLOCK_SIZE               (1024 * 1024)

#define COS_THETA_LUT_SIZE              (100000)

#define SPEED_VACUUM                    (300000000.0f) //(299792458.0f)

// File extensions used
#define RD_EXTENSION                    (".dat")    // Readable extension
#define BN_EXTENSION                    (".rws")    // Binary extension
#define SW_EXTENSION                    (".json")   // Sweep extension

// Option descriptions
#define DOC_P   "set phase function to <arg>"
#define DOC_g   "set g parameter of phase functions to <value>"
#define DOC_i   "set input file to <arg>"
#define DOC_r   "set output to readable format"
#define DOC_z   "set z position of the receptor to <value>"
#define DOC_x   "set x position of the receptor to <value>"
#define DOC_y   "set y position of the receptor to <value>"
#define DOC_c   "set attenuation of the medium to <value>"
#define DOC_w   "set albedo of the medium to <value>"
#define DOC_F   "set full angle of the source to <value> (degrees)"
#define DOC_b   "set beam waist of the source to <value> (mm)"
#define DOC_d   "set divergence of the source to <value> (mrad)"
#define DOC_a   "set aperture of the receptor to <value> (cm)"
#define DOC_f   "set field-of-view of the receptor to <value> (degrees)"
#define DOC_701 "set rouletting threshold to <value>"
#define DOC_n   "set number of photons to <value>"
#define DOC_t   "set time of simulation to <value> (ns)"
#define DOC_T   "set source type to <arg>"
#define DOC_C   "Continue simulation contained in <arg>"
#define DOC_o   "set output file to <arg>"
#define DOC_e   "Set theta angle of the receptor to <value> (degrees)"
#define DOC_h   "Set phi angle of the receptor to <value> (degrees)"
#define DOC_E   "Set theta angle of the source to <value> (degrees)"
#define DOC_H   "Set phi angle of the source to <value> (degrees)"
#define DOC_D   "Set the sea surface y position to <value> (meters)"
#define DOC_s    "Set standard deviation of theta angle of the receptor to <value> (degrees)"
#define DOC_W   "Set standard deviation of surface inclination to <value> (degrees)"
#define DOC_700 "Set water refractive index to <value>"
#define DOC_m   "Set maximun number of events to <value>"
#define DOC_M   "Set reception mode to <arg>"
#define DOC_X   "Set receptor offset in x axis of the receptor plane"
#define DOC_Y   "Set receptor offset in y axis of the receptor plane"
#define DOC_702 "Set air refractive index outside the surface to <value>"
#define DOC_703 "Set air refractive index inside the detector"
#define DOC_S   "Set standard deviation of theta angle of the source to <value> (degrees)"
#define DOC_704 "Set photon block size for parallelization to <value>"
#define DOC_l   "Set length factor, l,  to <value>. This factor set the events threshold to l times the scattering length. The scattering length is the average of the number of events suffered by a photon when this travels the link length"
#define DOC_q   "Hide configuration and results output"
#define DOC_L   "Set number of layers to <value>" //número de capas

// Option IDs
#define ID_P    0x01
#define ID_g    0x02
#define ID_z    0x03
#define ID_x    0x04
#define ID_y    0x05
#define ID_c    0x06
#define ID_w    0x07
#define ID_b    0x08
#define ID_d    0x09
#define ID_a    0x0A
#define ID_f    0x0B
#define ID_701  0x0C
#define ID_n    0x0D
#define ID_T    0x0E
#define ID_F    0x0F
#define ID_e    0x10
#define ID_h    0x11
#define ID_E    0x12
#define ID_H    0x13
#define ID_D    0x14
#define ID_s    0x15
#define ID_W    0x16
#define ID_700  0x17
#define ID_m    0x18
#define ID_M    0x19
#define ID_X    0x1A
#define ID_Y    0x1B
#define ID_702  0x1C
#define ID_703  0x1D
#define ID_S    0x1E
#define ID_l    0x1F
#define ID_L    0x20    //identificador numero de capas


#endif