// $Id: Panyimain.cpp,v 1.1.1.1 2016/12/30 06:03:09 zjcao Exp $
// Ref PRD 96, 044028 (2017) 
#ifdef newc
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>
#else
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#endif
#include "DIY_SEOBNRE_ALL.h"


int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in );
int output_waveform(REAL8TimeSeries * h_plus, REAL8TimeSeries * h_cross, int ampPhase);

const char * usage =
"Generate a simulation using the lalsimulation library\n\n"
"The following options can be given (will assume a default value if omitted):\n"
//"--domain DOM               'TD' for time domain (default) or 'FD' for frequency\n"
//"                           domain; not all approximants support all domains\n"
"--amp-phase                If given, will output:\n"
"                           |h+ - i hx|, Arg(h+ - i hx) (TD) or\n"
"                           |h+(f)|, Arg(h+(f)), |hx(f)|, Arg(hx(f)) (FD)\n"
"                           If not given, will output h+ and hx (TD and FD)\n"

"                           NOTE: Other approximants may be available if the\n"
"                           developer forgot to edit this help message\n"
//"--phase-order ORD          Twice PN order of phase (default ORD=7 <==> 3.5PN)\n"
"--amp-order ORD            Twice PN order of amplitude (default 0 <==> Newt.)\n"
"--phiRef phiRef            Phase at the reference frequency (default 0)\n"
//"--fRef FREF                Reference frequency in Hz\n"
"                           (default: 0)\n"
"--sample-rate SRATE        Sampling rate of TD approximants in Hz (default 4096)\n"
//"--deltaF DF                Frequency bin size for FD approximants in Hz (default 1/8)\n"
"--m1 M1                    Mass of the 1st object in solar masses (default 10)\n"
"--m2 M2                    Mass of the 2nd object in solar masses (default 1.4)\n"
"--inclination IOTA         Angle in radians between line of sight (N) and \n"
"                           orbital angular momentum (L) at the reference\n"
"                           (default 0, face on)\n"
"--spin1x S1X               Vector components for spin of mass1 (default all 0)\n"
"--spin1y S1Y               z-axis=line of sight, L in x-z plane at reference\n"
"--spin1z S1Z               Kerr limit: s1x^2 + s1y^2 + s1z^2 <= 1\n"
"--spin2x S2X               Vector components for spin of mass2 (default all 0)\n"
"--spin2y S2Y               z-axis=line of sight, L in x-z plane at reference\n"
"--spin2z S2Z               Kerr limit: s2x^2 + s2y^2 + s2z^2 <= 1\n"
//"--tidal-lambda1 L1         (tidal deformability of mass 1) / (mass of body 1)^5\n"
//"                           (~128-2560 for NS, 0 for BH) (default 0)\n"
//"--tidal-lambda2 L2         (tidal deformability of mass 2) / (mass of body 2)^5\n"
//"                           (~128-2560 for NS, 0 for BH) (default 0)\n"
//"--spin-order ORD           Twice PN order of spin effects\n"
//"                           (default ORD=-1 <==> All spin effects)\n"
//"--tidal-order ORD          Twice PN order of tidal effects\n"
//"                           (default ORD=-1 <==> All tidal effects)\n"
"--f-min FMIN               Lower frequency to start waveform in Hz (default 40)\n"
"--f-max FMAX               Frequency at which to stop waveform in Hz\n"
"                           (default: generate as much as possible)\n"
"--e0 e0                    eccentricity at FMIN (default 0)\n"
"--distance D               Distance in Mpc (default 100)\n"
//"--axis AXIS                for PhenSpin: 'View' (default), 'TotalJ', 'OrbitalL'\n"
"--longAscNodes             Longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation\n"
"--output                   Do not out put to file\n"
"--outname FNAME            Output to file FNAME (default 'simulation.dat')\n"
"--verbose                  If included, add verbose output\n"
"--jobtag                   If not equal to \"None\", print dynamics as XX_jobtag.txt\n"
"                           (default None: no output)\n"
"                           Formation: \"time, omega, r, pr, phi, pphi \"\n"
;

/* Parse command line, sanity check arguments, and return a newly
 * allocated GSParams object */

GSParams *parse_args(unsigned int argc, char **argv) {
    unsigned int i;
    GSParams *params;
    params = (GSParams *) malloc(sizeof(GSParams));
    memset(params, 0, sizeof(GSParams));

    /* Set default values to the arguments */
    params->ampO = 0;
    params->phiRef = 0.;
    params->deltaT = 1./4096.;
    params->m1 = 10 * LAL_MSUN_SI;
    params->m2 = 1.4 * LAL_MSUN_SI;
    params->f_min = 40;
    params->f_max = 0.; /* Generate as much as possible */
    params->e0 = 0.;
    params->distance = 100. * 1e6 * LAL_PC_SI;
    params->inclination = 0.;
    params->s1x = 0.;
    params->s1y = 0.;
    params->s1z = 0.;
    params->s2x = 0.;
    params->s2y = 0.;
    params->s2z = 0.;
    params->longAscNodes = 0.;
    params->verbose = 0; /* No verbosity */
    params->output = 1;
    strncpy(params->outname, "simulation.dat", 256); /* output to this file */
    params->ampPhase = 0; /* output h+ and hx */
    strncpy(params->jobtag, "None", 256);

    /* consume command line */
    for (i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            goto fail;
        
        } else if (strcmp(argv[i], "--amp-order") == 0) {
            params->ampO = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--sample-rate") == 0) {
            params->deltaT = 1./atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]) * LAL_MSUN_SI;
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]) * LAL_MSUN_SI;
        } else if (strcmp(argv[i], "--spin1x") == 0) {
            params->s1x = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin1y") == 0) {
            params->s1y = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin1z") == 0) {
            params->s1z = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2x") == 0) {
            params->s2x = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2y") == 0) {
            params->s2y = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2z") == 0) {
            params->s2z = atof(argv[++i]);
        } else if (strcmp(argv[i], "--f-min") == 0) {
            params->f_min = atof(argv[++i]);
            fprintf(stderr,"f_0 = %eHz (%e[1/M])\n", params->f_min,params->f_min*((params->m1+params->m2)/LAL_MSUN_SI*LAL_MTSUN_SI));
        } else if (strcmp(argv[i], "--f-max") == 0) {
            params->f_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--e0") == 0) {
            params->e0 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]) * 1e6 * LAL_PC_SI;
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 ) {
            params->output = 0;
        } else if (strcmp(argv[i], "--outname") == 0) {
            strncpy(params->outname, argv[++i], 256);
        } else if (strcmp(argv[i], "--amp-phase") == 0) {
            params->ampPhase = 1;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            params->verbose = 1;
        } else if (strcmp(argv[i], "--jobtag") == 0) {
            strncpy(params->jobtag, argv[++i], 256);
        } else if (strcmp(argv[i], "--longAscNodes") == 0) {
            params->longAscNodes = atof(argv[++i]);
        } else {
            fprintf(stderr,"Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    //printf("M=m1+m2=%e(s), M/R=%e\n",(params->m1+params->m2)/LAL_MSUN_SI*LAL_MTSUN_SI,(params->m1+params->m2)/LAL_MSUN_SI*LAL_MRSUN_SI/params->distance);

    return params;

    fail:
    fprintf(stderr,"%s", usage);
    exit(1);
}
int dump_TD(FILE *f, REAL8TimeSeries *hplus, REAL8TimeSeries *hcross) {
    size_t i;
    REAL8 t0 = XLALGPSGetREAL8(&(hplus->epoch));
    if (hplus->data->length != hcross->data->length) {
        fprintf(stderr,"Error: hplus and hcross are not the same length\n");
        return 1;
    } else if (hplus->deltaT != hcross->deltaT) {
        fprintf(stderr,"Error: hplus and hcross do not have the same sample rate\n");
        return 1;
    }

    fprintf(f, "# t hplus hcross\n");
    for (i=0; i < hplus->data->length; i++)
        fprintf(f, "%.16e %.16e %.16e\n", t0 + i * hplus->deltaT, 
                hplus->data->data[i], hcross->data->data[i]);

    return 0;
}

int dump_convertTD(GSParams *params, REAL8TimeSeries *hplus, REAL8TimeSeries *hcross) {
    size_t i;
    if (hplus->data->length != hcross->data->length) {
        fprintf(stderr, "Error: hplus and hcross are not the same length\n");
        return 1;
    } else if (hplus->deltaT != hcross->deltaT) {
        fprintf(stderr,"Error: hplus and hcross do not have the same sample rate\n");
        return 1;
    }

    const REAL8 UM = (params->m1+params->m2)/LAL_MSUN_SI*LAL_MTSUN_SI;
    const REAL8 UMoR = (params->m1+params->m2)/LAL_MSUN_SI*LAL_MRSUN_SI/params->distance;

    FILE *f;
    f = fopen(params->outname, "w");
    fprintf(f, "# t hplus hcross\n");
    for (i=0; i < hplus->data->length; i++)
        // convert to [M] for t and [M/R] for h
        fprintf(f, "%.16e %.16e %.16e\n", (i * hplus->deltaT)/UM, 
                (hplus->data->data[i])/UMoR, (hcross->data->data[i])/UMoR);
    
    fclose(f);

    return 0;
}

void genwaveform (double * outhp,double * outhc, double * t0, int * datalength,\
    double phiRef, double deltaT, double m1, double m2, double s1z, double s2z,\
    double f_min, double e0, double distance, double inclination, double long_asc_nodes) {

    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    int i=0;
    
    m1 = m1 * LAL_MSUN_SI;
    m2 = m2 * LAL_MSUN_SI;
    distance = distance * 1.e6 * LAL_PC_SI;
    // in Panyi_elip.cpp 
    XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phiRef, 
                    deltaT, m1, m2, 0, 
                    0, s1z, 0, 0,  
                    s2z, f_min, e0, 
                    distance, inclination, long_asc_nodes, "None"
                    );
    //printf("%d\n",hplus->data->length);

    for (i=0;i<hplus->data->length;i++){
        *(outhp+i) = hplus->data->data[i];
        *(outhc+i) = hcross->data->data[i];
    }

    *(datalength) = hplus->data->length;
    if(*(datalength)>100000)
    fprintf(stderr,"ERROR: Data length exceeds the maximum allowed size\n");

    * (t0) = XLALGPSGetREAL8(&(hplus->epoch));
    //printf("t0 is %e, hplus epoch is %e, \n",* t0, XLALGPSGetREAL8(&(hplus->epoch)));

}

int output_waveform(REAL8TimeSeries * h_plus, REAL8TimeSeries * h_cross, int ampPhase)
{
    double t0, t_end;
    int j;
    t0 = XLALGPSGetREAL8(&(h_plus->epoch));
    t_end = t0 + (h_plus->data->length - 1) * h_plus->deltaT;
    if (ampPhase)
    {
        REAL8Vector *amp;
        REAL8Vector *phi;
        amp->length = h_plus->data->length;
        phi->length = h_plus->data->length;
        double phi0;
        
        for (j = 0; j < h_plus->data->length; ++j)
        {
            double complex z = h_plus->data->data[j] - I * h_cross->data->data[j];
            amp->data[j] = cabs(z);
            phi->data[j] = carg(z);
        }
        
        XLALREAL8VectorUnwrapAngle(phi, phi);
        
        phi0 = 2 * phi->data[phi->length - 1] - phi->data[phi->length - 2];
        phi0 -= fmod(phi0 + copysign(LAL_PI, phi0), 2.0 * LAL_PI) - copysign(LAL_PI, phi0);
        for ( j = 0; j < phi->length; ++j )
            phi->data[j] -= phi0;
        
        fprintf(stdout, " # time (s)\th_abs (strain)\t h_arg (rad)\n ");
        for ( j=0; j < h_plus->data->length; ++j )
            fprintf(stdout, "%.9f\t%e\t%e\n", t0 + j * h_plus->deltaT, amp->data[j], phi->data[j] );
        
        XLALDestroyREAL8Vector(phi);
        XLALDestroyREAL8Vector(amp);
    } else
    {
        fprintf(stdout, "# time (s)\th_+ (strain)\th_x (strain)\n" );
        for ( j=0; j < h_plus->data->length; ++j )
            fprintf(stdout, "%.9f\t%e\t%e\n" , t0 + j * h_plus->deltaT, h_plus->data->data[j], h_cross->data->data[j] );
    }
    return 0;
}

int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in )
{
    REAL8 prev;
    REAL8 diff;
    INT4  wrap;
    UINT4 i;
    if ( ! out || ! in )
        return XLAL_EFAULT;
    if ( ! out->data || ! in->data || in->length == 0 )
        return XLAL_EINVAL;
    if ( out->length != in->length )
        return XLAL_EBADLEN;
    wrap = 0;
    prev = out->data[0] = in->data[0];
    for ( i = 1; i < in->length; ++i ) {
        diff  = in->data[i] - prev;
        prev  = in->data[i];
        wrap += (diff < -LAL_PI) - (diff > LAL_PI);
        out->data[i] = in->data[i] + wrap * LAL_TWOPI;
    }
    return 0;
}

int XLALSimInspiralChooseTDWaveform(
                                    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
                                    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
                                    REAL8 phiRef,                               /**< reference orbital phase (rad) */
                                    REAL8 deltaT,                               /**< sampling interval (s) */
                                    REAL8 m1,                                   /**< mass of companion 1 (kg) */
                                    REAL8 m2,                                   /**< mass of companion 2 (kg) */
                                    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
                                    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
                                    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
                                    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
                                    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
                                    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
                                    REAL8 f_min,                                /**< starting GW frequency (Hz) */
                                    REAL8 e0,                                   /**< eccentricity at starting GW frequency (Hz) */
                                    REAL8 r,                                    /**< distance of source (m) */
                                    REAL8 i,                                    /**< inclination of source (rad) */
                                    REAL8 longAscNodes,                   /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
                                    char *jobtag
)
{
    
    int ret;
    /* N.B. the quadrupole of a spinning compact body labeled by A is
     * Q_A = - quadparam_A chi_A^2 m_A^3 (see gr-qc/9709032)
     * where quadparam = 1 for BH ~= 4-8 for NS.
     * This affects the quadrupole-monopole interaction.
     * For now, hardcode quadparam1,2 = 1.
     * Will later add ability to set via LALSimInspiralTestGRParam
     */
    REAL8 v0 = 1., quadparam1 = 1., quadparam2 = 1.;
    
    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    
    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        fprintf(stderr,"XLAL Warning - : Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", deltaT);
    if( deltaT < 1./16385. )
        fprintf(stderr,"XLAL Warning - : Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        fprintf(stderr,"XLAL Warning - : Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        fprintf(stderr,"XLAL Warning - : Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        fprintf(stderr,"XLAL Warning - : Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( S1x*S1x + S1y*S1y + S1z*S1z > 1.000001 )
        fprintf(stderr,"XLAL Warning - : S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", S1x, S1y, S1z);
    if( S2x*S2x + S2y*S2y + S2z*S2z > 1.000001 )
        fprintf(stderr,"XLAL Warning - : S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", S2x, S2y, S2z);
    if( f_min < 1. )
        fprintf(stderr,"XLAL Warning - : Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n",f_min);
    if( f_min > 40.000001 )
        fprintf(stderr,"XLAL Warning - : Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", f_min);
    if( e0 < 0 || e0 > 1 )
        fprintf(stderr,"XLAL Warning - : unphysical e0 = %e.\n", e0);
    
    
    REAL8 polariz=longAscNodes;
    //polariz+=-LAL_PI/2.;
    

    /* Call the waveform driver routine */
    ret = XLALSimSEOBNRE(hplus, hcross, phiRef,
                         deltaT, m1, m2, f_min, e0, r, i, S1z, S2z, jobtag); // in current file
    
    
    if (polariz && (*hplus) && (*hcross) ) {
      REAL8 tmpP,tmpC;
      REAL8 cp=cos(2.*polariz);
      REAL8 sp=sin(2.*polariz);
      for (UINT4 idx=0;idx<(*hplus)->data->length;idx++) {
        tmpP=(*hplus)->data->data[idx];
        tmpC=(*hcross)->data->data[idx];
        (*hplus)->data->data[idx] =cp*tmpP+sp*tmpC;
        (*hcross)->data->data[idx]=cp*tmpC-sp*tmpP;
      }
    }
    
    return ret;
}

int main (int argc , char **argv) {

    FILE *f;

    int start_time;
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    GSParams *params;

    /* parse commandline */
    params = parse_args(argc, argv);
    
    /* generate waveform */
    start_time = time(NULL);

    // in Panyi_elip.cpp 
    XLALSimInspiralChooseTDWaveform(&hplus, &hcross, params->phiRef, 
                    params->deltaT, params->m1, params->m2, params->s1x, 
                    params->s1y, params->s1z, params->s2x, params->s2y, 
                    params->s2z, params->f_min, params->e0, 
                    params->distance, params->inclination, params->longAscNodes, params->jobtag
                    );
          
    if (params->verbose)
        fprintf(stderr,"Generation took %.0f seconds\n",
                difftime(time(NULL), start_time));

    /* dump file */
    if(params->output)
    {
      f = fopen(params->outname, "w");
      dump_TD(f, hplus, hcross);
      fclose(f);
    }
    else
    {
        output_waveform(hplus, hcross, params->ampPhase);
    }

    return 1;
}