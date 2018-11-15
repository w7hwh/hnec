using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    /// <summary>
    /// class to contain all output data that NEC produces for one frequency
    /// </summary>
    class NecOutput
    {
        double f;          // frequency in MHz
                           //Radpattern* rp;    // pointer to radiation pattern data, if available; NULL otherwise
                           //public fixed double d[(int)NECO.neco_numdat];
        int rpgpovalid;    /* flag: !=0 if the gpo field of *rp already contains valid data for the present settings */
                           //Currents* cu;      /* pointer to currents data, if available; NULL otherwise */
                           //Nearfield* nf;     /* pointer to nearfield data, if available; NULL otherwise */
        double maxe, maxh;  /* largest values of E and H in nearfield data */

  
        /// <summary>
        /// indices in the element 'd' of the NECoutput structure defined below
        /// </summary>
        public enum NECO
        {
            neco_zr, neco_zi,       /* index of real and imaginary part of impedance */
            neco_zphi, neco_zabs,   /* phi and absolute value of impedance */
            neco_swr,               /* index of SWR */
            neco_maxgain, neco_fb,  /* index of maximum gain and corresponding front/back ratio */
            neco_phi, neco_theta,   /* index direction of maximum gain */
            neco_vgain, neco_vfb,   /* index of gain in viewing direction and corresponding f/b ratio */
            neco_vgain2, neco_vfb2, /* same, but taking polarization into account; fb/fb2 as fb1/fb2 below */
            neco_maxgain_hor, neco_fb_hor1, neco_fb_hor2, neco_phi_hor, neco_theta_hor,
            neco_maxgain_vert, neco_fb_vert1, neco_fb_vert2, neco_phi_vert, neco_theta_vert,
            neco_maxgain_lhcp, neco_fb_lhcp1, neco_fb_lhcp2, neco_phi_lhcp, neco_theta_lhcp,
            neco_maxgain_rhcp, neco_fb_rhcp1, neco_fb_rhcp2, neco_phi_rhcp, neco_theta_rhcp,
            neco_numdat             /* last element of the enum: number of elements in the array */
        };

        /* since (maximum) gain and corresponding quantities are not available in all NECoutput records,
           the following macro is needed to tell whether a quantity is available always, or only if
           radiation pattern data is available in the record: */
        //#define ONLY_IF_RP(i) ((i)>=neco_maxgain)
        /* a constant for more systematic access to the polarization-specific fields: */
        //#define Neco_gsize 5     /* group size */
        //#define Neco_polgain (neco_maxgain_hor-Neco_gsize)
        //#define Neco_polfb1 (neco_fb_hor1-Neco_gsize)
        //#define Neco_polfb2 (neco_fb_hor2-Neco_gsize)
        //#define Neco_polphi (neco_phi_hor-Neco_gsize)
        //#define Neco_poltheta (neco_theta_hor-Neco_gsize)
        /* now, we can write say  neco_phi_lhcp  as  Neco_gsize*POLlhcp+Neco_polphi  */
        /* note: fb...1 is f/b ratio for this polarization; fb...2 is same but for back total power is used */

        //extern NECoutput * neco;
        //extern int numneco;
        //extern int maxfreqs;

        //extern double globalmaxdb; /* maximum gain, global over all frequencies and all output data files */

        //extern int rp_index;       /* index of the entry in neco[] whose radiation pattern is being shown in the 3D plot */
    }
}
