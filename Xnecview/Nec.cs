/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 1998-2006,2011, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 *
 *  Main module, which mostly deals with initial setup, (re)reading of
 *  files, and parsing command line options.
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// The static Nec class holds all the information from the NEC input 
/// and/or output file(s) that's needed by multiple classes.
/// It also holds the "command line" options for the run.
/// </summary>
namespace Necview
{
    public static class Nec
    {
        #region Point
        /// <summary>
        /// class to describe a point in 3D space
        /// </summary>
        public class Point
        {
            public double x;
            public double y;
            public double z;

            public override string ToString() => string.Format("Point({0},{1},{2})", x, y, z);
        };
        #endregion
        #region Wire
        /// <summary>
        /// class to describe a straight wire (segment)
        /// </summary>
        public class Wire
        {
            public Point p0;          /* begin point */
            public Point p1;          /* end point */
            public double rad;        /* radius */
            public int itg;           /* tag number */
            public int ns;            /* number of segments for multi-segment wire; -segmentnumber for single-segment wire that's part of a larger wire (GA and GH cards)  */
            public int abs;           /* absolute segment number of the first segment of this wire */
            public bool dispnr;       /* flag: true on any multi-segment wire, and on (a) middle segment for a wire that was split into multiple single-segment wires (GA and GH) */

            public Wire()
            {
                this.p0 = new Nec.Point();
                this.p1 = new Nec.Point();
            }

            public override string ToString() => string.Format("Wire = {0}, {1}, {2}, {3}, {4}, {5}, {6}", itg, ns, p0.ToString(), p1.ToString(), rad, abs, dispnr);
        };

        public static List<Wire> wires = new List<Wire>(25);

        public static void DumpWires()
        {
            Console.WriteLine("------- DumpWires() dumping {0} wires -------", wires.Count);
            foreach (Wire w in wires)
            {
                Console.WriteLine(w.ToString());
            }
            Console.WriteLine("----------------------------------------");
        }
        #endregion
        #region Surface
        /// <summary>
        /// class to describe a surface patch
        /// </summary>
        public class Surface
        {
            public SU type;     /* see below */
            public Point p0;    /* corner */
            public Point p1;    /* corner */
            public Point p2;    /* corner */
            public Point p3;    /* corner */
            public Point pc;    /* center */
            public Point pa;    /* tip of outward normal arrow */

            public Surface()
            {
                this.p0 = new Nec.Point();
                this.p1 = new Nec.Point();
                this.p2 = new Nec.Point();
                this.p3 = new Nec.Point();
                this.pc = new Nec.Point();
                this.pa = new Nec.Point();
            }

            public override string ToString() => string.Format("Surface = {0}, {1}, {2}, {3}, {4}, {5}, {6}", 
                        type, p0.ToString(), p1.ToString(), p2.ToString(), p3.ToString(), pc.ToString(), pa.ToString());
        };

        public enum SU
        {
            SU_rect = 1,    /* rectangular shape (p2,p3 ignored) */
            SU_arb = 2,     /* "arbitrary" shape, shown as a set of 8 crossing lines (p2,p3 ignored) */
            SU_tri = 3,     /* triangle (p3 ignored) */
            SU_quad = 4     /* quadrilateral shape */
        }

        public static List<Surface> surfaces = new List<Surface>();

        public static void DumpSurfaces()
        {
            Console.WriteLine("------- DumpSurfaces() dumping {0} surfaces -------", surfaces.Count);
            foreach (Surface s in surfaces)
            {
                Console.WriteLine(s.ToString());
            }
            Console.WriteLine("----------------------------------------");
        }
        #endregion
        #region Exci
        /// <summary>
        /// class to describe an excitation (only voltage sources, EX card)
        /// </summary>
        public class Exci
        {
            // TODO fix exci
            //Wire* wire;        /* pointer to Wire data */
            public int seg;           /* segment number within wire */

            public override string ToString() => string.Format("Exci = {0}", seg);
        };

        public static List<Exci> excis = new List<Exci>();

        public static void DumpExci()
        {
            Console.WriteLine("------- DumpExci() dumping {0} excis -------", excis.Count);
            foreach (Exci e in excis)
            {
                Console.WriteLine(e.ToString());
            }
            Console.WriteLine("----------------------------------------");
        }
        #endregion
        #region Load
        /* class to describe a load */
        public class Load
        {
            //Wire* wire;        /* pointer to Wire data */
            public int firstseg;      /* number of first loaded segment within wire */
            public int lastseg;       /* number of last loaded segment within wire */
        };               /* note: if an LD card specifies a loading which extends over several wires, it is represented here by one Load struct per wire */

        public static int numloads;
        public static int maxloads;
        public static List<Load> loads = new List<Load>();
        #endregion
        #region Netw
        /* class to describe a "network" in the antenna structure; 
         * a transmission line is a special case of this */
        public class Netw
        {
            //Wire* wire1;       /* pointer to Wire data of one end of the network / transmission line */
            public int seg1;          /* segment number within wire */
            //Wire* wire2;       /* same for the other end of the network / transmission line */
            public int seg2;
            public NETW_TYPE type;

            public override string ToString() => string.Format("Netw = {0}, {1}, {2}", seg1, seg2, type);
            };

        public enum NETW_TYPE
        {
            generic = 0,            // generic network
            noPhaseReversal = 1,    // transmission line without phase reversal
            phaseReversal = -1      // transmission line with phase reversal
        }
        public static List<Netw> netws = new List<Netw>();

        public static void DumpNetws()
        {
            Console.WriteLine("------- DumpNetws() dumping {0} netws -------", netws.Count);
            foreach (Netw n in netws)
            {
                Console.WriteLine(n.ToString());
            }
            Console.WriteLine("----------------------------------------");
        }
        #endregion

        public static double extr_str;	/* size of the biggest dimension of the structure data to be shown; used for (initial) scaling */
        public static double extr;	/* size of the biggest dimension of all data to be shown; used for (initial) scaling */

        public static int quick;         /* flag to indicate that drawing should be rough but quick */
        public static int redraw;        /* flag which signifies need for redrawing */
        public static int dragging;      /* flag to indicate that user is dragging the structure */

        /* these define the projection: */
        public static double phi, theta;      /* angles defining direction of eye */
        public static double zoom;           /* zoom factor */
        //public double trx,try;        /* 2D-translation, as a fraction of winsize */
        public static int winsizex, winsizey; /* size of window in pixels */
        public static int winsize;          /* min(winsizex,winsizey) */

        public static GS gainscale;	     /* type of scaling of gain: */
        /// <summary>
        /// possible values for gainscale
        /// </summary>
        public enum GS
        {
            GSlinpower,
            GSlinvolt,
            GSarrl,
            GSlog,
            numGS
        };
        public static string[] GSnames;         /* array with names of gainscales */
        public static string[][] GSnumbers;     /* array with descriptions of levels of scale lines in plot */

        public static GP gainplot;         /* type of radiation plot: */
        public enum GP
        {
            GPnone,
            GPslice,
            GPframe,
            GPopaque,
            GPnearfield,
            numGP           // TODO remove
        };

        public static SP structplot;       /* type of structure plot: */
        public enum SP
        {
            SPnone,
            SPstruct,
            SPtags,
            SPcurrents,
            SPanim,
            numSP           // TODO remove
        };

        public static int scaleplot;        /* flag: plot gain scale? */

        public static POL polarization;     /* polarization for gain and currents displays */
        public enum POL    // note: if this order changes, find_fb() may also need to be changed
        {
            POLnone,
            POLhor,
            POLvert,
            POLlhcp,
            POLrhcp,
            POLcolour
        };

        public static int phasedirlock;     /* flag: is the direction used for calculation phases locked? */
        public static double phasephi, phasetheta;   /* if so, this is that direction; otherwise, it is identical to phi and theta */

        public static double phaseoffset;   /* phase offset for currents display */
        public static int distcor;          /* flag: correct phase for distance to observer? */
        public static int maxthickness;     /* line thickness used for maximum current in currents display */

        public static double animphase;     /* phase (in rad) of the "animation" of currents and E and H fields */
        public static double animfreq;      /* frequency of the animation: number of displayed cycles per second */
        public static double anim_update_interval;  /* time (in milliseconds) between updates of the animation display */

        public static double escale;        /* scale factor for the E field */
        public static double hscale;        /* scale factor for the H field */
        public static int show_p;           /* flag: display Poynting vector? */
        public static double qscale;        /* scale factor for the currents */
        public static double iscale;        /* scale factor for the charges */


        public static int win2sizex, win2sizey; /* size of freqplot window in pixels */
        public static int plot2_swr, plot2_z, plot2_z2, plot2_maxgain, plot2_dir, plot2_vgain; /* flags for the each of the frequency plots in window 2 */


        /* structure to contain, for one direction, the radiated power and polarization information */
        public class Gain
        {
            float total;   /* total radiated power, dBi */
            float axial;   /* axial ratio, with sign indicating the polarity: + left-hand circular, - right-hand */
            float tilt;    /* tilt angle of major axis, radians; defined as in NEC output: 0 = vertical, >0 topleft/bottomright */
        };

        public class Radpattern
        {    /* structure to contain (far field) radiation pattern for one frequency */
             //    double* gtheta;  /* pointer to array of theta values */
             //    double* gphi;    /* pointer to array of phi values */
             //    Gain** gain;     /* pointer to an array of pointers (one for every phi) to an array of Gain structs (one for every theta at that phi) */
             //    Point** gpo;     /* pointer to an array of pointers (one for every phi) to an array of Points (one for every theta at that phi) */
            int numtheta;    /* number of theta values */
            int numphi;      /* number of phi values */
                             //    struct RAdpattern *next;    /* pointer to the next such structure, used if more than one set of radiation pattern data is available for one frequency */
                             //   double* sin_gtheta;   /* pointer to an array containing the sines of the theta values */
                             //    double* sin_gphi;     /* ... etc. for phi */
                             //    double* cos_gtheta;   /* ... and cosines */
                             //    double* cos_gphi;
        };

        #region Segcur
        /// <summary>
        /// structure to contain the current (phase, magnitude, location, direction) information of one segment
        /// </summary>
        /// <remarks>
        /// note: part of the below information is in principle also available in the
        /// arrays describing the antenna's structure as read from the NEC input file;
        /// however, it seems more consistent to read this structure information now from the
        /// NEC output file, so we're sure we have the right data even if the input and output
        /// files are (accidentally) from different models.
        ///
        /// note2: the current in surface patches is also represented by the below structure;
        /// basically, a surface patch current is transformed into two equivalent segments,
        /// although in many practical cases one of them carries almost no current and can
        /// therefore be omitted.
        ///
        /// note3: for the animations, it is more convenient to have real and imaginary
        /// vectors. For segment currents, these obviously point in the same (or opposite)
        /// directions, namely along the segment. For surface currents, they are non-
        /// trivial; however, we don't need a second segment here (in contrast to the
        /// other representation), so only the first 'numanimseg' (see below, Currents
        /// structure) segments contain valid data in re[] and im[].
        ///
        /// note4: although the name suggests otherwise, this structure can also contain
        /// charge information.
        /// </remarks>
        public class Segcur
        {
            Point p0;       /* coordinate of segment endpoint */
            Point p1;       /* coordinate of segment endpoint */
            Point c;        /* coordinate of segment center point */
            float a;        /* amplitude of the current */
            float phase;    /* phase of the current (degrees) */
            float dist;     /* distance of center of segment from the plane through the origin, perpendicular to the viewing direction; positive if in front of this plane, negative otherwise */
            float area;     /* area of a surface patch */
            float[] re;     /* real vector */
            float[] im;     /* imaginary vector */
            float[] q0;     /* two unity vectors, orthogonal to each other and to the direction of the segment; used for drawing charges */
            float[] q1;     /* two unity vectors, orthogonal to each other and to the direction of the segment; used for drawing charges */
            float qre;      /* real value of the charge */
            float qim;      /* imaginary value of the charge */

            public Segcur()
                {
                this.p0 = new Nec.Point();
                this.p1 = new Nec.Point();
                this.c = new Nec.Point();
                this.re = new float[3];
                this.im = new float[3];
                this.q0 = new float[3];
                this.q1 = new float[3];
            }
        };
        #endregion
        #region Currents
        /// <summary>
        /// class to contain the current distribution for all segments at one frequency
        /// </summary>
        public class Currents
        {
            //   public Segcur* s;
            public int maxseg;       /* current size of s[] */
            public int numseg;       /* number of valid entries in s[] */
            public int numanimseg;   /* number of those that are needed for the animations; i.e., excluding the extra segments caused by the decomposition of the surface patch currents */
            public int numrealseg;   /* number of those that are actually segments; so s[numrealseg...numseg-1] are surface patches */
            public double maxI;      /* highest segment current */
            public double maxQ;      /* highest segment charge */
        };
        #endregion
        #region Nearfield
        /// <summary>
        /// structure to contain the electric and magnetic field at some point (near the antenna)
        /// </summary>
        public class NF
        {
            Point p;            /* location of the point */
            bool evalid;         /* flag to indicate presence of E data */
            bool hvalid;         /* flag to indicate presence of H data */
            float[] ere;        /* real component of E vector */
            float[] eim;        /* imaginary component of E vector */
            float[] hre;        /* real component of H vector */
            float[] him;        /* imaginary component of H vector */
            float[] poy;        /* time-averaged Poynting vector */
            //struct NF *next;   /* pointer to next datapoint */

            public NF()
            {
                this.p = new Nec.Point();
                this.evalid = false;
                this.hvalid = false;
                this.ere = new float[3];
                this.eim = new float[3];
                this.hre = new float[3];
                this.him = new float[3];
                this.poy = new float[3];
            }
        };

        public static List<NF> Nearfield = new List<NF>();
        #endregion

        #region Nec Output
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

        /// <summary>
        /// class to contain all output data that NEC produces for one frequency
        /// </summary>
        public class NECoutput
        {
            double f;          /* frequency in MHz */
            //Radpattern* rp;    /* pointer to radiation pattern data, if available; NULL otherwise */
            //public fixed double d[(int)NECO.neco_numdat];
            int rpgpovalid;    /* flag: !=0 if the gpo field of *rp already contains valid data for the present settings */
            //Currents* cu;      /* pointer to currents data, if available; NULL otherwise */
            //Nearfield* nf;     /* pointer to nearfield data, if available; NULL otherwise */
            double maxe, maxh;  /* largest values of E and H in nearfield data */
        };

        //extern NECoutput * neco;
        //extern int numneco;
        //extern int maxfreqs;

        //extern double globalmaxdb; /* maximum gain, global over all frequencies and all output data files */

        //extern int rp_index;       /* index of the entry in neco[] whose radiation pattern is being shown in the 3D plot */
        #endregion

        //extern double r0;          /* reference impedance for SWR calculations */

        //extern int window1open;
        //extern int window2open;

        //extern char* inputfilename;

        //extern int fontheight;

        //#ifndef GdkColor
        //   #include <gdk/gdk.h>
        //#endif
        //#define c_bg    col[0]
        //#define c_axis  col[1]
        //#define c_wire  col[2]
        //#define c_surf  col[3]
        //#define c_back  col[4]
        //#define c_gain  col[5]
        //#define c_scale col[6]
        //#define c_exci  col[7]
        //#define c_load  col[8]
        //#define c_netw  col[9]
        //#define c_inactive col[10]
        //#define c_efield   col[11]
        //#define c_hfield   col[12]
        //#define c_poynting col[13]
        //#define c_qpos     col[14]
        //#define c_qneg     col[15]
        //#define c_gain_lin   col[16]
        //#define c_gain_lhcp  col[17]
        //#define c_gain_rhcp  col[18]
        //#define c_currents (col+19)
        //#define NC NC_phase+19
        //extern GdkColor col[NC];

        //typedef struct {
        //   void (* SetLineAttributes) (unsigned int,int,int,int);
        //   void (* DrawLine) (double, double, double, double);
        //   void (* SetForeground) (GdkColor*);
        //   void (* ClearWindow) ();
        //   void (* DrawString) (double, double, char*, double, double);
        //   void (* Complete) ();
        //   void (* SetClipRectangle) (double, double, double, double);
        //   void (* ClearRectangle) (double, double, double, double);
        //} Outdev;
        //extern Outdev *out;   /* pointer to Outdev struct, determines destination of drawing commands issued through the below macros: */

        //#define SetLineAttributes(a,b,c,d) out->SetLineAttributes((a),(b),(c),(d))
        //#define DrawLine(a,b,c,d) out->DrawLine((a),(b),(c),(d))
        //#define SetForeground(a) out->SetForeground(a)
        //#define ClearWindow() out->ClearWindow()
        //#define DrawStringLL(a,b,c) DrawString(a,b,c,0,0)     /* draw text specifying lower left corner */
        //#define DrawStringUL(a,b,c) DrawString(a,b,c,0,1)     /* draw text specifying upper left corner */
        //#define DrawString(a,b,c,d,e) out->DrawString(a,b,c,d,e) /* draw text specifying a corner defined by d and e:
        //                                                            d = 0...1 for left...right; e = 0...1 for lower...upper*/
        //#define SetClipRectangle(a,b,c,d) out->SetClipRectangle(a,b,c,d)  /* set clip rectangle to (a,b)..(c,d) */
        //#define ClearRectangle(a,b,c,d) out->ClearRectangle(a,b,c,d)  /* clear part of window */
        //void process_nec_output(NECoutput* ne);     /* transform the gain data array into an array of points in 3D space */
        //void calc_vgain(void);                      /* update the vgain records in accordance with the viewing direction */
        //void mark_gpo_invalid(void);                /* mark gpo fields as invalid; must be called after every change of gain scale */
        //void calcswr(NECoutput* ne);                /* calculate the SWR from zr and zi */
        //void calcphiabs(NECoutput* ne);             /* calculate phi and abs of impedance from zr and zi */
        //void initX(int* argcp, char** argv);         /* initialize the X-window stuff */
        //void mainX(int really);                     /* main loop (wait for event, process event, redisplay) */
        //int Pending(void);                          /* test whether events are pending; returns !=0 if present (drawing) action should be interrupted */
        //void reread(void);                          /* reread the input file(s) */
        //void process_optionstring(char* s);         /* process a set of options, contained in the string s */
        //void draw_all(int ie);                      /* draw antenna structure and/or gain pattern; ie=interupt_enable */
        //void calcproj(void);                        /* calculate the projection matrix etc. */
        //int write_png(int which,const char* filename);    /* write window 1 or 2 to an PNG file; return !=0 in case of error */
        //int write_postscript(const char* filename, void (* drawfn)(int), int xsize, int ysize);
        ///* write postscript to file; return !=0 in case of error */
        //void draw_all2(int dummy);                  /* draw plot(s) vs. frequency */
        //double xfreq(int x);                        /* return frequency corresponding to X coordinate in frequency plot */
        //int freqx(double f);                        /* return X coordinate corresponding to frequency */
        //int freqindex(double f);                    /* return index in neco[] of data closest to frequency f */
        //void* mymalloc(size_t n);                   /* like malloc(), but exits with error message if memory not available; implementation in parse_output.c */
        //void* myrealloc(void*, size_t n);           /* like realloc(), but exits with error message if memory not available; implementation in parse_output.c */
        //void draw_opaque(Radpattern*);             /* draw radiation pattern as an opaque surface */
        //void setpolcolour(double axial);            /* set the colour corresponding to this axial ratio */
        //char* opaque_impossible(Radpattern* rp);    /* returns NULL if opaque drawing of this Radpattern is possible; otherwise, pointer to string with explanation */

        #region drem
        // if your system doesn't have the function drem() available, uncomment
        // the following line and recompile (after deleting *.o)
        //  #define drem(x,y) ( (x) - (y)*rint((x)/(y)) ) 
        // TODO remove this temporary function and use % where it's called
        public static double drem(double x, double y)
        {
            return x % y;
        }
        #endregion

        // --------------------- end of .h file -------------------------

        //        double extr = 0;
        //        double extr_str = 0;

        //        int gainscale = 0;
        //        int gainplot;
        //        int structplot;
        //        int scaleplot = 0;
        //        int phasedirlock = 0;
        //        double phasephi, phasetheta;
        //        int polarization = POLnone;
        //        int distcor = 1;
        //        double phaseoffset = 0;
        //        int maxthickness = 10;
        //        double animphase = 0;
        //        double animfreq = 1.0;
        //        double anim_update_interval = Default_anim_update_interval;
        //        double escale = 1;
        //        double hscale = 1;
        //        int show_p = 1;
        //        double qscale = 1;
        //        double iscale = 1;

        //        int g_argc;
        //        char** g_argv;

        //        char* inputfilename;

        //        int window1open = 0;
        //        int window2open = 0;

        //        char* expeps = NULL;
        //        char* exppng = NULL;

        //        /*---------- global variables, together containing the antenna structure -------------*/
        //        Wire* wires = NULL;
        //        Surface* surfaces;
        //        Exci* excis = NULL;
        //        Load* loads = NULL;
        //        Netw* netws = NULL;
        //        /*---------- global variables, together containing the gain distribution -------------*/

        //        int numneco = 0;
        //        int maxfreqs = 1;
        //        NECoutput* neco = NULL;

        //        int rp_index = -1;       /* index of the entry in neco[] whose radiation pattern is being shown in the 3D plot */

        //        /*------------------------------------------------------------------------------------*/

        public static void Init()
        {
            //    wires = mymalloc(maxwires * sizeof(Wire));
            //    surfaces = mymalloc(maxsurfaces * sizeof(Surface));
            //    excis = mymalloc(maxexcis * sizeof(Exci));
            //    loads = mymalloc(maxloads * sizeof(Load));
            //    netws = mymalloc(maxnetws * sizeof(Netw));
            //    neco = mymalloc(maxfreqs * sizeof(NECoutput));
        }

        public static void Info()
        {
            Console.WriteLine("Copyright (C) 1998-2003 Pieter-Tjerk de Boer -- pa3fwm@amsat.org");
            Console.WriteLine("Xnecview comes with ABSOLUTELY NO WARRANTY. This is free software, and you are");
            Console.WriteLine("welcome to redistribute it under the conditions of version 2 of the GNU GPL.");
            Console.WriteLine("For more information: see the files README and COPYING which should accompany");
            Console.WriteLine("this software, or ask the author.");
        }

        /// <summary>
        /// Reads and processes a NEC input or output file.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public static Boolean Readfile(string filename)
        {
            try
            {
                System.IO.StreamReader file = new System.IO.StreamReader(filename);
                ReadNecOutput rno = new ReadNecOutput();
                if (!rno.Read_nec_output(file)) // try interpreting as output data first
                {
                    file.DiscardBufferedData();
                    file.BaseStream.Seek(0, System.IO.SeekOrigin.Begin);
                    ReadNecInput rni = new ReadNecInput();
                    rni.Read_nec_input(file);   /* if data is not output data, try interpreting it as input data */
                }
            }
            catch
            {
                System.Console.WriteLine("Can't open \"{0}\" for reading", filename);
                //    usage();
                return false;
            }

            return true;
        }

        /// <summary>
        /// Clears all data storage and tries to re-read the file
        /// </summary>
        public static void Reread()
        {
            int i;
            int old_rp_index;

            extr = 0;
            //globalmaxdb = -1e30;
            //            while (numneco > 0)
            //            {  /* delete all old output data and free the associated Radpattern arrays */
            //                Radpattern* rp,*rp1;
            //                Nearfield* nf,*nf1;
            //                numneco--;
            //                rp = neco[numneco].rp;
            //                while (rp)
            //                {
            //                    for (i = 0; i < rp->numphi; i++)
            //                    {
            //                        free(rp->gain[i]);
            //                        free(rp->gpo[i]);
            //                    }
            //                    free(rp->gain);
            //                    free(rp->gpo);
            //                    free(rp->gphi);
            //                    free(rp->gtheta);
            //                    free(rp->sin_gphi);
            //                    free(rp->cos_gphi);
            //                    free(rp->sin_gtheta);
            //                    free(rp->cos_gtheta);
            //                    rp1 = rp;
            //                    rp = rp->next;
            //                    free(rp1);
            //                }
            //                if (neco[numneco].cu)
            //                {
            //                    free(neco[numneco].cu->s);
            //                    free(neco[numneco].cu);
            //                }
            //                nf = neco[numneco].nf;
            //                while (nf)
            //                {
            //                    nf1 = nf;
            //                    nf = nf->next;
            //                    free(nf1);
            //                }
            //            }
            //            old_rp_index = rp_index;
            //            rp_index = -1;

            //            for (i = 0; i < g_argc; i++) readfile(g_argv[i]);
            //            mark_gpo_invalid();
            //            if (rp_index >= 0)
            //            {
            //                if (old_rp_index >= 0 && old_rp_index < numneco
            //                      && (neco[old_rp_index].rp || neco[old_rp_index].cu || neco[old_rp_index].nf))
            //                    rp_index = old_rp_index;
            //                process_nec_output(neco + rp_index);
            //            }
            //            calcproj();
            //            calc_vgain();
        }
    }
}
