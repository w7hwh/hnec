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

            public Point()
            {
                x = y = z = 0;
            }

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
                rad = 0.0;
                itg = ns = abs = 0;
                dispnr = false;
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
                this.type = SU.SU_none;
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
            SU_none = 0,    /* not initialized yet */
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
            public int wireNo;   // pointer to Wire data
            public int seg;      // segment number within wire

            public Exci()
            {
                wireNo = 0;
                seg = 0;
            }

            public Exci(int newWire, int newSeg)
            {
                wireNo = newWire;
                seg = newSeg;
            }

            public override string ToString() => string.Format("Exci = {0}, {1}", seg, wireNo);
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
        /// <summary>
        /// class to describe a load
        /// </summary>
        /// <remarks>
        ///  note: if an LD card specifies a loading which extends over several wires,
        ///  it is represented here by one Load struct per wire
        /// </remarks>
        public class Load
        {
            public int wireNo;   // pointer to Wire data
            public int firstseg; // number of first loaded segment within wire
            public int lastseg;  // number of last loaded segment within wire

            /// <summary>
            /// 
            /// </summary>
            public Load()
            {
                wireNo = 0;
                firstseg = lastseg = 0;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public override string ToString() => string.Format("Load = {0}, {1}, {2}", wireNo, firstseg, lastseg);
        };

        /// <summary>
        /// 
        /// </summary>
        public static List<Load> loads = new List<Load>();

        /// <summary>
        /// 
        /// </summary>
        public static void DumpLoads()
        {
            Console.WriteLine("------- DumpLoads() dumping {0} loads -------", loads.Count);
            foreach (Load l in loads)
            {
                Console.WriteLine(l.ToString());
            }
            Console.WriteLine("----------------------------------------");
        }
        #endregion
        #region Netw
         /// <summary>
         /// class to describe a "network" in the antenna structure; 
         /// a transmission line is a special case of this
         /// </summary>
        public class Netw
        {
            public int wireNo1;       // pointer to Wire data of one end of the network / transmission line
            public int seg1;          // segment number within wire
            public int wireNo2;       // same for the other end of the network / transmission line
            public int seg2;          // segment number within wire
            public NETW_TYPE type;

            /// <summary>
            /// 
            /// </summary>
            public Netw()
            {
                wireNo1 = wireNo2 = 0;
                seg1 = seg1 = 0;
                type = NETW_TYPE.generic;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="w1"></param>
            /// <param name="s1"></param>
            /// <param name="w2"></param>
            /// <param name="s2"></param>
            /// <param name="t"></param>
            public Netw(int w1, int s1, int w2, int s2, NETW_TYPE t)
            {
                wireNo1 = w1;
                seg1 = s1;
                wireNo2 = w2;
                seg2 = s2;
                type = t;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public override string ToString() => string.Format("Netw = {0}, {1}, {2}, {3}, {4}", wireNo1, seg1, wireNo2, seg2, type);
        };

        public enum NETW_TYPE
        {
            generic = 0,            // generic network
            noPhaseReversal = 1,    // transmission line without phase reversal
            phaseReversal = -1      // transmission line with phase reversal
        }
        public static List<Netw> netws = new List<Netw>();

        /// <summary>
        /// 
        /// </summary>
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

            /// <summary>
            /// 
            /// </summary>
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

 
        //extern double r0;          /* reference impedance for SWR calculations */

        //extern int window1open;
        //extern int window2open;

        //extern char* inputfilename;

        //extern int fontheight;

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
