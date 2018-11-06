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
#undef HAVE_GETOPT
#undef HAVE_LIBPNG

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Necview
{
    public class Program
    {
        public static readonly string VERSION = "1.36";

        private static int Maxopts = 100;

        //        char* myargv[Maxopts];
        //        int myargc;
        //        double frequency = -1;

        private static void Info()
        {
            Console.WriteLine("Copyright (C) 1998-2003 Pieter-Tjerk de Boer -- pa3fwm@amsat.org");
            Console.WriteLine("Xnecview comes with ABSOLUTELY NO WARRANTY. This is free software, and you are");
            Console.WriteLine("welcome to redistribute it under the conditions of version 2 of the GNU GPL.");
            Console.WriteLine("For more information: see the files README and COPYING which should accompany");
            Console.WriteLine("this software, or ask the author.");
        }

        private static void Usage()
        {
            Console.Write("Usage: xnecview [options] filenames [options]\n"
                            + "filenames   : names of NEC2 input and/or output files to be displayed\n"
                            + "options:\n"
                            + "  -h or --help : show this information\n"
#if (HAVE_GETOPT)
                   + "  --struct     : set structure view to 'struct'\n"
                           + "  --tags       : set structure view to 'struct+tags'\n"
                           + "  --currents   : set structure view to 'currents'\n"
                            + "  --animation  : set structure view to 'animation'\n"
                            + "  --slice      : set radiation view to 'slice'\n"
                            + "  --frame      : set radiation view to 'frame'\n"
                            + "  --opaque     : set radiation view to 'opaque'\n"
                            + "  --near       : set radiation view to 'near field'\n"
                            + "  --linpower   : set radiation scale linear in power\n"
                            + "  --linvoltage : set radiation scale linear in voltage\n"
                            + "  --arrl       : set radiation scale to ARRL style\n"
                            + "  --log        : set radiation scale to logarithmic\n"
                            + "  --qscale num : set charges scale (animation)\n"
                            + "  --iscale num : set currents scale (animation)\n"
                            + "  --escale num : set electric field scale\n"
                            + "  --hscale num : set magnetic field scale\n"
                            + "  --hidepoynting : hide Poynting vector in near field display\n"
                            + "  --afreq num  : set animation frequency (Hz)\n"
                            + "  --aphase num : set animation phase (degrees)\n"
                            + "  --aupdate num : set time between animation updates (milliseconds, default 100)\n"
                            + "  --freq num   : set frequency (MHz)\n"
                            + "  --z0 num     : set reference impedance (ohm)\n"
                            + "  --expeps filename : only export picture to .eps-file\n"
#if (HAVE_LIBPNG)
                   + "  --exppng filename : only export picture to .png-file\n"
#endif
                    + "  --view phi,theta,zoom,trx,try : set viewing direction and zoom\n"
                            + "Note: typing 'v' prints the current values for all of these settings,\n"
                           + "for easy copying into scripts.\n"
#endif
                   );
            //     exit(1);
        }



        private void Process_options()
        {
#if (HAVE_GETOPT)
            //   struct option longopts[]={
            //      { "near",     0, &gainplot, GPnearfield
            //    },
            //      { "slice",    0, &gainplot, GPslice
            //},
            //      { "frame",    0, &gainplot, GPframe },
            //      { "opaque",   0, &gainplot, GPopaque },

            //      { "struct",    0, &structplot, SPstruct },
            //      { "tags",      0, &structplot, SPtags },
            //      { "currents",  0, &structplot, SPcurrents },
            //      { "animation", 0, &structplot, SPanim },

            //      { "linpower",  0, &gainscale, GSlinpower },
            //      { "linvoltage",  0, &gainscale, GSlinvolt },
            //      { "arrl", 0, &gainscale, GSarrl },
            //      { "log",  0, &gainscale, GSlog },

            //      { "pol", 1, NULL, 40 },

            //      { "iscale", 1, NULL, 13 },
            //      { "qscale", 1, NULL, 14 },
            //      { "escale", 1, NULL, 15 },
            //      { "hscale", 1, NULL, 16 },
            //      { "hidepoynting", 0, &show_p, 0 },

            //      { "view", 1, NULL, 17 },
            //      { "freq", 1, NULL, 18 },
            //      { "afreq", 1, NULL, 19 },
            //      { "aphase", 1, NULL, 20 },
            //      { "aupdate", 1, NULL, 21 },

            //      { "expeps", 1, NULL, 30 },
#if (HAVE_LIBPNG)
            //      { "exppng", 1, NULL, 31 },
#endif
            //      { "z0", 1, NULL, 32 },

            //      { "help", 0, NULL, 'h' },

            //      { NULL, 0, NULL, 0 }
            //   };

            //   while (1) {
            //      int c;
            //c = getopt_long(myargc, myargv, "-h", longopts, NULL);
            //      if (c==-1) break;
            //      switch (c) {
            //         case 1: 
            //            g_argv[g_argc++]=optarg;
            //            readfile(optarg);
            //            if (!inputfilename) inputfilename=optarg;
            //            break;
            //         case 11: phi=atof(optarg); break;
            //         case 12: theta=atof(optarg); break;
            //         case 13: iscale=atof(optarg); break;
            //         case 14: qscale=atof(optarg); break;
            //         case 15: escale=atof(optarg); break;
            //         case 16: hscale=atof(optarg); break;
            //         case 17: sscanf(optarg,"%lg,%lg,%lg,%lg,%lg",&phi,&theta,&zoom,&trx,&try); break;
            //         case 18: frequency=atof(optarg); break;
            //         case 19: animfreq=atof(optarg); break;
            //         case 20: animphase=atof(optarg); animphase*=M_PI/180; break;
            //         case 21: anim_update_interval=atof(optarg); break;
            //         case 30: expeps=optarg; break;
#if (HAVE_LIBPNG)
            //         case 31: exppng=optarg; break;
#endif
            //         case 32: r0=atof(optarg); break;
            //         case 40: if (!strcmp(optarg,"total")) polarization=POLnone;
            //                  else if (!strcmp(optarg,"hor")) polarization=POLhor;
            //                  else if (!strcmp(optarg,"vert")) polarization=POLvert;
            //                  else if (!strcmp(optarg,"lhcp")) polarization=POLlhcp;
            //                  else if (!strcmp(optarg,"left")) polarization=POLlhcp;
            //                  else if (!strcmp(optarg,"rhcp")) polarization=POLrhcp;
            //                  else if (!strcmp(optarg,"right")) polarization=POLrhcp;
            //                  else if (!strcmp(optarg,"colour")) polarization=POLcolour;
            //                  else if (!strcmp(optarg,"color")) polarization=POLcolour;
            //                  break;

            //         case 'h': info(); usage(); break;
            //         case '?': printf("Try `%s -h' for help.\n", myargv[0]); break;
            //      }
            //   }
            //   while (optind!=myargc) {
            //      g_argv[g_argc++]=myargv[optind];
            //      readfile(myargv[optind]);
            //      if (!inputfilename) inputfilename=myargv[optind];
            //      optind++;
            //   }
#else  // fallback in case we don't have getopt.h:
            //   int i;
            //   if (myargc==1) { info(); usage(); return; }
            //   for (i=1;i<myargc;i++) {
            //      if (myargv[i][0]=='-') { info(); usage(); break; }
            //      g_argv[g_argc++]=myargv[i];
            //      readfile(myargv[i]);
            //      if (!inputfilename) inputfilename=myargv[i];
            //   }
#endif
        }

        public static void Process_optionstring(string s)
        {
#if (HAVE_GETOPT)
            //    /* we would like to just call getopt_long() to process these options,
            //       but that doesn't work: process_optionstring() is called from
            //       inside read_nec_input(), which on its turn is called while processing
            //       (with getopt_long()) the real command line. It seems getopt_long()
            //       uses quite a few hidden state variables, making such a recursive use
            //       impossible. Therefore, we insert the new options into the original
            //       command line option list.
            //    */
            //    int ac;
            //    char** av;
            //    char* ar;
            //    ar = mymalloc(strlen(s) + 1);
            //    strcpy(ar, s);
            //    s = ar;
            //    av = mymalloc((1 + strlen(s)) * sizeof(char*));
            //    ac = 0;
            //    s = strtok(s, " \n\r");
            //    while (s)
            //    {
            //        av[ac++] = s;
            //        s = strtok(NULL, " \n\r");
            //    }
            //    if (ac + myargc <= Maxopts)
            //    {
            //        memmove(myargv + optind + ac, myargv + optind, (myargc - optind) * sizeof(char*));
            //        memcpy(myargv + optind, av, ac * sizeof(char*));
            //        myargc += ac;
            //    }
            //    free(av);
#endif
        }

        private static void Exit()
        {
            if (System.Windows.Forms.Application.MessageLoop)
            {
                // WinForms app
                System.Windows.Forms.Application.Exit();
            }
            else
            {
                // Console app
                System.Environment.Exit(1);
            }
        }

        public static void Main(params string[] args)
        {
            int i;
            int argc = args.Length;

            Console.WriteLine("XNECVIEW {0}", VERSION);

            if (argc == 1)
            {
                Info();
                Usage();
                Exit();   
            }
            //    initX(&argc, argv);
            //    inputfilename = NULL;

            //    setlocale(LC_NUMERIC, "C");

            Nec.Init();     // initialize the empty Lists and variables

            //    structplot = -1;
            //    gainplot = -1;
            //    g_argc = 0;
            //    g_argv = mymalloc(argc * sizeof(char*));
            //    for (i = 0; i < argc; i++) myargv[i] = argv[i];
            //    myargc = argc;
            //    process_options();

            //    if (numneco == 0 && numwires == 0 && numsurfaces == 0)
            //    {
            //        fprintf(stderr, "No data!\n");
            //        return 1;
            //    }

            //    if (rp_index >= 0)
            //    {
            //        if (frequency >= 0) rp_index = freqindex(frequency);
            //        process_nec_output(neco + rp_index);
            //        if (gainplot < 0)
            //        {
            //            if (neco[rp_index].rp == NULL && neco[rp_index].nf != NULL) gainplot = GPnearfield;
            //            else
            //            {
            //                for (i = 0; i < numneco; i++)
            //                    if (neco[i].rp || neco[i].nf)
            //                    {
            //                        gainplot = 1111;
            //                        break;
            //                    }
            //                if (gainplot < 0) gainplot = GPnone;
            //            }
            //        }
            //    }
            //    else
            //    {
            //        gainplot = GPnone;
            //    }

            //    if (numwires + numsurfaces > 0)
            //    {
            //        window1open = 1;
            //        if (structplot < 0)
            //        {
            //            if (rp_index >= 0) structplot = SPstruct;
            //            else structplot = SPtags;
            //        }
            //    }
            //    else
            //    {
            //        for (i = 0; i < numneco; i++)
            //            if (neco[i].cu)
            //            {
            //                if (structplot < 0 && gainplot < 0)
            //                {
            //                    structplot = SPcurrents;
            //                    window1open = 1;
            //                }
            //                break;
            //            }
            //        if (rp_index >= 0) window1open = 1;
            //        if (structplot < 0) structplot = SPnone;
            //    }

            //    if (numneco > 1)
            //    {
            //        window2open = 1;
            //        if (neco[rp_index].rp == NULL) { plot2_maxgain = 0; plot2_z = 1; }
            //    }

            //    if (numneco == 1 && !window1open)
            //    {
            //        fprintf(stderr, "No data available which can be displayed graphically.\n");
            //        return 0;
            //    }

            //    if (gainplot == GPopaque || gainplot == 1111)
            //    {
            //        char* s;
            //        Radpattern* rp;
            //        rp = neco[rp_index].rp;
            //        s = opaque_impossible(rp);
            //        if (s)
            //        {
            //            fprintf(stderr, "Opaque gain plot not possible; reason: %s.\n", s);
            //            gainplot = GPframe;
            //        }
            //        else if (gainplot == 1111)
            //        {
            //            if (rp->next)
            //            {
            //                gainplot = GPframe;
            //                fprintf(stderr, "Opaque gain plot off by default, because of multiple RP cards.\n");
            //            }
            //            else
            //            {
            //                gainplot = GPopaque;
            //            }
            //        }
            //    }

            //    calcproj();
            //    calc_vgain();

            //    mainX(!expeps && !exppng);

            //    if (expeps)
            //    {
            //        int err;
            //        if (!window1open && window2open) err = write_postscript(expeps, Postscript.PSDRAW.TWO, win2sizex, win2sizey);
            //        else err = write_postscript(expeps, Postscript.PSDRAW.TWO, winsizex, winsizey);
            //        if (err)
            //        {
            //            fprintf(stderr, "Error writing postscript\n");
            //            return 1;
            //        }
            //    }
# if (HAVE_LIBPNG)
            //    if (exppng)
            //    {
            //        int err;
            //        if (!window1open && window2open)
            //        {
            //            err = write_png(2, exppng);
            //        }
            //        else
            //        {
            //            err = write_png(1, exppng);
            //        }
            //        if (err)
            //        {
            //            fprintf(stderr, "Error writing PNG file\n");
            //            return 1;
            //        }
            //    }
#endif

            //    return 0;
        }

    }
}
