/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 1998-2006, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 *
 *  This module contains all low level postscript drawing routines.
 */
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    class Postscript : IOutdevInterface
    {
        public enum PSDRAW
        {
            ONE,
            TWO
        };
        //        FILE* psfile = NULL;

        int needstroke;
        int needmoveto;
        Color xcprev;

        public Postscript()
        {
        }

        public void DrawLine(double x1, double y1, double x2, double y2)
        {
            double lastx = -1, lasty = -1;
            //            if (x1 != lastx || y1 != lasty || needmoveto)
            //            {
            //                if (needstroke) fprintf(psfile, "stroke\n");
            //                fprintf(psfile, "%g %g m\n", x1, y1);
            //                needmoveto = 0;
            //            }
            //            fprintf(psfile, "%g %g l\n", x2, y2);
            //            lastx = x2; lasty = y2;
            //            needstroke = 1;
        }

        public void SetLineAttributes(uint line_width, int line_style, int cap_style, int join_style)
        {
            //            if (needstroke) { fprintf(psfile, "stroke\n"); needstroke = 0; }
            //            if (line_width == 0) line_width++;   /* Zero width lines for X11 speed of drawing but in
            //					 postscript don't really want the thinnest lines
            //					 which are possible on the device! */
            //            fprintf(psfile, "%u setlinewidth ", line_width);
            //            switch (line_style)
            //            {
            //                case GDK_LINE_SOLID:
            //                    fprintf(psfile, "[] 0 setdash\n");
            //                    break;
            //                case GDK_LINE_ON_OFF_DASH:
            //                    fprintf(psfile, "[5 5] 0 setdash\n");
            //                    break;
            //            }
            //            if (cap_style != GDK_CAP_ROUND || join_style != GDK_JOIN_ROUND || line_style == GDK_LINE_DOUBLE_DASH)
            //            {
            //                puts("Error in ps_SetLineAttributes!");
            //            }
            //            needmoveto = 1;
        }

        public void SetForeground(Color color)
        {
            //if (xc == xcprev) return;
            //xcprev = xc;
            //            if (needstroke) { fprintf(psfile, "stroke\n"); needstroke = 0; }
            //            fprintf(psfile, "%g %g %g setrgbcolor\n", xc->red / 65535., xc->green / 65535., xc->blue / 65535.);
            //            needmoveto = 1;
        }

        public void ClearWindow()
        {
        }

        public void DrawString(double a, double b, string s, double d, double e)
        {
            //            fprintf(psfile, "gsave (%s) dup stringwidth pop %g mul %g add\n", s, -d, a);
            //            fprintf(psfile, "%g descent sub ascent descent add %g mul add m 1 -1 scale show grestore\n", b, e);
        }

        public void Complete()
        {
        }

        public void SetClipRectangle(double x1, double y1, double x2, double y2)
        {
            //            if (needstroke) { fprintf(psfile, "stroke\n"); needstroke = 0; }
            //            xcprev = NULL;
            //            fprintf(psfile, "grestore gsave\n");
            //            fprintf(psfile, "newpath %g %g moveto %g %g lineto %g %g lineto %g %g lineto closepath clip newpath\n",
            //                                    x1, y1, x1, y2, x2, y2, x2, y1);
        }

        public void ClearRectangle(double a, double b, double c, double d)
        {
        }

        private static int PSsize = 400;
        private static int PSllx = 50;
        private static int PSlly = 150;

        int Write_postscript(string filename, PSDRAW drawfn, int xsize, int ysize)
        {
            //   Outdev* oldout;
            //        time_t ti;

            //   psfile=fopen(filename,"w");
            //   if (!psfile) return 1;
            //   oldout=out;
            //   out=&out_ps;

            //   time(&ti);
            //        fprintf(psfile,"%%!PS-Adobe-2.0 EPSF-2.0\n"
            //                  "%%%%DocumentFonts: %s\n"
            //                  "%%%%Title: Xnecview output\n"
            //                  "%%%%Creator: Xnecview %s\n"
            //                  "%%%%CreationDate: %s"
            //                  "%%%%Pages: 1\n"
            //                  "%%%%BoundingBox: %g %g %g %g\n"
            //                  "%%%%EndComments\n",
            //                       PSFONT,
            //                       VERSION,
            //                       ctime(&ti),
            //                  (double) PSllx,(double)PSlly,(double) PSllx+xsize,(double) PSlly+ysize);

            //   fprintf(psfile,"/xnecview 10 dict def\n xnecview begin \n");   /* use a private dictionary as recommended in EPSF-3.0 standard */
            //        fprintf(psfile,"/l /lineto load def\n");        /* define abbreviations l and m for lineto and moveto */
            //        fprintf(psfile,"/m /moveto load def\n");
            //        fprintf(psfile,"end\n");                  /* end use (definition) of local dictionary */
            //        fprintf(psfile,"%%%%EndProlog\n");
            //        fprintf(psfile,"%%%%Page: 1 1\n");
            //        fprintf(psfile,"xnecview begin\n");                 /* use local dictionary again */
            //        fprintf(psfile,"%g %g translate 1 -1 scale\n",     /* set scaling and translation such that our Xwindow coordinates can be used */
            //           (double) PSllx,
            //           (double) PSlly+ysize );
            //        fprintf(psfile,"newpath 0 0 moveto 0 %g lineto %g %g lineto %g 0 lineto closepath clip newpath\n",
            //                       (double) ysize, (double) xsize, (double) ysize, (double) xsize);
            //        fprintf(psfile,"1 setlinejoin 1 setlinecap\n");
            //        fprintf(psfile,"/%s findfont %i scalefont setfont\n", PSFONT, fontheight);    /* prepare the font */
            //        fprintf(psfile,"gsave 0 0 moveto (10ALIZgjy) true charpath pathbbox grestore /ascent exch def pop -1 mul /descent exch def pop\n");  /* estimate ascent and descent of font */
            //        fprintf(psfile,"gsave\n");  /* this gsave will be grestore'd by the first SetClipRectangle call */

            //        needstroke=0;
            //   needmoveto=0;
            //   xcprev=NULL;

            switch (drawfn)
            {
                case PSDRAW.ONE:
                    Draw.Draw_all(0);
                    break;

                case PSDRAW.TWO:
                    //FreqPlot.Draw_all2(0);
                    break;
            }

            //#if 0
            //   fprintf(psfile,"save 0 setgray 10 %i moveto 1 -1 scale (Antenna: %s) show restore\n",
            //      fontheight,inputfilename);
            //   if (gainplot!=GPnone) {
            //      fprintf(psfile,"save 0 setgray 10 %i moveto 1 -1 scale (Scale: %s%s%s) show restore\n",
            //         2*fontheight,
            //         GSnames[gainscale],
            //         scaleplot ? ", " : "",
            //         scaleplot ? GSnumbers[gainscale] : "");
            //   }
            //#endif

            //   if (needstroke) fprintf(psfile,"stroke\n");
            //        fprintf(psfile,"grestore\n");
            //        fprintf(psfile,"end\n");         /* end use of local dictionary */
            //        fprintf(psfile,"showpage\n");

            //   out=oldout;
            //   return ferror(psfile) | fclose(psfile);
            return 0;   //TODO fix so it returns real status
        }
    }
}
