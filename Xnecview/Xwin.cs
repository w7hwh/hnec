/*  XNECVIEW - a program for visualizing NEC2 input and output data
*
*  Copyright (C) 1998-2006,2011, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
*
*  Distributed on the conditions of version 2 of the GPL: see the files
*  README and COPYING, which accompany this source file.
*
*  This module contains all X-windows related stuff; this includes all
*  event handlers, and thus most of the interaction with the user.
*/
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    class Xwin : IOutdevInterface
    {
        //# include "icon.xbm"

        //        /* note: gdraw2, ggc2, etc. refer to the second window, i.e., the window containing the plots of SWR etc. vs. frequency */
        //        GtkWidget* gdraw1 = NULL,*gdraw2=NULL;
        //GtkWidget* vbox1;    /* the (GTK) vertical box for window 1, containing the top row of buttons, the drawing itself, and the bottom row of buttons */
        //        GtkWidget* botrow_curr;   /* the (GTK) bottom row of buttons in window 1 for settings of the currents display */
        //        GtkWidget* botrow_anim;   /* the (GTK) bottom row of buttons in window 1 for settings of the (animated) near field display */
        //        GdkGC* ggc,*ggc2;
        //GdkFont* gfont;

        //        GdkPixmap* gbackg = NULL;       /* picture buffer used for double buffering code; it will serve as the background pixmap for the window */
        //        GtkLabel* msgwidget;     /* label widget in top right corner, used for several messages */
        //        int fontheight;
        //        GdkPixmap* gbackg2 = NULL;      /* for window2, we draw into this pixmap */
        //        GdkColormap* gcm;
        //        int depth;

        //        int redraw = 1;        /* flag which signifies need for redrawing of struct/gain plot */
        //        int dragging = 0;      /* flag to indicate that user is dragging the struct/gain plot */
        //        int quick;           /* flag to indicate that drawing should be rough but quick */
        //        int quit = 0;          /* flag to indicate that the user wants to quit */
        //        int redraw2 = 1;       /* flag which signifies need for redrawing of frequency plots; extra for speedup: if it ==2, then only the vgain plot needs to be redrawn */
        //        int animbusy = 0;      /* flag to indicate that we're still busy processing the previous animation picture, to slow down the animation sufficiently on slow machines and/or complicated models */

        //        double old_animfreq = 0;   /* storage for animation frequency when animation is temporarily suspended using the 'z' key */

        //        Outdev*out=NULL;


        ///* color values for several items: */
        //GdkColor col[NC];
        //        GdkColor gcol[NC];


        //        /* -------------------- X windows stuff: drawing ---------------------------------------- */

        public void SetLineAttributes(uint a, int b, int c, int d)
        {
            //gdk_gc_set_line_attributes(ggc, a, b, c, d);
        }

        public void DrawLine(double a, double b, double c, double d)
        {
            //#if 0
            //   if (a>32000 || a<-32000 
            //      || b>32000 || b<-32000 
            //      || c>32000 || c<-32000 
            //      || d>32000 || d<-32000 
            //      )
            //#endif
            //            if (a > winsizex || a < 0
            //               || b > winsizey || b < 0
            //               || c > winsizex || c < 0
            //               || d > winsizey || d < 0
            //               )
            //            {
            //                /* clipping is needed, otherwise we run into problems when feeding these extremely large coordinates to the X server; apparently, the X11 protocol uses only 16 bits for the coordinates */
            //                double h;

            //                if (a > winsizex && c > winsizex) return;
            //                if (a < 0 && c < 0) return;
            //                if (b > winsizey && d > winsizey) return;
            //                if (b < 0 && d < 0) return;

            //                if (c < a) { h = a; a = c; c = h; h = b; b = d; d = h; }
            //                if (a < 0)
            //                {
            //                    b -= (d - b) * a / (c - a);
            //                    a = 0;
            //                }
            //                if (c > winsizex)
            //                {
            //                    d -= (d - b) * (c - winsizex) / (c - a);
            //                    c = winsizex;
            //                }

            //                if (d < b) { h = a; a = c; c = h; h = b; b = d; d = h; }
            //                if (b < 0)
            //                {
            //                    a -= (c - a) * b / (d - b);
            //                    b = 0;
            //                }
            //                if (d > winsizey)
            //                {
            //                    c -= (c - a) * (d - winsizey) / (d - b);
            //                    d = winsizey;
            //                }
            //            }
            //            gdk_draw_line(gbackg, ggc, (int)(a + 0.5), (int)(b + 0.5), (int)(c + 0.5), (int)(d + 0.5));
        }

        public void SetForeground(Color xc)
        {
            //gdk_gc_set_foreground(ggc, xc);
        }

        public void ClearWindow()
        {
            //X_SetForeground(&c_bg);
            //gdk_draw_rectangle(gbackg, ggc, TRUE, 0, 0, winsizex, winsizey);
        }


        public void DrawString(double a, double b, string s, double d, double e)    /* draw string */
        {
            //            if (d > 0) a -= d * gdk_string_width(gfont, s);
            //            b += -gfont->descent + e * (gfont->descent + gfont->ascent);
            //            gdk_draw_string(gbackg, gfont, ggc, (int)(a + 0.5), (int)(b + 0.5), s);
        }


        public void Complete()
        {
            //            gdk_window_set_back_pixmap(gdraw1->window, gbackg, 0);
            //            gdk_window_clear(gdraw1->window);
        }

        public void SetClipRectangle(double x1, double y1, double x2, double y2)
        {
            //            GdkRectangle gr;
            //            gr.x = (int)(x1 + 0.5);
            //            gr.y = (int)(y1 + 0.5);
            //            gr.width = (int)(x2 - gr.x + 1.5);
            //            gr.height = (int)(y2 - gr.y + 1.5);
            //            gdk_gc_set_clip_rectangle(ggc2, &gr);
        }

        public void ClearRectangle(double x1, double y1, double x2, double y2)
        {
            //            X2_SetClipRectangle(0, 0, win2sizex, win2sizey);
            //            X2_SetForeground(&c_bg);
            //            gdk_draw_rectangle(gbackg2, ggc2, TRUE, (int)(x1 + 0.5), (int)(y1 + 0.5), (int)(x2 - x1 + 1.5), (int)(y2 - y1 + 1.5));
        }

        //        Outdev outX = { 
        //              X_SetLineAttributes, 
        //              X_DrawLine, 
        //              X_SetForeground, 
        //              X_ClearWindow, 
        //              X_DrawString, 
        //              X_Complete, 
        //              NULL };


        //        /* -------------------- support function for querying the image ---------------------------- */


        //        int query_pixmap(GdkPixmap* w, int x, int y)
        //        /* returns index in phase colours array of colour at (x,y); or -1 if the colour is c_inactive, or -2 if the colour is yet another one */
        //        {
        //            int i;
        //            guint32 u;
        //            GdkImage* im;
        //            int width, height;

        //            gdk_window_get_size(w, &width, &height);
        //            if (x < 0 || y < 0 || x >= width || y >= height) return -2;

        //            im = gdk_image_get(w, x, y, 1, 1);
        //            u = gdk_image_get_pixel(im, 0, 0);
        //            gdk_image_destroy(im);
        //            if (u == c_bg.pixel) return -2;  /* most common case first :-) */
        //            for (i = 0; i < NC_phase; i++)
        //                if (c_currents[i].pixel == u) return i;
        //            if (u == c_inactive.pixel) return -1;
        //            return -2;
        //        }



        //        /* -------------------- export a widget as a PNG file ---------------------------------------- */

        //# ifdef HAVE_LIBPNG
        //        int write_png(int which,const char* filename)
        //{
        //   GdkImage* image;
        //        GdkWindow* w;
        //        FILE* f;
        //        int width, height;
        //        int x, y;
        //        png_structp pp;
        //        png_infop ip;
        //        unsigned char* buf;
        //        png_color pal[NC];
        //        int i;

        //   if (which==1) {
        //      if (!out) { out=&outX; draw_all(0);
        //    }
        //    w=gbackg;
        //   } else {
        //      if (!out) { out=&outX2; draw_all2(0); }
        //      w=gbackg2;
        //   }
        //   gdk_window_get_size(w,&width,&height);

        //f=fopen(filename,"wb");
        //   if (!f) return 1;

        //   image = gdk_image_get(w, 0, 0, width, height);
        //   if (!image) {
        //      fclose(f);
        //      return 1;
        //   }

        //   pp=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        //   if (!pp) {
        //      fclose(f);
        //      return 1;
        //   }
        //   ip=png_create_info_struct(pp);
        //   if (!ip) {
        //      png_destroy_write_struct(&pp, (png_infopp) NULL);
        //fclose(f);
        //      return 1;
        //   }
        //   if (setjmp(pp->jmpbuf)) {
        //      png_destroy_write_struct(&pp,&ip);
        //fclose(f);
        //gdk_image_destroy(image);
        //      return 1;
        //   }
        //   png_init_io(pp, f);
        //png_set_IHDR(pp, ip, width, height,8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
        //   for (i=0;i<NC;i++) {
        //      pal[i].red=  (col[i].red>>8);
        //      pal[i].green=(col[i].green>>8);
        //      pal[i].blue= (col[i].blue>>8);
        //   }
        //   png_set_PLTE(pp, ip, pal, NC);
        //png_write_info(pp, ip);

        //buf=mymalloc(width);
        //   for (y=0;y<height;y++) {
        //      for (x=0;x<width;x++) {
        //         guint32 u;
        //u = gdk_image_get_pixel(image, x, y);
        //         for(i=0;i<NC;i++) if (col[i].pixel==u) break;
        //         buf[x]=i;
        //      }
        //      png_write_row(pp, buf);
        //   }
        //   free(buf);

        //png_write_end(pp, NULL);
        //png_destroy_write_struct(&pp,&ip);

        //fclose(f);
        //gdk_image_destroy(image);
        //   return 0;
        //}
        //#endif




        ///* -------------------- event handlers (and their supporting code) for window 1 (struct/gain plot) ---------------------- */

        //char GSnumbers[][32]={" 0 -1 -3 -6 -10 dB",
        //                      " 0 -1 -3 -6 -10 dB",
        //                      " 0 -3 -10 -20 -30 dB",
        //                      " 0 -10 -20 -30 dB"};

        //void upd_msg(void)
        //{
        //    char s[32];
        //    return;
        //    if (scaleplot && (gainplot == GPslice || gainplot == GPframe || gainplot == GPopaque))
        //    {
        //        gtk_label_set_text(msgwidget, GSnumbers[gainscale]);
        //    }
        //    else
        //    {
        //        sprintf(s, "phi=%g  theta=%g", rint(phi), rint(theta));
        //        gtk_label_set_text(msgwidget, s);
        //    }
        //}


        //void setrpfreq(void)
        //{
        //    int x;
        //    redraw = 1;
        //    process_nec_output(neco + rp_index);
        //    if (window2open)
        //    {
        //        gdk_window_clear(gdraw2->window);
        //        gdk_gc_set_foreground(ggc2, &c_back);
        //        x = freqx(neco[rp_index].f);
        //        gdk_draw_line(gdraw2->window, ggc2, x, 0, x, win2sizey);
        //    }
        //}




        //gint resize_event(GtkWidget* w, GdkEventConfigure* ev, gpointer dummy)
        //{
        //    if (!gbackg) return TRUE;
        //    if (winsizex == ev->width && winsizey == ev->height) return TRUE;
        //    winsizex = ev->width;
        //    winsizey = ev->height;
        //    gdk_pixmap_unref(gbackg);
        //    gbackg = gdk_pixmap_new(w->window, winsizex, winsizey, gdk_drawable_get_depth(w->window));
        //    calcproj();
        //    redraw = 1;
        //#if 0
        //   w->requisition.width=winsizex;
        //   w->requisition.height=winsizey;
        //#endif
        //    return TRUE;
        //}


        //int lastx = 0;
        //int lasty = 0;
        //int origx = 0;
        //int origy = 0;

        //gint buttonpress_event(GtkWidget* w, GdkEventButton* ev, gpointer dummy)
        //{
        //    origx = lastx = ev->x;
        //    origy = lasty = ev->y;
        //    return TRUE;
        //}


        //gint buttonrelease_event(GtkWidget* w, GdkEventButton* ev, gpointer dummy)
        //{
        //    if (dragging)
        //    {
        //        dragging = 0;
        //        if (quick)
        //        {
        //            quick = 0;
        //            redraw = 1;
        //        }
        //        if (structplot == SPcurrents) redraw = 1;  /* this redraw is needed to grab a correct copy of the phaseplot */
        //        return TRUE;
        //    }

        //    if (ev->button == 1)
        //    {
        //        zoom *= 1.4142;
        //        trx -= (double)(ev->x - winsizex / 2) / zoom / winsize;
        //        try-= (double)(ev->y - winsizey / 2) / zoom / winsize;
        //        calcproj();
        //        redraw = 1;
        //        calc_vgain();
        //        if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //        }

        //   if (ev->button == 2)
        //        {
        //            zoom = ini_zoom;
        //            phi = ini_phi;
        //            theta = ini_theta;
        //            scaleplot = 0;
        //            trx = ini_trx;
        //            try= ini_try;
        //            upd_msg();
        //            calcproj();
        //            redraw = 1;
        //            calc_vgain();
        //            if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //            }

        //   if (ev->button == 3)
        //            {
        //                zoom /= 1.4142;
        //                calcproj();
        //                redraw = 1;
        //                calc_vgain();
        //                if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //            }
        //            return TRUE;
        //        }


        //        gint motion_event(GtkWidget* w, GdkEventMotion* ev, gpointer dummy)
        //        {
        //            int dx, dy;
        //            int pos_x, pos_y;
        //            unsigned int keys_buttons;

        //            /* Get the pointer position*/
        //            gdk_window_get_pointer(w->window, &pos_x, &pos_y, &keys_buttons);

        //            /* if we have a currents display, try to find the phase of the point under the cursor */
        //            if (structplot == SPcurrents && !dragging)
        //            {
        //                int i;
        //                double phase;
        //                int r, xc, yc, x, y;
        //                r = winsize / 10;
        //                xc = yc = r + r / 2 + 1;
        //                i = query_pixmap(gbackg, pos_x, pos_y);
        //                if (i >= 0)
        //                {
        //                    phase = ((i * 360.0 / NC_phase) - phaseoffset) * M_PI / 180.0;
        //                    x = xc + r * cos(phase);
        //                    y = yc + r * sin(phase);
        //                    gdk_window_clear(gdraw1->window);
        //                    gdk_gc_set_foreground(ggc, &c_axis);   /* black */
        //                    gdk_draw_line(gdraw1->window, ggc, xc, yc, x, y);
        //                }
        //                else if (i == -1)
        //                {
        //                    gdk_window_clear(gdraw1->window);
        //                }
        //            }

        //            /* Don't allow motion which is too small*/
        //            if (abs(pos_x - origx) + abs(pos_y - origy) < 2) return TRUE;

        //            /* Calculate motion vector */
        //            dx = pos_x - lastx;
        //            dy = pos_y - lasty;

        //            if (ev->state & GDK_BUTTON1_MASK)
        //            {
        //                /* drag with button1: rotate */
        //                theta -= 180.* (double)dy / (double)winsizey;
        //                if (theta < 0) theta = 0;
        //                if (theta > 180) theta = 180;
        //                phi -= 180.* (double)dx / (double)winsizex;
        //                scaleplot = 0;
        //                calc_vgain();
        //                if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //            }
        //            else if (ev->state & GDK_BUTTON2_MASK)
        //            {
        //                /* drag with button2: zoom in/out */
        //                double z;
        //                z = 5 * (double)(dx - dy) / winsize;
        //                zoom *= z + 1;
        //                trx -= (double)(origx - winsizex / 2) / winsize * z / zoom;
        //                try-= (double)(origy - winsizey / 2) / winsize * z / zoom;
        //                } else if (ev->state & GDK_BUTTON3_MASK)
        //                {
        //                    /* drag with button3: translate */
        //                    trx += (double)(dx) / zoom / winsize;
        //                    try+= (double)(dy) / zoom / winsize;
        //                    } else return TRUE;

        //                    quick = (ev->state & GDK_CONTROL_MASK);
        //                    upd_msg();
        //                    calcproj();
        //                    redraw = 1;

        //                    dragging = 1;

        //                    lastx = pos_x;
        //                    lasty = pos_y;

        //                    return TRUE;
        //                }


        //                gint keypress_event(GtkWidget* w, GdkEventKey* ev, gpointer dummy)
        //                {
        //                    int i;

        //                    switch (ev->keyval)
        //                    {
        //                        case GDK_Up:
        //                        case GDK_KP_Up:
        //                            theta += 5;
        //                            if (theta > 180) theta = 180;
        //                            upd_msg();
        //                            calcproj();
        //                            redraw = 1;
        //                            calc_vgain();
        //                            if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //                            break;
        //                        case GDK_Down:
        //                        case GDK_KP_Down:
        //                            theta -= 5;
        //                            if (theta < 0) theta = 0;
        //                            upd_msg();
        //                            calcproj();
        //                            redraw = 1;
        //                            calc_vgain();
        //                            if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //                            break;
        //                        case GDK_Left:
        //                        case GDK_KP_Left:
        //                            phi += 5;
        //                            upd_msg();
        //                            calcproj();
        //                            redraw = 1;
        //                            calc_vgain();
        //                            if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //                            break;
        //                        case GDK_Right:
        //                        case GDK_KP_Right:
        //                            phi -= 5;
        //                            upd_msg();
        //                            calcproj();
        //                            redraw = 1;
        //                            calc_vgain();
        //                            if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //                            break;
        //                        case GDK_Page_Down:
        //                        case GDK_KP_Page_Down:
        //                            i = rp_index;
        //                            do
        //                            {
        //                                i--;
        //                                if (i < 0) return TRUE;
        //                            } while (!neco[i].rp && !neco[i].cu && !neco[i].nf);
        //                            rp_index = i;
        //                            setrpfreq();
        //                            break;
        //                        case GDK_Page_Up:
        //                        case GDK_KP_Page_Up:
        //                            i = rp_index;
        //                            do
        //                            {
        //                                i++;
        //                                if (i >= numneco) return TRUE;
        //                            } while (!neco[i].rp && !neco[i].cu && !neco[i].nf);
        //                            rp_index = i;
        //                            setrpfreq();
        //                            break;
        //                        case '>':
        //                            animphase += M_PI / 8;
        //                            redraw = 1;
        //                            break;
        //                        case '<':
        //                            animphase -= M_PI / 8;
        //                            redraw = 1;
        //                            break;
        //                        case 'v':
        //                            switch (structplot)
        //                            {
        //                                case SPstruct: printf(" --struct"); break;
        //                                case SPtags: printf(" --tags"); break;
        //                                case SPcurrents: printf(" --currents"); break;
        //                                case SPanim: printf(" --animation"); break;
        //                            }
        //                            switch (gainplot)
        //                            {
        //                                case GPnearfield: printf(" --near"); break;
        //                                case GPslice: printf(" --slice"); break;
        //                                case GPframe: printf(" --frame"); break;
        //                                case GPopaque: printf(" --opaque"); break;
        //                            }
        //                            if (gainplot == GPslice || gainplot == GPframe || gainplot == GPopaque) switch (gainscale)
        //                                {
        //                                    case GSlinpower: printf(" --linpower"); break;
        //                                    case GSlinvolt: printf(" --linvolt"); break;
        //                                    case GSlog: printf(" --log"); break;
        //                                    case GSarrl: printf(" --arrl"); break;
        //                                }
        //                            if (gainplot == GPslice || gainplot == GPframe || gainplot == GPopaque || structplot == SPcurrents) switch (polarization)
        //                                {
        //                                    case POLnone: printf(" --pol=total"); break;
        //                                    case POLhor: printf(" --pol=hor"); break;
        //                                    case POLvert: printf(" --pol=vert"); break;
        //                                    case POLlhcp: printf(" --pol=lhcp"); break;
        //                                    case POLrhcp: printf(" --pol=rhcp"); break;
        //                                    case POLcolour: printf(" --pol=colour"); break;
        //                                }
        //                            printf("\n");
        //                            if (structplot != SPnone || gainplot != SPnone) printf(" --view=%g,%g,%g,%g,%g", phi, theta, zoom, trx,try);
        //        if (numneco > 1) printf(" --freq=%g", neco[rp_index].f);
        //        printf("\n");
        //        if (structplot == SPanim || gainplot == GPnearfield) printf(" --afreq=%g --aphase=%g", animfreq, animphase * 180 / M_PI);
        //        if (gainplot == GPnearfield) printf(" --escale=%g --hscale=%g%s", escale, hscale, show_p ? "" : " --hidepoynting");
        //        if (structplot == SPanim) printf(" --iscale=%g --qscale=%g", iscale, qscale);
        //        printf("\n");
        //        break;
        //      case 'z':
        //         if (animfreq == old_animfreq) break;
        //        if (animfreq == 0) { animfreq = old_animfreq; old_animfreq = 0; }
        //        else { old_animfreq = animfreq; animfreq = 0; }
        //        break;
        //        default:
        //         return FALSE;
        //    }
        //    return TRUE;
        //}





        //void cmd_reload(void)
        //{
        //    reread();
        //    redraw = redraw2 = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //    gtk_widget_grab_focus(gdraw2);
        //}

        //void cmd_X(void)
        //{
        //    phi = 0; theta = 90; redraw = 1;
        //    scaleplot = 1;
        //    calcproj();
        //    calc_vgain();
        //    if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //    upd_msg();
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_Y()
        //{
        //    phi = 90; theta = 90; redraw = 1;
        //    scaleplot = 1;
        //    calcproj();
        //    calc_vgain();
        //    if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //    upd_msg();
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_Z()
        //{
        //    phi = 0; theta = 0; redraw = 1;
        //    scaleplot = 1;
        //    calcproj();
        //    calc_vgain();
        //    if (plot2_vgain && redraw2 == 0) redraw2 = 2;
        //    upd_msg();
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //char* GSnames[] = { "lin.P", "lin.V", "arrl", "log", NULL };

        //void cmd_setgs(char** p)
        //{
        //    gainscale = p - GSnames;
        //    mark_gpo_invalid();
        //    if (rp_index >= 0)
        //    {
        //        extr /= GAINSIZE;
        //        process_nec_output(neco + rp_index);
        //    }
        //    upd_msg();
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //char* GPnames[] = { "none", "slice", "frame", "opaque", "near E/H", NULL };

        //void cmd_setgp(char** p)
        //{
        //    int o = gainplot;
        //    gainplot = p - GPnames;
        //    if (gainplot == GPnearfield)
        //    {
        //        if (o != GPnearfield && structplot != SPanim)
        //        {
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //            gtk_widget_show(botrow_anim);
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //        }
        //    }
        //    else if (o == GPnearfield && structplot != SPanim)
        //    {
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //        gtk_widget_hide(botrow_anim);
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //    }
        //    upd_msg();
        //    redraw = 1;
        //}

        //char* SPnames[] = { "none", "struct", "+tags", "currents", "anim.", NULL };

        //void cmd_setsp(char** p)
        //{
        //    int o = structplot;
        //    structplot = p - SPnames;
        //    if (structplot == SPcurrents)
        //    {
        //        if (o != SPcurrents)
        //        {
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //            gtk_widget_show(botrow_curr);
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //        }
        //    }
        //    else if (o == SPcurrents)
        //    {
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //        gtk_widget_hide(botrow_curr);
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //    }
        //    if (structplot == SPanim)
        //    {
        //        if (o != SPanim && gainplot != GPnearfield)
        //        {
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //            gtk_widget_show(botrow_anim);
        //            gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //        }
        //    }
        //    else if (o == SPanim && gainplot != GPnearfield)
        //    {
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, FALSE, FALSE, 0, GTK_PACK_START);
        //        gtk_widget_hide(botrow_anim);
        //        gtk_box_set_child_packing(GTK_BOX(vbox1), gdraw1, TRUE, TRUE, 0, GTK_PACK_START);
        //    }
        //    redraw = 1;
        //}


        //void cmd_togglelock(GtkWidget* w, gpointer* f)
        //{
        //    phasedirlock = GTK_TOGGLE_BUTTON(w)->active;
        //    if (phasedirlock)
        //    {
        //        phasephi = phi;
        //        phasetheta = theta;
        //    }
        //    else
        //    {
        //        if (structplot == SPcurrents)
        //        {
        //            calcproj();
        //            redraw = 1;
        //        }
        //    }
        //    gtk_widget_grab_focus(gdraw1);
        //}


        //char* POLnames[] = { "total", "hor.", "vert.", "lhcp", "rhcp", "colour", NULL };

        //void cmd_setpol(char** p)
        //{
        //    polarization = p - POLnames;
        //    mark_gpo_invalid();
        //    if (rp_index >= 0)
        //    {
        //        extr /= GAINSIZE;
        //        process_nec_output(neco + rp_index);
        //    }
        //    redraw = 1;
        //    calc_vgain();
        //    redraw2 = 1;
        //}

        //char* DCnames[] = { "no", "yes", NULL };

        //void cmd_setdc(char** p)
        //{
        //    distcor = p - DCnames;
        //    redraw = 1;
        //}


        //void cmd_phaseoffset(GtkAdjustment* adj)
        //{
        //    phaseoffset = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_maxthickness(GtkAdjustment* adj)
        //{
        //    maxthickness = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_escale(GtkAdjustment* adj)
        //{
        //    escale = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_hscale(GtkAdjustment* adj)
        //{
        //    hscale = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_qscale(GtkAdjustment* adj)
        //{
        //    qscale = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_iscale(GtkAdjustment* adj)
        //{
        //    iscale = adj->value;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void cmd_animfreq(GtkAdjustment* adj)
        //{
        //    animfreq = adj->value;
        //    gtk_widget_grab_focus(gdraw1);
        //}


        //void cmd_export_eps_ok(GtkWidget* w, GtkFileSelection* fs)
        //{
        //    const char* s;
        //    s = gtk_file_selection_get_filename(fs);
        //    if (write_postscript(s, draw_all, winsizex, winsizey)) fprintf(stderr, "Error writing postscript\n");
        //    gtk_widget_destroy(GTK_WIDGET(fs));
        //}

        //#ifdef HAVE_LIBPNG
        //void cmd_export_png_ok(GtkWidget* w, GtkFileSelection* fs)
        //{
        //    const char* s;
        //    s = gtk_file_selection_get_filename(fs);
        //    if (write_png(1, s)) fprintf(stderr, "Error writing PNG\n");
        //    gtk_widget_destroy(GTK_WIDGET(fs));
        //}
        //#endif


        //int do_animation(void)
        //{
        //    if (animfreq == 0) return 1;
        //    if (animbusy) return 1;
        //    if (structplot != SPanim && gainplot != GPnearfield) return 1;
        //    animphase += 0.002 * anim_update_interval * M_PI * animfreq;
        //    if (animphase > 2 * M_PI) animphase -= 2 * M_PI;
        //    redraw = 1;
        //    animbusy = 1;
        //    return 1;
        //}


        ///* -------------------- event handlers (and their supporting code) for window 2 (frequency plot) ---------------------- */


        //gint expose_event2(GtkWidget* w, GdkEventExpose* ev, gpointer dummy)
        //{
        //    int x;
        //    if (ev->count != 0) return TRUE;
        //    if (rp_index < 0) return TRUE;
        //    x = freqx(neco[rp_index].f);
        //    gdk_gc_set_foreground(ggc2, &c_back);
        //    gdk_draw_line(gdraw2->window, ggc2, x, 0, x, win2sizey);
        //    return TRUE;
        //}


        //gint resize_event2(GtkWidget* w, GdkEventConfigure* ev, gpointer dummy)
        //{
        //    if (!gbackg2) return FALSE;
        //    win2sizex = ev->width;
        //    win2sizey = ev->height;
        //    gdk_pixmap_unref(gbackg2);
        //    gbackg2 = gdk_pixmap_new(w->window, win2sizex, win2sizey, gdk_drawable_get_depth(w->window));

        //    redraw2 = 1;
        //    return TRUE;
        //}



        //gint buttonpress_event2(GtkWidget* w, GdkEventButton* ev, gpointer dummy)
        //{
        //    int i;
        //    i = freqindex(xfreq(ev->x));
        //    if (i < 0) return TRUE;
        //    if (rp_index == i) return TRUE;
        //    rp_index = i;
        //    setrpfreq();
        //    return TRUE;
        //}



        //gint motion_event2(GtkWidget* w, GdkEventMotion* ev, gpointer dummy)
        //{
        //    int pos_x, pos_y;
        //    unsigned int keys_buttons;

        //    /* Get the pointer position*/
        //    gdk_window_get_pointer(w->window, &pos_x, &pos_y, &keys_buttons);

        //    ev->x = pos_x;
        //    ev->y = pos_y;
        //    buttonpress_event2(w, (GdkEventButton*)ev, dummy);
        //    return TRUE;
        //}


        //gint keypress_event2(GtkWidget* w, GdkEventKey* ev, gpointer dummy)
        //{
        //    int i;

        //    switch (ev->keyval)
        //    {
        //        case GDK_Page_Down:
        //        case GDK_KP_Page_Down:
        //        case GDK_Down:
        //        case GDK_KP_Down:
        //        case GDK_Left:
        //        case GDK_KP_Left:
        //            i = rp_index;
        //            do
        //            {
        //                i--;
        //                if (i < 0) return TRUE;
        //            } while (neco[i].rp == NULL);
        //            rp_index = i;
        //            setrpfreq();
        //            break;
        //        case GDK_Page_Up:
        //        case GDK_KP_Page_Up:
        //        case GDK_Up:
        //        case GDK_KP_Up:
        //        case GDK_Right:
        //        case GDK_KP_Right:
        //            i = rp_index;
        //            do
        //            {
        //                i++;
        //                if (i >= numneco) return TRUE;
        //            } while (neco[i].rp == NULL);
        //            rp_index = i;
        //            setrpfreq();
        //            break;
        //        default:
        //            return FALSE;
        //    }
        //    return TRUE;
        //}



        //void cmd_export2_eps_ok(GtkWidget* w, GtkFileSelection* fs)
        //{
        //    const char* s;
        //    s = gtk_file_selection_get_filename(fs);
        //    if (write_postscript(s, draw_all2, win2sizex, win2sizey)) fprintf(stderr, "Error writing postscript\n");
        //    gtk_widget_destroy(GTK_WIDGET(fs));
        //}

        //#ifdef HAVE_LIBPNG
        //void cmd_export2_png_ok(GtkWidget* w, GtkFileSelection* fs)
        //{
        //    const char* s;
        //    s = gtk_file_selection_get_filename(fs);
        //    if (write_png(2, s)) fprintf(stderr, "Error writing PNG\n");
        //    gtk_widget_destroy(GTK_WIDGET(fs));
        //}
        //#endif

        //void cmd_export(GtkWidget* ww, void(* f)(), char* ext)
        //{
        //    GtkWidget* w;
        //    char* s,*p;

        //    w = gtk_file_selection_new("Export PostScript");
        //    gtk_signal_connect(GTK_OBJECT(GTK_FILE_SELECTION(w)->ok_button), "clicked",
        //       (GtkSignalFunc)f, w);
        //    gtk_signal_connect_object(GTK_OBJECT(GTK_FILE_SELECTION(w)->cancel_button), "clicked",
        //       (GtkSignalFunc)gtk_widget_destroy, GTK_OBJECT(w));

        //    s = mymalloc(strlen(inputfilename) + 8);
        //    strcpy(s, inputfilename);
        //    p = strrchr(s, '.');
        //    if (!p) p = s + strlen(s);
        //    strcpy(p, ext);
        //    gtk_file_selection_set_filename(GTK_FILE_SELECTION(w), s);
        //    free(s);

        //    gtk_widget_show(w);
        //    gtk_widget_grab_focus(gdraw1);
        //    gtk_widget_grab_focus(gdraw2);
        //}

        //void cmd_export_eps(GtkWidget* ww, void(* f)())
        //{
        //    cmd_export(ww, f, ".eps");
        //}

        //#ifdef HAVE_LIBPNG
        //void cmd_export_png(GtkWidget* ww, void(* f)())
        //{
        //    cmd_export(ww, f, ".png");
        //}
        //#endif



        //void cmd_setZ0(GtkEntry* w)
        //{
        //    double r;
        //    const char* q,*p;
        //    char s[20];

        //    q = gtk_entry_get_text(w);
        //    p = strchr(q, '=');
        //    if (p) p++;
        //    else p = q;
        //    r = atof(p);
        //    if (r > 0)
        //    {
        //        int i;
        //        r0 = r;
        //        for (i = 0; i < numneco; i++) calcswr(neco + i);
        //        if (plot2_swr) redraw2 = 1;
        //    }
        //    sprintf(s, "Z0=%g", r0);
        //    gtk_entry_set_text(w, s);
        //}


        ///* -------------------- X and GTK initialisation ------------------------------------------------ */

        //void getcolor(char* name, GdkColor* xc)
        //{
        //    if (depth == 1)
        //    {          /* on 1bpp displays, choose black for everything except the background */
        //        if (strcmp(name, C_BG)) name = "black";
        //        else name = "white";
        //    }
        //    gdk_color_parse(name, xc);
        //    if (!gdk_colormap_alloc_color(gcm, xc, TRUE, FALSE))
        //    {
        //        /* if allocation failed, use black */
        //        gdk_color_black(gcm, xc);
        //    }
        //}


        ///* Function for conversion from hue/lightness/saturation to red/green/blue.
        //   This function has been copied literally from GTK+ (gtkstyle.c).
        //   However, I doubt that this has the property of "perceptual uniformity";
        //   see e.g. Poynton's Color FAQ at
        //   http://www.inforamp.net/~poynton/notes/colour_and_gamma/ColorFAQ.html .
        //*/
        //static void
        //hls_to_rgb(gdouble* h,
        //            gdouble* l,
        //            gdouble* s)
        //{
        //    gdouble hue;
        //    gdouble lightness;
        //    gdouble saturation;
        //    gdouble m1, m2;
        //    gdouble r, g, b;

        //    lightness = *l;
        //    saturation = *s;

        //    if (lightness <= 0.5)
        //        m2 = lightness * (1 + saturation);
        //    else
        //        m2 = lightness + saturation - lightness * saturation;
        //    m1 = 2 * lightness - m2;

        //    if (saturation == 0)
        //    {
        //        *h = lightness;
        //        *l = lightness;
        //        *s = lightness;
        //    }
        //    else
        //    {
        //        hue = *h + 120;
        //        while (hue > 360)
        //            hue -= 360;
        //        while (hue < 0)
        //            hue += 360;

        //        if (hue < 60)
        //            r = m1 + (m2 - m1) * hue / 60;
        //        else if (hue < 180)
        //            r = m2;
        //        else if (hue < 240)
        //            r = m1 + (m2 - m1) * (240 - hue) / 60;
        //        else
        //            r = m1;

        //        hue = *h;
        //        while (hue > 360)
        //            hue -= 360;
        //        while (hue < 0)
        //            hue += 360;

        //        if (hue < 60)
        //            g = m1 + (m2 - m1) * hue / 60;
        //        else if (hue < 180)
        //            g = m2;
        //        else if (hue < 240)
        //            g = m1 + (m2 - m1) * (240 - hue) / 60;
        //        else
        //            g = m1;

        //        hue = *h - 120;
        //        while (hue > 360)
        //            hue -= 360;
        //        while (hue < 0)
        //            hue += 360;

        //        if (hue < 60)
        //            b = m1 + (m2 - m1) * hue / 60;
        //        else if (hue < 180)
        //            b = m2;
        //        else if (hue < 240)
        //            b = m1 + (m2 - m1) * (240 - hue) / 60;
        //        else
        //            b = m1;

        //        *h = r;
        //        *l = g;
        //        *s = b;
        //    }
        //}

        //void get_currents_colors(void)
        //{
        //    int i;
        //    double r, g, b;

        //    for (i = 0; i < NC_phase; i++)
        //    {
        //        r = i * 360./ NC_phase;
        //        g = 0.5;
        //        b = 1;
        //        hls_to_rgb(&r, &g, &b);
        //        c_currents[i].red = r * 65535;
        //        c_currents[i].green = g * 65535;
        //        c_currents[i].blue = b * 65535;
        //        gdk_colormap_alloc_color(gcm, c_currents + i, TRUE, TRUE);
        //    }
        //}


        //GtkWidget* make_option_menu(char** list, void (* f)(char**), int v)
        //{
        //    GtkWidget* menu, *button;
        //    GtkWidget* w;
        //    char** p;

        //    menu = gtk_menu_new();
        //    p = list;
        //    while (*p)
        //    {
        //        w = gtk_menu_item_new_with_label(*p);
        //        gtk_menu_append(GTK_MENU(menu), w);
        //        gtk_widget_show(w);
        //        gtk_signal_connect_object(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(f), (gpointer)(p));
        //        p++;
        //    }
        //    gtk_menu_set_active(GTK_MENU(menu), v);

        //    button = gtk_option_menu_new();
        //    gtk_option_menu_set_menu(GTK_OPTION_MENU(button), menu);

        //    return button;
        //}

        //void cmd_quit(GtkWidget* w, gpointer* dummy)
        //{
        //    redraw = redraw2 = 0;
        //    quit = 1;
        //}

        //void toggle1helper(GtkWidget* w, gpointer* f)
        //{
        //    *(int*)f = GTK_TOGGLE_BUTTON(w)->active;
        //    redraw = 1;
        //    gtk_widget_grab_focus(gdraw1);
        //}

        //void toggle2helper(GtkWidget* w, gpointer* f)
        //{
        //    *(int*)f = GTK_TOGGLE_BUTTON(w)->active;
        //    redraw2 = 1;
        //    gtk_widget_grab_focus(gdraw2);
        //}


        //void maininitX(int really)
        //{
        //    GtkWidget* win1, *win2;
        //    GtkWidget* vbox, *toprow;
        //    GtkWidget* w;
        //    GtkAccelGroup* acc;
        //    GdkBitmap* icon;
        //    GtkTooltips* tt;
        //    GtkObject* adj;

        //    tt = gtk_tooltips_new();
        //#define SET_TT(s) gtk_tooltips_set_tip(tt,w,s,NULL)

        //    {
        //        GdkVisual* vi;
        //        vi = gdk_visual_get_best();
        //        depth = vi->depth;
        //        gcm = gdk_colormap_new(vi, FALSE);
        //    }

        

        

        //    /* ----------- initialize some variables etc. that are needed for direct Xlib drawing ------------ */

        //    getcolor("white", &c_bg);
        //    getcolor(C_AXIS, &c_axis);
        //    getcolor(C_WIRE, &c_wire);
        //    getcolor(C_SURF, &c_surf);
        //    getcolor(C_BACK, &c_back);
        //    getcolor(C_GAIN, &c_gain);
        //    getcolor(C_SCALE, &c_scale);
        //    getcolor(C_EXCI, &c_exci);
        //    getcolor(C_LOAD, &c_load);
        //    getcolor(C_NETW, &c_netw);
        //    getcolor(C_INACTIVE, &c_inactive);
        //    getcolor(C_EFIELD, &c_efield);
        //    getcolor(C_HFIELD, &c_hfield);
        //    getcolor(C_POYNTING, &c_poynting);
        //    getcolor(C_QPOS, &c_qpos);
        //    getcolor(C_QNEG, &c_qneg);
        //    getcolor(C_GAIN_LIN, &c_gain_lin);
        //    getcolor(C_GAIN_LHCP, &c_gain_lhcp);
        //    getcolor(C_GAIN_RHCP, &c_gain_rhcp);
        //    get_currents_colors();

        //    gfont = gdk_font_load(XFONT);
        //    if (!gfont) gfont = gdk_font_load("fixed");
        //    if (!gfont) { fprintf(stderr, "Can't load fonts '%s' or 'fixed'.\n", XFONT); exit(1); }
        //    fontheight = gfont->ascent + gfont->descent;

        //    if (window1open) gdk_gc_set_font(ggc, gfont);
        //    if (window2open) gdk_gc_set_font(ggc2, gfont);

        //    if (gdraw1 == NULL) gdraw1 = gdraw2;
        //    else if (gdraw2 == NULL) gdraw2 = gdraw1;

        //    g_timeout_add(anim_update_interval, (GtkFunction)do_animation, NULL);
        //}


        ///* -------------------- main loop (wait for event, process event, redisplay) ------------ */

        //int* present_redraw = NULL;

        public static int Pending()
        {
            //    if (present_redraw == NULL) return 0;
            //    gtk_main_iteration_do(FALSE);
            //    return *present_redraw | quit;
            return 0;       //TODO return real value
        }


        //void mainX(int really)
        //{
        //    maininitX(really);
        //    if (!really) return;
        //    while (!quit)
        //    {
        //        if (window1open && redraw)
        //        {
        //            present_redraw = &redraw;
        //            redraw = 0;
        //         out = &outX;
        //            draw_all(1);
        //        }
        //        if (window2open && redraw2)
        //        {
        //            int r = redraw2;
        //            present_redraw = &redraw2;
        //            redraw2 = 0;
        //         out = &outX2;
        //            draw_all2(r == 2);
        //            if (rp_index >= 0)
        //            {
        //                int x;
        //                x = freqx(neco[rp_index].f);
        //                gdk_gc_set_foreground(ggc2, &c_back);
        //                gdk_draw_line(gdraw2->window, ggc2, x, 0, x, win2sizey);
        //            }
        //        }
        //        animbusy = 0;
        //        gtk_main_iteration_do(!((window1open && redraw) || (window2open && redraw2)));
        //    }
        //}

        //void initX(int* argcp, char** argv)
        //{
        //    gtk_init(argcp, &argv);
        //    winsize = winsizex = winsizey = ini_winsize;
        //    win2sizex = win2sizey = ini_winsize;
        //}
    }

    class Xwin2 : IOutdevInterface
    {
        public void SetLineAttributes(uint a, int b, int c, int d)
        {
            //gdk_gc_set_line_attributes(ggc2, a, b, c, d);
        }

        public void DrawLine(double a, double b, double c, double d)
        {
            //            if (b < -32000 || d < -32000 || b > 32000 || d > 32000) return;
            //            gdk_draw_line(gbackg2, ggc2, (int)(a + 0.5), (int)(b + 0.5), (int)(c + 0.5), (int)(d + 0.5));
        }

        public void SetForeground(Color xc)
        {
            //gdk_gc_set_foreground(ggc2, xc);
        }

        public void SetClipRectangle(double x1, double y1, double x2, double y2)
        {
            //            GdkRectangle gr;
            //            gr.x = (int)(x1 + 0.5);
            //            gr.y = (int)(y1 + 0.5);
            //            gr.width = (int)(x2 - gr.x + 1.5);
            //            gr.height = (int)(y2 - gr.y + 1.5);
            //            gdk_gc_set_clip_rectangle(ggc2, &gr);
        }

        public void ClearWindow()
        {
            //            X2_SetClipRectangle(0, 0, win2sizex, win2sizey);
            //            X2_SetForeground(&c_bg);
            //            gdk_draw_rectangle(gbackg2, ggc2, TRUE, 0, 0, win2sizex, win2sizey);
        }

        public void ClearRectangle(double x1, double y1, double x2, double y2)
        {
            //            X2_SetClipRectangle(0, 0, win2sizex, win2sizey);
            //            X2_SetForeground(&c_bg);
            //            gdk_draw_rectangle(gbackg2, ggc2, TRUE, (int)(x1 + 0.5), (int)(y1 + 0.5), (int)(x2 - x1 + 1.5), (int)(y2 - y1 + 1.5));
        }


        public void DrawString(double a, double b, string s, double d, double e)    /* draw string */
        {
            //            if (d > 0) a -= d * gdk_string_width(gfont, s);
            //            b += -gfont->descent + e * (gfont->descent + gfont->ascent);
            //            gdk_draw_string(gbackg2, gfont, ggc2, (int)(a + 0.5), (int)(b + 0.5), s);
        }

        public void Complete()
        {
            //            gdk_window_set_back_pixmap(gdraw2->window, gbackg2, 0);
            //            gdk_window_clear(gdraw2->window);
        }

        //        Outdev outX2 ={
        //   X2_SetLineAttributes,
        //   X2_DrawLine,
        //   X2_SetForeground,
        //   X2_ClearWindow,
        //   X2_DrawString,
        //   X2_Complete,
        //   X2_SetClipRectangle,
        //   X2_ClearRectangle,
        //};
    }
}