using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Necview
{
    /// <summary>
    /// window 2 (plots of gain, SWR etc. vs frequency)
    /// </summary>
    public partial class PlotVsFrequency : Form
    {
        /// <summary>
        /// initialize window 2 (plots of gain, SWR etc. vs frequency)
        /// </summary>
        public PlotVsFrequency()
        {
            InitializeComponent();

            //    /* -----------  ----------------------------------------------- */

            //        acc = gtk_accel_group_new();

            //        toprow = gtk_hbox_new(FALSE, 0);

            //        w = gtk_button_new_with_label("reload"); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
            //        gtk_signal_connect(GTK_OBJECT(w), "clicked", GTK_SIGNAL_FUNC(cmd_reload), (gpointer)NULL);
            //        gtk_widget_add_accelerator(w, "clicked", acc, GDK_R, 0, 0);
            //        gtk_widget_add_accelerator(w, "clicked", acc, GDK_period, 0, 0);
            //        SET_TT("Reread the data files");

            //        {
            //            /* create the 'export' menu; and to get that one into the top bar, we need to create a menu_bar too */
            //            GtkWidget* m,*v;

            //            v = gtk_menu_item_new_with_label("export"); gtk_widget_show(v);
            //            m = gtk_menu_new(); gtk_menu_item_set_submenu(GTK_MENU_ITEM(v), m); gtk_widget_show(m);
            //            w = gtk_menu_item_new_with_label("PostScript"); gtk_menu_append(GTK_MENU(m), w); gtk_widget_show(w);
            //            gtk_signal_connect(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(cmd_export_eps), (gpointer)cmd_export2_eps_ok);
            //# ifdef HAVE_LIBPNG
            //            w = gtk_menu_item_new_with_label("PNG"); gtk_menu_append(GTK_MENU(m), w); gtk_widget_show(w);
            //            gtk_signal_connect(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(cmd_export_png), (gpointer)cmd_export2_png_ok);
            //#endif
            //            w = gtk_menu_bar_new(); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
            //            gtk_menu_bar_append(GTK_MENU_BAR(w), v);
            //        }

            //        w = gtk_entry_new_with_max_length(7); gtk_box_pack_start(GTK_BOX(toprow), w, TRUE, TRUE, 0); gtk_widget_show(w);
            //        gtk_widget_set_usize(w, 70, 0);
            //        gtk_signal_connect(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(cmd_setZ0), (gpointer)NULL);
            //        gtk_signal_connect(GTK_OBJECT(w), "leave_notify_event", GTK_SIGNAL_FUNC(cmd_setZ0), NULL);
            //        cmd_setZ0(GTK_ENTRY(w));
            //        SET_TT("Reference impedance for SWR");

            //#define ADDTOGGLE2(label,p,s)  \
            //        w = gtk_toggle_button_new_with_label(label); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);  \
            //   gtk_signal_connect(GTK_OBJECT(w), "toggled", GTK_SIGNAL_FUNC(toggle2helper), (gpointer)p);  \
            //   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), *p);  \
            //   SET_TT(s);

            //        gdraw2 = w;  /* necessary since in the following toggle2helper() will be called, so gdraw2 should be valid */
            //        ADDTOGGLE2("vgain", &plot2_vgain, "Toggle plot of gain toward viewer");
            //        ADDTOGGLE2("maxgain", &plot2_maxgain, "Toggle plot of maximum gain");
            //        ADDTOGGLE2("dir", &plot2_dir, "Toggle plot of direction of maximum gain");
            //        ADDTOGGLE2("SWR", &plot2_swr, "Toggle plot of SWR");
            //        ADDTOGGLE2("Re/Im", &plot2_z, "Toggle plot of impedance (Re/Im)");
            //        ADDTOGGLE2("phi/abs", &plot2_z2, "Toggle plot of impedance (phi/abs)");

            //        vbox = gtk_vbox_new(FALSE, 0);
            //        gtk_box_pack_start(GTK_BOX(vbox), toprow, FALSE, FALSE, 0);
            //        gdraw2 = w = gtk_drawing_area_new();
            //        gtk_drawing_area_size(GTK_DRAWING_AREA(w), win2sizex, win2sizey);
            //        gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 0);
            //        gtk_widget_show(w);
            //        GTK_WIDGET_SET_FLAGS(w, GTK_CAN_FOCUS);

            //        gtk_signal_connect(GTK_OBJECT(w), "configure_event", GTK_SIGNAL_FUNC(resize_event2), NULL);
            //        gtk_signal_connect(GTK_OBJECT(w), "button_press_event", GTK_SIGNAL_FUNC(buttonpress_event2), NULL);
            //        gtk_signal_connect(GTK_OBJECT(w), "motion_notify_event", GTK_SIGNAL_FUNC(motion_event2), NULL);
            //        gtk_signal_connect_after(GTK_OBJECT(w), "key_press_event", GTK_SIGNAL_FUNC(keypress_event2), NULL);
            //        gtk_widget_set_events(w, GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK |
            //           GDK_BUTTON_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK |
            //           GDK_KEY_PRESS_MASK);

            //        gtk_container_add(GTK_CONTAINER(win2), vbox);
            //        gtk_widget_show(toprow);
            //        gtk_widget_show(vbox);
            //        if (really) gtk_widget_show(win2);
            //        else gtk_widget_realize(gdraw2);

            //        gtk_widget_grab_focus(gdraw2);

            //        gtk_window_add_accel_group(GTK_WINDOW(win2), acc);

            //        ggc2 = gdk_gc_new(w->window);

            //        gbackg2 = gdk_pixmap_new(w->window, win2sizex, win2sizey, gdk_drawable_get_depth(w->window));

            //        gtk_signal_connect(GTK_OBJECT(w), "expose_event", GTK_SIGNAL_FUNC(expose_event2), NULL);

            //        icon = gdk_bitmap_create_from_data(win2->window, icon_bits, icon_width, icon_height);
            //        gdk_window_set_icon(win2->window, NULL, icon, NULL);
        }

        int win2sizex, win2sizey;       // size of window in pixels

        bool plot2_swr = true;     // show the SWR graph?
        bool plot2_maxgain = true;     // show the maxgain and front/back graph?
        bool plot2_vgain = false;    // show the vgain graph?
        bool plot2_z = false;    // show the impedance graph?
        bool plot2_z2 = false;    // show the phi(z)/abs(z) graph?
        bool plot2_dir = false;    // show the direction-of-maximum-gain graph?

        double r0 = Necview.Config.R0;  // reference impedance for SWR calculation

        /// <summary>
        /// mi and ma point to minimum and maximum value of a range of values to
        /// be plotted, and np to the maximum acceptable number of subdivision.
        /// This function tries to modify the minimum and maximum and the number
        /// of subdivision such that the resulting grid lines are at "round" numbers.
        /// </summary>
        /// <param name="mi"></param>
        /// <param name="ma"></param>
        /// <param name="np"></param>
        void fixrange1(ref double mi, ref double ma, ref int np)
        {
            double d, e;
            double a;
            double newmin, newmax;
            int i;
            int n = np;
            double[] acceptable = { 10, 5.0, 2.0, 1.0, -1 };

            if (ma == mi)
            {
                if (mi > 0)
                {
                    mi = 0;
                    ma = 2 * ma;
                }
                else if (mi < 0)
                {
                    mi = 2 * mi;
                    ma = 0;
                }
                else
                {
                    mi = -10;
                    ma = 10;
                }
            }
            d = (ma - mi) / n;
            e = 1.0;
            while (e < d)
                e *= 10;
            while (e > d)
                e /= 10;
            a = d / e;
            i = 0;
            while (acceptable[i] > a)
                i++;
            if (acceptable[i] == -1)
                i--;
            i++;
            do
            {
                i--;
                if (i < 0)
                {
                    e *= 10;
                    i = 0;
                    while (acceptable[i + 1] > 0)
                        i++;
                }
                a = acceptable[i];
                d = a * e;
                newmin = d * Math.Floor(mi / d);
                newmax = d * Math.Ceiling(ma / d);
                n = (int)((newmax - newmin) / d + 0.5);
            } while (n > np);
            np = n;
            mi = newmin;
            ma = newmax;
        }

        /// <summary>
        /// like fixrange2(), but for two (vertical) axes simultaneously
        /// </summary>
        /// <param name="mi1"></param>
        /// <param name="ma1"></param>
        /// <param name="mi2"></param>
        /// <param name="ma2"></param>
        /// <param name="np"></param>
        void fixrange2(ref double mi1, ref double ma1, ref double mi2, ref double ma2, ref int np)
        {
            double[] acceptable = { 100.0, 50.0, 25.0, 20.0, 10.0, 5.0, 2.5, 2.0, 1.0, 0.5, 0.25, 0.2, 0.1, 0.05, 0.025, 0.02, 0.01, -1 };
            double a, d, e1, e2, s;
            int n = np;
            int i1, i2;
            int i, j;
            int ibest, jbest;
            int[] n1 = new int[5];
            int[] n2 = new int[5];
            double[] x1 = new double[4];
            double[] x2 = new double[4];
            double best;

            if (ma1 == mi1)
            {
                if (mi1 > 0)
                {
                    mi1 = 0;
                    ma1 = 2 * ma1;
                }
                else if (mi1 < 0)
                {
                    mi1 = 2 * mi1;
                    ma1 = 0;
                }
                else
                {
                    mi1 = -10;
                    ma1 = 10;
                }
            }
            d = (ma1 - mi1) / n;    // d is the ideal, but usually not acceptable, stepsize, for axis 1 
            ma1 -= 0.00001 * d;     // prevent rounding errors from causing a boundary of say 1000 to be seen 
                                    // as slightly larger than say 10 steps of 100 each
            mi1 += 0.00001 * d;     // idem
            d -= 0.00001 * d;
            e1 = 1.0;
            while (e1 < d)
                e1 *= 10;
            while (e1 > d)
                e1 /= 10;           // e1 is the appropriate power of 10 to scale the steps, for axis 1
            a = d / e1;
            i1 = 0;
            while (acceptable[i1 + 1] >= a)
                i1++;               // i1 is the index in the acceptable[] array of the highest acceptable stepsize, for axis 1
            for (i = 0; i < 4; i++)
            {
                // consider this and the next 3 lower stepsizes:
                s = e1 * acceptable[i1 - i];
                n1[i] = (int)(Math.Ceiling(ma1 / s) - Math.Floor(mi1 / s));     // minimum number of acceptable steps
                x1[i] = (ma1 - mi1) / s;                                        // "usage factor": how many of these steps does the data cover?
            }

            // same calculations for axis 2
            if (ma2 == mi2)
            {
                if (mi2 > 0)
                {
                    mi2 = 0;
                    ma2 = 2 * ma2;
                }
                else if (mi2 < 0)
                {
                    mi2 = 2 * mi2;
                    ma2 = 0;
                }
                else
                {
                    mi2 = -10;
                    ma2 = 10;
                }
            }
            d = (ma2 - mi2) / n;
            ma2 -= 0.00001 * d;
            mi2 += 0.00001 * d;
            d -= 0.00001 * d;
            e2 = 1.0;
            while (e2 < d)
                e2 *= 10;
            while (e2 > d)
                e2 /= 10;
            a = d / e2;
            i2 = 0;
            while (acceptable[i2 + 1] >= a)
                i2++;
            for (i = 0; i < 4; i++)
            {
                s = e2 * acceptable[i2 - i];
                n2[i] = (int)(Math.Ceiling(ma2 / s) - Math.Floor(mi2 / s));
                x2[i] = (ma2 - mi2) / s;
            }

            // search for best combination: the combination for which 
            // the data covers as large a fraction of both axes as possible
            best = 0;
            ibest = jbest = 0;
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    double x;
                    int nn;
                    nn = n1[i];
                    if (n2[j] > nn) nn = n2[j];
                    x = (x1[i] / nn) * (x2[j] / nn);
                    if (x > best * 1.1 || (x > best && nn >= np))
                    {
                        best = x;
                        ibest = i;
                        jbest = j;
                        np = nn;
                    }
                }

            n = np;
            i1 -= ibest;
            i2 -= jbest;
            s = e1 * acceptable[i1];
            mi1 = s * Math.Floor(mi1 / s);
            ma1 = mi1 + n * s;
            s = e2 * acceptable[i2];
            mi2 = s * Math.Floor(mi2 / s);
            ma2 = mi2 + n * s;
        }


        double minf, maxf;
        int xleft, xright;

        //#define idxOK(idx,ne) (idx>=0 && ( !ONLY_IF_RP(idx) || ne->rp ) && ne->d[idx]>-DBL_MAX)

        void freqplot(
           int idx1,                        /* index in neco[].d[] of quantity for left axis */
           int idx2,                        /* index in neco[].d[] of quantity for right axis */
           int idx1a,                       /* index in neco[].d[] of second quantity for left axis (dotted line) */
           int idx2a,                       /* index in neco[].d[] of second quantity for right axis (dotted line) */
           string title1, string title2,      /* titles for left and right */
           string title,                     /* center title */
           Color color1, Color color2,  /* colours for both curves */
           double ybotf, double ytopf       /* vertical position; 0...1 = top...bottom of window */
        )
        {
            int ybot, ytop;
            int i;
            double min1, max1, min2, max2;
            //            NECoutput* ne;
            int ntx, nty;
            int xx1, xx2, yy1, yy2;
            int xx1a, xx2a, yy1a, yy2a;

            // choose the corner points of the graph area
            //            ybot = ybotf * win2sizey - fontheight;
            //            ytop = ytopf * win2sizey + fontheight;
            //            xleft = 5 * fontheight;
            //            xright = win2sizex - 5 * fontheight;

            //            /* find the ranges */
            //            minf = maxf = neco[0].f;
            //            min1 = min2 = DBL_MAX;
            //            max1 = max2 = -DBL_MAX;
            //            for (i = 0, ne = neco; i < numneco; i++, ne++)
            //            {
            //                if (ne->f < minf) minf = ne->f;
            //                if (ne->f > maxf) maxf = ne->f;
            //                if (idxOK(idx1, ne))
            //                {
            //                    if (ne->d[idx1] < min1) min1 = ne->d[idx1];
            //                    if (ne->d[idx1] > max1) max1 = ne->d[idx1];
            //                }
            //                if (idxOK(idx1a, ne))
            //                {
            //                    if (ne->d[idx1a] < min1) min1 = ne->d[idx1a];
            //                    if (ne->d[idx1a] > max1) max1 = ne->d[idx1a];
            //                }
            //                if (idxOK(idx2, ne))
            //                {
            //                    if (ne->d[idx2] < min2) min2 = ne->d[idx2];
            //                    if (ne->d[idx2] > max2) max2 = ne->d[idx2];
            //                }
            //                if (idxOK(idx2a, ne))
            //                {
            //                    if (ne->d[idx2a] < min2) min2 = ne->d[idx2a];
            //                    if (ne->d[idx2a] > max2) max2 = ne->d[idx2a];
            //                }
            //            }
            //            if (min1 > max1) { idx1 = -1; idx1a = -1; }
            //            if (min2 > max2) { idx2 = -1; idx2a = -1; }

            //            /* extend the ranges to have 'round' numbers at each division */
            //            ntx = (xright - xleft) / fontheight / 3;
            //            fixrange1(&minf, &maxf, &ntx);
            //            nty = (ybot - ytop) / fontheight;
            //            if (idx1 >= 0)
            //            {
            //                if (idx1 == neco_zr)
            //                {
            //                    if (max1 > 20 * r0) max1 = 20 * r0;
            //                }
            //            }
            //            if (idx2 >= 0)
            //            {
            //                if (idx2 == neco_zi || idx2 == neco_zabs)
            //                {
            //                    if (max2 > 20 * r0 && min2 < 20 * r0) max2 = 20 * r0;
            //                    if (min2 < -20 * r0 && max2 > -20 * r0) min2 = -20 * r0;
            //                }
            //            }
            //            if (idx1 == neco_swr) { min1 = 0; if (max1 > 10) max1 = 9; else max1 -= 1; }
            //            if (idx2 == neco_swr) { min2 = 0; if (max2 > 10) max2 = 9; else max2 -= 1; }
            //            if (idx1 < 0 && idx1a < 0)
            //            {
            //                if (idx2 < 0 && idx2a < 0) return;
            //                fixrange1(&min2, &max2, &nty);
            //            }
            //            else
            //            {
            //                if (idx2 < 0 && idx2a < 0) fixrange1(&min1, &max1, &nty);
            //                else fixrange2(&min1, &max1, &min2, &max2, &nty);
            //            }
            //            if (idx1 == neco_swr) { min1 += 1; max1 += 1; }
            //            if (idx2 == neco_swr) { min2 += 1; max2 += 1; }



            // macros for converting from "real" values to screen coordinates
            //#define sx(f) (((f)-minf)/(maxf-minf)*(xright-xleft)+xleft)
            //#define sy1(f) (((f)-min1)/(max1-min1)*(ytop-ybot)+ybot)
            //#define sy2(f) (((f)-min2)/(max2-min2)*(ytop-ybot)+ybot)

            //            SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //            /* vertical grid lines and associated labels */
            //            for (i = 0; i <= ntx; i++)
            //            {
            //                double f;
            //                int x;
            //                char s[20];
            //                f = minf + (maxf - minf) * ((double)i) / ntx;
            //                x = sx(f);
            //                if (i > 0 && i < ntx)
            //                {
            //                    SetForeground(&c_scale);
            //                    DrawLine(x, ybot, x, ytop);
            //                }
            //                SetForeground(&c_axis);
            //                sprintf(s, "%g", f);
            //                DrawString(x, ybot + 1, s, 0.5, 1);
            //                if (i == ntx)
            //                {
            //                    int l;
            //                    l = strlen(s);
            //                    sprintf(s, "                MHz" + (15 - (l + 1) / 2));
            //                    DrawString(x, ybot + 1, s, 0, 1);
            //                }
            //            }
            //            if (idx1 < 0) { min1 = 1; max1 = 2; }
            //            /* horizontal grid lines and associated labels */
            //            for (i = 0; i <= nty; i++)
            //            {
            //                double f;
            //                int y;
            //                char s[20];
            //                f = min1 + (max1 - min1) * ((double)i) / nty;
            //                if (fabs(f / (max1 - min1)) < 0.1 / nty) f = 0;
            //                y = sy1(f);
            //                if (i > 0 && i < nty)
            //                {
            //                    SetForeground(&c_scale);
            //                    DrawLine(xleft, y, xright, y);
            //                }
            //                if (idx1 >= 0)
            //                {
            //                    sprintf(s, "%g  ", f);
            //                    SetForeground(color1);
            //                    DrawString(xleft, y, s, 1, 0.5);
            //                }
            //                if (idx2 >= 0)
            //                {
            //                    f = min2 + (max2 - min2) * ((double)i) / nty;
            //                    if (fabs(f / (max2 - min2)) < 0.1 / nty) f = 0;
            //                    y = sy2(f);
            //                    sprintf(s, "  %g", f);
            //                    SetForeground(color2);
            //                    DrawString(xright, y, s, 0, 0.5);
            //                }
            //            }
            //            SetForeground(&c_axis);
            //            /* border around the graph */
            //            DrawLine(xleft, ybot, xright, ybot);
            //            DrawLine(xleft, ytop, xright, ytop);
            //            DrawLine(xleft, ybot, xleft, ytop);
            //            DrawLine(xright, ybot, xright, ytop);

            //            /* title(s) */
            //            if (title)
            //            {
            //                SetForeground(&c_axis);
            //                DrawString((xleft + xright) / 2, ytop - 1, title, 0.5, 0);
            //            }
            //            if (title1)
            //            {
            //                SetForeground(color1);
            //                DrawString(xleft - fontheight, ytop - 1, title1, 0, 0);
            //            }
            //            if (title2)
            //            {
            //                SetForeground(color2);
            //                DrawString(xright + fontheight, ytop - 1, title2, 1, 0);
            //            }

            //            /* the actual data points and connecting lines */
            //            SetClipRectangle(xleft - 2, ytop, xright + 2, ybot);
            //            xx1 = xx2 = yy1 = yy2 = -1;
            //            xx1a = xx2a = yy1a = yy2a = -1;
            //            for (i = 0, ne = neco; i < numneco; i++, ne++)
            //            {
            //                int x, y;
            //                x = sx(ne->f);
            //                if (idx1a >= 0 || idx2a >= 0) SetLineAttributes(0, GDK_LINE_ON_OFF_DASH, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                if (idxOK(idx1a, ne))
            //                {
            //                    y = sy1(ne->d[idx1a]);
            //                    SetForeground(color1);
            //                    if (numneco < win2sizex / 4)
            //                    {
            //                        DrawLine(x - 2, y - 2, x + 2, y - 2);
            //                        DrawLine(x - 2, y + 2, x + 2, y + 2);
            //                        DrawLine(x - 2, y - 2, x - 2, y + 2);
            //                        DrawLine(x + 2, y - 2, x + 2, y + 2);
            //                    }
            //                    if (xx1a != -1) DrawLine(xx1a, yy1a, x, y);
            //                    xx1a = x; yy1a = y;
            //                }
            //                if (idxOK(idx2a, ne))
            //                {
            //                    y = sy2(ne->d[idx2a]);
            //                    SetForeground(color2);
            //                    if (numneco < win2sizex / 4)
            //                    {
            //                        DrawLine(x - 3, y, x, y - 3);
            //                        DrawLine(x - 3, y, x, y + 3);
            //                        DrawLine(x + 3, y, x, y - 3);
            //                        DrawLine(x + 3, y, x, y + 3);
            //                    }
            //                    if (xx2a != -1) DrawLine(xx2a, yy2a, x, y);
            //                    xx2a = x; yy2a = y;
            //                }
            //                if (idx1a >= 0 || idx2a >= 0) SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                if (idxOK(idx1, ne))
            //                {
            //                    y = sy1(ne->d[idx1]);
            //                    SetForeground(color1);
            //                    if (numneco < win2sizex / 4)
            //                    {
            //                        DrawLine(x - 2, y - 2, x + 2, y - 2);
            //                        DrawLine(x - 2, y + 2, x + 2, y + 2);
            //                        DrawLine(x - 2, y - 2, x - 2, y + 2);
            //                        DrawLine(x + 2, y - 2, x + 2, y + 2);
            //                    }
            //                    if (xx1 != -1) DrawLine(xx1, yy1, x, y);
            //                    xx1 = x; yy1 = y;
            //                }
            //                if (idxOK(idx2, ne))
            //                {
            //                    y = sy2(ne->d[idx2]);
            //                    SetForeground(color2);
            //                    if (numneco < win2sizex / 4)
            //                    {
            //                        DrawLine(x - 3, y, x, y - 3);
            //                        DrawLine(x - 3, y, x, y + 3);
            //                        DrawLine(x + 3, y, x, y - 3);
            //                        DrawLine(x + 3, y, x, y + 3);
            //                    }
            //                    if (xx2 != -1) DrawLine(xx2, yy2, x, y);
            //                    xx2 = x; yy2 = y;
            //                }
            //            }
            //            SetClipRectangle(0, 0, win2sizex, win2sizey);
        }

        double xfreq(int x)
        {
            return ((double)x - xleft) / (xright - xleft) * (maxf - minf) + minf;
        }

        int freqx(double f)
        {
            return 0;   // TODO FIX THIS sx(f);
        }

        private void aboutToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }

        int freqindex(double f)
        {
            double d, dbest;
            int i, ibest;

            ibest = -1;
            dbest = Double.MaxValue;

            //            for (i = 0; i < numneco; i++)
            //            {
            //                if (!neco[i].rp && !neco[i].cu && !neco[i].nf) continue;
            //                d = fabs(f - neco[i].f);
            //                if (d < dbest)
            //                {
            //                    ibest = i;
            //                    dbest = d;
            //                }
            //            }
            return ibest;
        }



        public static void Draw_all2(int onlyvgain)
        {
            int n;
            double size, step, y;

            //            n = plot2_swr + plot2_z + plot2_z2 + plot2_maxgain + plot2_dir + plot2_vgain;
            //            size = 1.0 / (n + (n - 1) * 0.05);
            //            step = 1.05 * size;
            //            y = 1.0;

            //            if (!onlyvgain) ClearWindow();
            //            if (plot2_z2)
            //            {
            //                if (!onlyvgain) freqplot(neco_zphi, neco_zabs, -1, -1, "phi(Z) [deg]", "[ohm] |Z|", "impedance", &c_exci, &c_load, y, y - size);
            //                y -= step;
            //            }
            //            if (plot2_z)
            //            {
            //                if (!onlyvgain) freqplot(neco_zr, neco_zi, -1, -1, "real [ohm]", "[ohm] imag", "impedance", &c_exci, &c_load, y, y - size);
            //                y -= step;
            //            }
            //            if (plot2_swr)
            //            {
            //                if (!onlyvgain) freqplot(neco_swr, -1, -1, -1, "SWR", NULL, NULL, &c_wire, NULL, y, y - size);
            //                y -= step;
            //            }
            //            if (plot2_dir)
            //            {
            //                if (!onlyvgain)
            //                {
            //                    if (polarization == POLnone || polarization == POLcolour) freqplot(neco_phi, neco_theta, -1, -1, "phi [deg]", "[deg] theta", "direction of maximum gain", &c_exci, &c_load, y, y - size);
            //                    else freqplot(Neco_polphi + Neco_gsize * polarization,
            //                                  Neco_poltheta + Neco_gsize * polarization,
            //                                  neco_phi, neco_theta, "phi", "theta", "direction of maximum gain", &c_exci, &c_load, y, y - size);
            //                }
            //                y -= step;
            //            }
            //            if (plot2_maxgain)
            //            {
            //                if (!onlyvgain)
            //                {
            //                    if (polarization == POLnone || polarization == POLcolour) freqplot(neco_maxgain, neco_fb, -1, -1, "gain [dB]", "[dB] f/b", "in direction of maximum gain", &c_gain, &c_surf, y, y - size);
            //                    else freqplot(Neco_polgain + Neco_gsize * polarization,
            //                                  Neco_polfb1 + Neco_gsize * polarization,
            //                                  neco_maxgain,
            //                                  Neco_polfb2 + Neco_gsize * polarization,
            //                                  "gain", "f/b", "in direction of maximum gain", &c_gain, &c_surf, y, y - size);
            //                }
            //                y -= step;
            //            }
            //            if (plot2_vgain)
            //            {
            //                if (onlyvgain) ClearRectangle(0, (y - size) * win2sizey, win2sizex, y * win2sizey);
            //                if (polarization == POLnone || polarization == POLcolour) freqplot(neco_vgain, neco_vfb, -1, -1, "gain [dB]", "[dB] f/b", "in direction toward viewer", &c_gain, &c_surf, y, y - size);
            //                else freqplot(neco_vgain2, neco_vfb, neco_vgain, neco_vfb2, "gain [dB]", "[dB] f/b", "in direction toward viewer", &c_gain, &c_surf, y, y - size);
            //                y -= step;
            //            }
            //   out->Complete();
        }

        private void exitToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

    }
}
