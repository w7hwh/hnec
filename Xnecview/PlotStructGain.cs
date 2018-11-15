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
    /// 
    /// </summary>
    public partial class PlotStructGain : Form
    {
        /// <summary>
        /// 
        /// </summary>
        public PlotStructGain()
        {
            InitializeComponent();
        }

        //    /* ----------- initialize window 1 (struct/gain plot) ------------------------------------------------------------- */

        //    if (window1open)
        //    {

        //        acc = gtk_accel_group_new();

        //        /* ----- create the top row of buttons ----- */

        //        toprow = gtk_hbox_new(FALSE, 0);

        //        w = gtk_button_new_with_label("quit"); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        gtk_signal_connect(GTK_OBJECT(w), "clicked", GTK_SIGNAL_FUNC(cmd_quit), (gpointer)NULL);
        //        gtk_widget_add_accelerator(w, "clicked", acc, GDK_Q, 0, 0);
        //        SET_TT("Quit the program");

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
        //            gtk_signal_connect(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(cmd_export_eps), (gpointer)cmd_export_eps_ok);
        //# ifdef HAVE_LIBPNG
        //            w = gtk_menu_item_new_with_label("PNG"); gtk_menu_append(GTK_MENU(m), w); gtk_widget_show(w);
        //            gtk_signal_connect(GTK_OBJECT(w), "activate", GTK_SIGNAL_FUNC(cmd_export_png), (gpointer)cmd_export_png_ok);
        //#endif
        //            w = gtk_menu_bar_new(); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //            gtk_menu_bar_append(GTK_MENU_BAR(w), v);
        //        }

        //        w = make_option_menu(SPnames, cmd_setsp, structplot); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        SET_TT("Choose type of antenna structure display");

        //        w = make_option_menu(GPnames, cmd_setgp, gainplot); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        SET_TT("Choose type of radiation pattern display");

        //        w = make_option_menu(GSnames, cmd_setgs, gainscale); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        SET_TT("Choose gain scale");

        //        w = make_option_menu(POLnames, cmd_setpol, polarization); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        SET_TT("Choose polarisation");

        //        w = gtk_button_new_with_label("X"); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        gtk_signal_connect(GTK_OBJECT(w), "clicked", GTK_SIGNAL_FUNC(cmd_X), (gpointer)NULL);
        //        SET_TT("View along X axis");

        //        w = gtk_button_new_with_label("Y"); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        gtk_signal_connect(GTK_OBJECT(w), "clicked", GTK_SIGNAL_FUNC(cmd_Y), (gpointer)NULL);
        //        SET_TT("View along Y axis");

        //        w = gtk_button_new_with_label("Z"); gtk_box_pack_start(GTK_BOX(toprow), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        gtk_signal_connect(GTK_OBJECT(w), "clicked", GTK_SIGNAL_FUNC(cmd_Z), (gpointer)NULL);
        //        SET_TT("View along Z axis");

        //        w = gtk_label_new(""); gtk_box_pack_start(GTK_BOX(toprow), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        msgwidget = GTK_LABEL(w);

        //        /* ----- create the bottom row of buttons for the currents display ----- */

        //        botrow_curr = gtk_hbox_new(FALSE, 0);

        //        w = make_option_menu(DCnames, cmd_setdc, distcor); gtk_box_pack_start(GTK_BOX(botrow_curr), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        SET_TT("Toggle distance compensation");

        //        w = gtk_toggle_button_new_with_label("lock"); gtk_box_pack_start(GTK_BOX(botrow_curr), w, FALSE, FALSE, 0); gtk_widget_show(w);
        //        gtk_signal_connect(GTK_OBJECT(w), "toggled", GTK_SIGNAL_FUNC(cmd_togglelock), (gpointer) & phasedirlock);
        //        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), phasedirlock);
        //        SET_TT("(Un)Lock direction for phase");

        //        adj = gtk_adjustment_new(phaseoffset, -180.0, 180.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_phaseoffset), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_curr), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set phase offset");

        //        adj = gtk_adjustment_new(maxthickness, 0.0, 50.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_maxthickness), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_curr), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set thickness");

        //        botrow_curr->requisition.width = 10;  /* make sure that the window can later be resized to something "too" small, if the user wishes so */

        //        /* ----- create the bottom row of buttons for the animated displays of currents, charges and near field ----- */

        //        botrow_anim = gtk_hbox_new(FALSE, 0);

        //        adj = gtk_adjustment_new(iscale, 0.0, 10.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_iscale), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_anim), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set scale factor for currents");

        //        adj = gtk_adjustment_new(qscale, 0.0, 25.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_qscale), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_anim), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set scale factor for charges");

        //        adj = gtk_adjustment_new(escale, 0.0, 25.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_escale), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_anim), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set electric field scale factor");

        //        adj = gtk_adjustment_new(hscale, 0.0, 25.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_hscale), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_anim), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set magnetic field scale factor");

        //#define ADDTOGGLE1anim(label,p,s)  \
        //        w = gtk_toggle_button_new_with_label(label); gtk_box_pack_start(GTK_BOX(botrow_anim), w, FALSE, FALSE, 0); gtk_widget_show(w);  \
        //   gtk_signal_connect(GTK_OBJECT(w), "toggled", GTK_SIGNAL_FUNC(toggle1helper), (gpointer)p);  \
        //   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), *p);  \
        //   SET_TT(s);

        //        gdraw1 = w;  /* necessary since in the following toggle1helper() will be called, so gdraw1 should be valid */
        //        ADDTOGGLE1anim("P", &show_p, "Toggle display of Poynting vector");

        //        adj = gtk_adjustment_new(animfreq, 0.0, 2.0, 0.0, 0.0, 0.0);
        //        gtk_signal_connect(adj, "value_changed", GTK_SIGNAL_FUNC(cmd_animfreq), (gpointer)NULL);
        //        w = gtk_hscale_new(GTK_ADJUSTMENT(adj)); gtk_box_pack_start(GTK_BOX(botrow_anim), w, TRUE, TRUE, 0); gtk_widget_show(w);
        //        SET_TT("Set animation frequency");



        //        /* ----- create vbox containing everything ----- */

        //        vbox1 = gtk_vbox_new(FALSE, 0);
        //        gtk_box_pack_start(GTK_BOX(vbox1), toprow, FALSE, FALSE, 0);
        //        gtk_box_pack_end(GTK_BOX(vbox1), botrow_curr, FALSE, FALSE, 0);
        //        gtk_box_pack_end(GTK_BOX(vbox1), botrow_anim, FALSE, FALSE, 0);
        //        gdraw1 = w = gtk_drawing_area_new();
        //        gtk_drawing_area_size(GTK_DRAWING_AREA(w), winsizex, winsizey);
        //        gtk_box_pack_start(GTK_BOX(vbox1), w, TRUE, TRUE, 0);
        //        gtk_widget_show(w);

        //        gtk_signal_connect(GTK_OBJECT(w), "configure_event", GTK_SIGNAL_FUNC(resize_event), NULL);
        //        gtk_signal_connect(GTK_OBJECT(w), "button_press_event", GTK_SIGNAL_FUNC(buttonpress_event), NULL);
        //        gtk_signal_connect(GTK_OBJECT(w), "button_release_event", GTK_SIGNAL_FUNC(buttonrelease_event), NULL);
        //        gtk_signal_connect(GTK_OBJECT(w), "motion_notify_event", GTK_SIGNAL_FUNC(motion_event), NULL);
        //        gtk_signal_connect_after(GTK_OBJECT(w), "key_press_event", GTK_SIGNAL_FUNC(keypress_event), NULL);
        //        gtk_widget_set_events(w, GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK |
        //           GDK_BUTTON_RELEASE_MASK | GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK |
        //           GDK_KEY_PRESS_MASK);
        //        GTK_WIDGET_SET_FLAGS(w, GTK_CAN_FOCUS);

        //        gtk_container_add(GTK_CONTAINER(win1), vbox1);
        //        gtk_widget_show(toprow);
        //        if (structplot == SPcurrents) gtk_widget_show(botrow_curr);
        //        if (gainplot == GPnearfield || structplot == SPanim) gtk_widget_show(botrow_anim);
        //        gtk_widget_show(vbox1);
        //        if (really) gtk_widget_show(win1);
        //        else gtk_widget_realize(gdraw1);

        //        gtk_widget_grab_focus(gdraw1);

        //        gtk_window_add_accel_group(GTK_WINDOW(win1), acc);

        //        ggc = gdk_gc_new(w->window);

        //        if (toprow->requisition.width > winsizex)
        //        {
        //            winsizex = toprow->requisition.width;
        //            calcproj();
        //        }
        //        gbackg = gdk_pixmap_new(w->window, winsizex, winsizey, gdk_drawable_get_depth(w->window));

        //        upd_msg();

        //        icon = gdk_bitmap_create_from_data(win1->window, icon_bits, icon_width, icon_height);
        //        gdk_window_set_icon(win1->window, NULL, icon, NULL);
        //    }
    }
}
