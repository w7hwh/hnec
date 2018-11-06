/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 1998-2002, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    public static class Config
    {
        /* some defines: */
        public static readonly string VERSION = "1.36";           //TODO pull this from assembly info

        public static readonly int ini_phi = 30;                  /* initial rotation around Z-axis */
        public static readonly int ini_theta = 70;                /* initial tilting of Z-axis */
        public static readonly double ini_zoom = 1.0;             /* initial zoom-factor */
        public static readonly int ini_trx = 0;                   /* initial horizontal displacement of origin w.r.t. center of window */
        public static readonly int ini_try = 0;                   /* initial vertical displacement of origin w.r.t. center of window */

        public static readonly int ini_winsize = 500;             /* initial size of window (pixels) */

        public static readonly int MAXWIRES = 2000;               /* maximum number of wires */
        public static readonly int MAXSURFACES = 1000;            /* maximum number of surfaces */
        public static readonly int MAXEXCIS = 100;                /* maximum number of excitations */
        public static readonly int MAXLOADS = 1000;               /* maximum number of loads */
        public static readonly int MAXNETWS = 30;                 /* maximum number of "networks" (including transmission lines) */

        public static readonly int MAXFREQS = 100;                /* maximum number of frequencies in the NEC output file */
        public static readonly int MAXSEGMENTS = 600;             /* maximum number of segments in the currents data in the NEC output file */
        public static readonly int MAXPATCHES = 600;              /* maximum number of surface patches in the NEC output file */

        public static readonly double Axislen = 1.1;              /* length of axes, compared to largest dimension of the antenna or gain grid */
        public static readonly double GAINSIZE = 3.0;             /* ratio between size of gain grid and largest dimension of the antenna */

        public static readonly string C_BG = "white";             /* color of the background */
        public static readonly string C_AXIS = "black";           /* color of the axes */
        public static readonly string C_WIRE = "blue";            /* color of the wires */
        public static readonly string C_SURF = "green3";          /* color of the front of surfaces */
        public static readonly string C_BACK = "cyan";            /* color of the back of surfaces */
        public static readonly string C_GAIN = "red";             /* color of the gain pattern */
        public static readonly string C_SCALE = "gray70";         /* color of the gain scale */
        public static readonly string C_EXCI = "orange";          /* color of wire segment containing excitation */
        public static readonly string C_LOAD = "brown";           /* color of loaded wire segment */
        public static readonly string C_NETW = "grey";            /* color of networks and transmission lines */
        public static readonly string C_INACTIVE = "gray86";      /* color in currents display of segments in which negligibly little current flows */
        public static readonly string C_EFIELD = "red";           /* color of electric field vectors */
        public static readonly string C_HFIELD = "green2";        /* color of magnetic field vectors */
        public static readonly string C_POYNTING = "yellow2";     /* color of Poynting vectors */
        public static readonly string C_QPOS = "cyan3";           /* color of positive charge */
        public static readonly string C_QNEG = "magenta";         /* color of negative charge */
        public static readonly string C_GAIN_LIN = "red";         /* color of gain when mostly linearly polarized (i.e., axial ratio < Polthr ) */
        public static readonly string C_GAIN_LHCP = "violet red"; /* color of gain when mostly left-hand circularly polarized */
        public static readonly string C_GAIN_RHCP = "orange";     /* color of gain when mostly right-hand circularly polarized */

        public static readonly int NC_phase = 64;                 /* number of different hues for the phase scale; must be a multiple of 4 */
        public static readonly int Default_anim_update_interval = 100; /* time (in milliseconds) between updates of the animation display */

        public static readonly string XFONT = "6x10";             /* font for text in the on-screen drawing */
        public static readonly string PSFONT = "helvetica";       /* font for postscript output (size is derived by scaling the X font) */

        public static readonly double R0 = 50.0;                  /* default reference impedance for SWR calculation */

        public static readonly double Polthr = Math.Sqrt(2) - 1;  /* threshold of axial ratio used in polarization-colouring */
    }
}
