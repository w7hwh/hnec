/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 1998-2006, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 *
 *  This module contains almost all 3D drawing related stuff, like coordinate
 *  transformation and drawing routines based on generic low-level graphics
 *  routines, as provided by the Outdev-struct pointed to by *out.
 *  The only 3D drawing code not included here, is the code for drawing the
 *  the gain pattern as an opaque surface, i.e., with hidden-line removal;
 *  that code is in draw_opaque.c, because of its size and complexity.
 *
 *  Note: the code for drawing the 2D plots of several quantities versus
 *  frequency is in 'freqplot.c'.
 *  
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    class Draw
    {
        /* ------------------------- coordinate conversion and related arithmetic -------------------------------------------------- */

        /* these define the projection: */
        double phi = Config.ini_phi;        /* angle defining direction of eye */
        double theta = Config.ini_theta;      /* angle defining direction of eye */
        double zoom = Config.ini_zoom;       /* zoom factor */
        double trx = Config.ini_trx;
        double Try = Config.ini_try;        /* 2D-translation, as a fraction of winsize */
        int winsizex, winsizey;               /* size of window in pixels */

        /* and these are calculated from them: */
        double Xx, Xy, Xz, Yx, Yy, Yz;       /* projection matrix */
        double Xtr, Ytr;                     /* 2D-translation */
        int winsize;                         /* min(winsizex,winsizey) */
        Nec.Point eyevec;               /* tip of arrow pointing toward observer */
        double Zx, Zy, Zz;                   /* weighting factors for calculating distance to observer (possibly in "locked" direction) */

        void Calcproj()     /* calculate the projection matrix etc. */
        {
            double ph, th;
            double scale;

            winsize = winsizey < winsizex ? winsizey : winsizex;
            scale = winsize * zoom / (3.0 * Nec.extr);
            ph = phi * (Math.PI / 180.0);
            th = theta * (Math.PI / 180.0);
            Xz = 0;
            Yz = -scale * Math.Sin(th);
            Xx = -scale * Math.Sin(ph);
            Xy = scale * Math.Cos(ph);
            Yx = scale * Math.Cos(ph) * Math.Cos(th);
            Yy = scale * Math.Sin(ph) * Math.Cos(th);
            Xtr = 0.5 + trx * zoom * winsize + 0.5 * winsizex;
            Ytr = 0.5 + Try * zoom * winsize + 0.5 * winsizey;
            eyevec.x = (float)(Math.Sin(th) * Math.Cos(ph));
            eyevec.y = (float)(Math.Sin(th) * Math.Sin(ph));
            eyevec.z = (float)Math.Cos(th);
            if (Nec.phasedirlock == 1)
            {
                ph = Nec.phasephi * (Math.PI / 180.0);
                th = Nec.phasetheta * (Math.PI / 180.0);
            }
            Zz = Math.Cos(th);
            Zx = Math.Sin(th) * Math.Cos(ph);
            Zy = Math.Sin(th) * Math.Sin(ph);
        }

        //void proj(Point* p, double* c)      /* perform the projection: *p defines a point in 3D-space, *c defines a point in the window */
        //            {
        //                c[0] = Xtr + Xx * p->x + Xy * p->y + Xz * p->z;
        //                c[1] = Ytr + Yx * p->x + Yy * p->y + Yz * p->z;
        //            }



        //            void interpolate(double* p, double f, double* a, double* b)   /* interpolate linearly between two points in the window (f=0..1 -> *p=*a..*b) */
        //            {
        //                p[0] = (1 - f) * a[0] + f * b[0];
        //                p[1] = (1 - f) * a[1] + f * b[1];
        //            }



        //            /* test whether a line through two points is "ascending" ( / ) as opposed to "descending" ( \ ) */
        //#define ascending(c1,c2) (((c2)[0]-(c1)[0])*((c2)[1]-(c1)[1])<0)

        //            /* ------------------------- antenna structure drawing ------------------------------------------------------------------- */

        void draw_wires()                /* draw the wires of the antenna */
        {
            //                Wire* wi;
            //                double c0[2], c1[2];
            //                char s[20];

            //                if (numwires == 0) return;

            //                SetForeground(&c_wire);
            //                SetLineAttributes(2, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                for (wi = wires; wi < wires + numwires; wi++)
            //                {
            //                    proj(&wi->p0, c0);
            //                    proj(&wi->p1, c1);
            //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //                    if (structplot == SPtags && !quick && wi->itg > 0 && wi->dispnr)
            //                    {
            //                        sprintf(s, "%i", wi->itg);
            //                        if (ascending(c0, c1))
            //                            DrawStringUL((c0[0] + c1[0]) / 2 + 2, (c0[1] + c1[1]) / 2, s);
            //                        else
            //                            DrawStringLL((c0[0] + c1[0]) / 2 + 2, (c0[1] + c1[1]) / 2, s);
            //                    }
            //                }
        }

        void draw_loads()                 /* (re)draw loaded wires */
        {
            //                Load* lo;
            //                double c0[2], c1[2], c2[2], c3[2];

            //                if (numloads == 0) return;
            //                SetForeground(&c_load);
            //                SetLineAttributes(3, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                for (lo = loads; lo < loads + numloads; lo++)
            //                {
            //                    proj(&lo->wire->p0, c0);
            //                    proj(&lo->wire->p1, c1);
            //                    if (lo->wire->ns > 0)
            //                    {
            //                        interpolate(c2, (double)(lo->firstseg - 1) / (lo->wire->ns), c0, c1);
            //                        interpolate(c3, (double)lo->lastseg / (lo->wire->ns), c0, c1);
            //                        DrawLine(c2[0], c2[1], c3[0], c3[1]);
            //                    }
            //                    else
            //                    {
            //                        DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //                    }
            //                }
        }

        void proj_segment(ref Nec.Wire wire, int seg, ref double c0, ref double c1)
        {
            if (wire.ns > 0)
            {
                //                    double d0[2], d1[2];
                //                    proj(&wire->p0, d0);
                //                    proj(&wire->p1, d1);
                //                    interpolate(c0, (double)(seg - 1) / (wire->ns), d0, d1);
                //                    interpolate(c1, (double)seg / (wire->ns), d0, d1);
            }
            else
            {
                //                    proj(&wire->p0, c0);
                //                    proj(&wire->p1, c1);
            }
        }

        void draw_excis()                 /* (re)draw wires with excitations */
        {
            //                Exci* ex;
            //                double c0[2], c1[2];

            //                if (numexcis == 0) return;
            //                SetForeground(&c_exci);
            //                SetLineAttributes(4, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                for (ex = excis; ex < excis + numexcis; ex++)
            //                {
            //                    proj_segment(ex->wire, ex->seg, c0, c1);
            //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //                }
        }

        void draw_netws()                 /* (re)draw networks and transmission lines */
        {
            //                Netw* ne;
            //                double c10[2], c11[2], c20[2], c21[2];  /* endpoints of relevant wires */

            //                if (Nec.netws.Count() == 0) return;
            //                SetForeground(&c_netw);
            //                for (ne = netws; ne < netws + Nec.netws.Count(); ne++)
            //                {
            //                    proj_segment(ne->wire1, ne->seg1, c10, c11);
            //                    proj_segment(ne->wire2, ne->seg2, c20, c21);
            //                    if (ne->type == 1)
            //                    {
            //                        SetLineAttributes(1, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                        DrawLine(c10[0], c10[1], c20[0], c20[1]);
            //                        DrawLine(c11[0], c11[1], c21[0], c21[1]);
            //                    }
            //                    else if (ne->type == -1)
            //                    {
            //                        SetLineAttributes(1, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                        DrawLine(c10[0], c10[1], c21[0], c21[1]);
            //                        DrawLine(c11[0], c11[1], c20[0], c20[1]);
            //                    }
            //                    else
            //                    {
            //                        double d1[2], d2[2];
            //                        interpolate(d1, 0.5, c10, c11);
            //                        interpolate(d2, 0.5, c20, c21);
            //                        SetLineAttributes(2, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                        DrawLine(d1[0], d1[1], d2[0], d2[1]);
            //                    }
            //                }
        }

        void draw_surfaces()
        {
            //                Surface* su;
            //                double c0[2], c1[2], c2[2], c3[2];
            //                double cc[2];

            //                if (numsurfaces == 0) return;
            //                SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                for (su = surfaces; su < surfaces + numsurfaces; su++)
            //                {
            //                    if (eyevec.x * (su->pa.x - su->pc.x) + eyevec.y * (su->pa.y - su->pc.y) + eyevec.z * (su->pa.z - su->pc.z) >= 0)
            //                        SetForeground(&c_surf);
            //                    else
            //                        SetForeground(&c_back);
            //                    proj(&su->pc, cc);
            //                    proj(&su->p0, c0);
            //                    proj(&su->p1, c1);
            //                    switch (su->type)
            //                    {
            //                        case SU_rect:
            //                        case SU_arb:
            //                            c2[0] = 2 * cc[0] - c0[0]; c2[1] = 2 * cc[1] - c0[1];
            //                            c3[0] = 2 * cc[0] - c1[0]; c3[1] = 2 * cc[1] - c1[1];
            //                            break;
            //                        case SU_quad:
            //                            proj(&su->p3, c3);
            //                        /* fallthrough */
            //                        case SU_tri:
            //                            proj(&su->p2, c2);
            //                            break;
            //                    }
            //                    if (su->type == SU_tri)
            //                    {
            //                        /* triangle; draw border: */
            //                        DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //                        DrawLine(c1[0], c1[1], c2[0], c2[1]);
            //                        DrawLine(c2[0], c2[1], c0[0], c0[1]);
            //                    }
            //                    else
            //                    {
            //                        /* joint drawing code for all other shapes, which are drawn as various quadrilaterals */
            //                        /* draw diagonals: */
            //                        DrawLine(c0[0], c0[1], c2[0], c2[1]);
            //                        DrawLine(c1[0], c1[1], c3[0], c3[1]);
            //                        if (su->type != SU_arb)
            //                        {
            //                            /* draw border: */
            //                            DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //                            DrawLine(c1[0], c1[1], c2[0], c2[1]);
            //                            DrawLine(c2[0], c2[1], c3[0], c3[1]);
            //                            DrawLine(c3[0], c3[1], c0[0], c0[1]);
            //                        }
            //                        else
            //                        {
            //                            /* draw additional diagonals */
            //                            DrawLine(0.7071 * (c0[0] + c1[0]) - 0.4142 * cc[0], 0.7071 * (c0[1] + c1[1]) - 0.4142 * cc[1], 0.7071 * (c2[0] + c3[0]) - 0.4142 * cc[0], 0.7071 * (c2[1] + c3[1]) - 0.4142 * cc[1]);
            //                            DrawLine(0.7071 * (c0[0] + c3[0]) - 0.4142 * cc[0], 0.7071 * (c0[1] + c3[1]) - 0.4142 * cc[1], 0.7071 * (c2[0] + c1[0]) - 0.4142 * cc[0], 0.7071 * (c2[1] + c1[1]) - 0.4142 * cc[1]);
            //                        }
            //                    }
            //                }
        }

        void draw_antenna(int ie)              /* draw the entire antenna */
        {
            draw_wires();
            draw_loads();
            draw_excis();
            draw_netws();
            if ((Xwin.Pending() > 0) && (ie > 0))
                return;
            draw_surfaces();
        }


        //            /* ------------------------- currents/charges drawing -------------------------------------------------------------------- */

        //            int cmpsegdist(const Segcur** s1, const Segcur** s2)
        //{
        //                double d;
        //                d = (*s1)->dist - (*s2)->dist;
        //                if (d < 0) return -1;
        //                if (d > 0) return 1;
        //                return 0;
        //            }


        //            void addvector(Point* pe, float a[3],double q)
        //{
        //                pe->x += q * a[0];
        //                pe->y += q * a[1];
        //                pe->z += q * a[2];
        //            }

        //            void draw_currents_anim(int ie, Currents* cu)      /* draw the currents and charge distribution as vectors and circles with animation */
        //            {
        //                Segcur* se;
        //                double si, co;
        //                double sc;
        //                double l;

        //                if (cu == NULL || cu->numseg == 0) return;
        //                l = 0;
        //                l += (cu->s[0].p1.x - cu->s[0].p0.x) * (cu->s[0].p1.x - cu->s[0].p0.x);
        //                l += (cu->s[0].p1.y - cu->s[0].p0.y) * (cu->s[0].p1.y - cu->s[0].p0.y);
        //                l += (cu->s[0].p1.z - cu->s[0].p0.z) * (cu->s[0].p1.z - cu->s[0].p0.z);
        //                l = sqrt(l);

        //                si = sin(-animphase);
        //                co = cos(-animphase);
        //                sc = iscale * 0.8 * l / cu->maxI;

        //                for (se = cu->s; se < cu->s + cu->numanimseg; se++)
        //                {
        //                    double c0[2], c1[2];
        //                    double q;
        //                    Point pe;
        //                    double ev[3];
        //                    int i;

        //                    for (i = 0; i < 3; i++) ev[i] = se->re[i] * co + se->im[i] * si;
        //                    q = qscale * l * 0.5 * (se->qre * co + se->qim * si) / cu->maxQ;

        //                    pe.x = se->c.x + sc * ev[0]; pe.y = se->c.y + sc * ev[1]; pe.z = se->c.z + sc * ev[2];
        //                    proj(&se->c, c0);
        //                    proj(&pe, c1);
        //                    SetLineAttributes(2, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
        //                    SetForeground(&c_wire);
        //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);

        //                    SetLineAttributes(1, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
        //                    if (q < 0) SetForeground(&c_qneg); else SetForeground(&c_qpos);
        //                    q = fabs(q);
        //                    pe = se->c;
        //                    addvector(&pe, se->q0, q); proj(&pe, c0);
        //                    addvector(&pe, se->q0, -q); addvector(&pe, se->q1, q); proj(&pe, c1); DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                    addvector(&pe, se->q0, -q); addvector(&pe, se->q1, -q); proj(&pe, c0); DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                    addvector(&pe, se->q0, q); addvector(&pe, se->q1, -q); proj(&pe, c1); DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                    addvector(&pe, se->q0, q); addvector(&pe, se->q1, q); proj(&pe, c0); DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                }
        //            }


        //            void draw_currents_noanim(int ie, Currents* cu)      /* draw the currents distribution as a static phase/magnitude display */
        //            {
        //                Segcur* se;
        //                Segcur** sse;
        //                int i;
        //                double c0[2], c1[2];
        //                double a = 0, phase;
        //                double dx, dy, l;
        //                double w;
        //# ifdef GAINCALC
        //                double totalr = 0, totali = 0;
        //#endif

        //                if (cu == NULL || cu->numseg == 0) return;
        //                sse = malloc(cu->numseg * sizeof(Segcur*));
        //                if (sse == NULL) return;

        //                /* first, calculate the distance of each segment (defined in xnecview.h),
        //                   and at the same time fill the array **sse that will be sorted next */
        //                se = cu->s;
        //                for (i = 0; i < cu->numseg; i++)
        //                {
        //                    se->dist = Zx * se->c.x + Zy * se->c.y + Zz * se->c.z;
        //                    sse[i] = se;
        //                    se++;
        //                }

        //                /* next, sort the array elements according to this distance, so the element closest to the observer will be drawn last */
        //                qsort(sse, cu->numseg, sizeof(Segcur*),
        //                   (int(*)(const void*, const void*))cmpsegdist);

        //            w = 360 / (299.8 / neco[rp_index].f);
        //            for (i = 0; i < cu->numseg; i++)
        //            {
        //                se = sse[i];
        //                proj(&se->p0, c0);
        //                proj(&se->p1, c1);

        //                dx = c1[0] - c0[0];
        //                dy = c1[1] - c0[1];
        //                if (!quick)
        //                {
        //                    l = sqrt(dx * dx + dy * dy);
        //                    if (l == 0) continue;
        //                    phase = se->phase;
        //                    switch (polarization)
        //                    {
        //                        case POLnone:
        //                        case POLcolour:
        //                            a = se->a;
        //                            break;
        //                        case POLhor:
        //                            a = se->a * dx / l;
        //                            break;
        //                        case POLvert:
        //                            a = se->a * dy / l;
        //                            break;
        //                        case POLlhcp:
        //                            a = se->a;
        //                            phase += 180 / M_PI * atan2(dy, dx);
        //                            break;
        //                        case POLrhcp:
        //                            a = se->a;
        //                            phase -= 180 / M_PI * atan2(dy, dx);
        //                            break;
        //                    }
        //                    if (a < 0) { phase += 180; a = -a; }
        //                    if (distcor) phase += w * se->dist;
        //                    phase += phaseoffset;
        //                    while (phase < 0) phase += 360;
        //                    while (phase >= 360) phase -= 360;
        //                }

        //                if (!quick && maxthickness * a >= 0.5)
        //                {
        //                    SetLineAttributes(maxthickness * a + 0.5, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
        //                    SetForeground(c_currents + (int)(phase * NC_phase / 360));
        //                }
        //                else
        //                {
        //                    SetLineAttributes(1, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
        //                    SetForeground(&c_inactive);
        //                }
        //                DrawLine(c0[0], c0[1], c1[0], c1[1]);

        //# ifdef GAINCALC
        //                totalr += a * l * cos(phase * M_PI / 180);
        //                totali += a * l * sin(phase * M_PI / 180);
        //#endif
        //            }
        //# ifdef GAINCALC
        //            totalr = sqrt(totalr * totalr + totali * totali);
        //            printf("-->  %12g %12g\n", totalr, 20 * log10(totalr));
        //#endif

        //            free(sse);
        //        }

        private void Draw_phasecircle()
        {
            int xc, yc;  /* coordinates of center */
            int r;      /* radius */
            int i, j, x, y;
            double phase, step;
            double rc, rs;
            //            double sintab[NC_phase / 4];
            //            double costab[NC_phase / 4];

            //            r = winsize / 10;
            //            xc = yc = r + r / 2 + 1;
            //            SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //            step = 1 / (2 * M_PI * r) * 360;

            //            for (i = 0; i < NC_phase / 4; i++)
            //            {
            //                phase = ((i * 360.0 / NC_phase) - phaseoffset) * M_PI / 180.0;
            //                sintab[i] = sin(phase);
            //                costab[i] = cos(phase);
            //            }

            //            j = 3 * (int)(2 * M_PI * r / NC_phase + 2);
            //            do
            //            {
            //                j /= 3;
            //                SetLineAttributes(j, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //                for (i = 0; i < NC_phase / 4; i++)
            //                {
            //                    rs = r * sintab[i];
            //                    rc = r * costab[i];
            //                    x = xc + rc;
            //                    y = yc + rs;
            //                    SetForeground(c_currents + i);
            //                    DrawLine(xc, yc, x, y);
            //                    x = xc - rs;
            //                    y = yc + rc;
            //                    SetForeground(c_currents + i + NC_phase / 4);
            //                    DrawLine(xc, yc, x, y);
            //                    x = xc - rc;
            //                    y = yc - rs;
            //                    SetForeground(c_currents + i + NC_phase / 2);
            //                    DrawLine(xc, yc, x, y);
            //                    x = xc + rs;
            //                    y = yc - rc;
            //                    SetForeground(c_currents + i + 3 * NC_phase / 4);
            //                    DrawLine(xc, yc, x, y);
            //                }
            //            } while (j > 0);
        }


        //        /* ------------------------- near field drawing -------------------------------------------------------------------------- */

        //        void draw_nearfield(int ie, NECoutput* ne)      /* draw the E and H field and Poynting vectors */
        //        {
        //            double f;
        //            Nearfield* nf;
        //            double sce, sch, scp;
        //            double si, co;

        //            nf = ne->nf;
        //            if (nf == NULL) return;

        //            if (nf->next == NULL)
        //            {
        //                f = 1;
        //            }
        //            else
        //            {
        //                f = (nf->next->p.x - nf->p.x) * (nf->next->p.x - nf->p.x);
        //                f += (nf->next->p.y - nf->p.y) * (nf->next->p.y - nf->p.y);
        //                f += (nf->next->p.z - nf->p.z) * (nf->next->p.z - nf->p.z);
        //                f = sqrt(f);
        //            }
        //            sce = f * escale / ne->maxe;
        //            sch = f * hscale / ne->maxh;
        //            scp = 2 * f * escale * hscale / (ne->maxe * ne->maxh);

        //            SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
        //            si = sin(-animphase);
        //            co = cos(-animphase);
        //            while (nf)
        //            {
        //                Point pe;
        //                double ev[3], hv[3], pv[3];
        //                double c0[2], c1[2];
        //                int valid;

        //                proj(&nf->p, c0);
        //                valid = 1;

        //                if (nf->evalid)
        //                {
        //                    int i;
        //                    for (i = 0; i < 3; i++) ev[i] = nf->ere[i] * co + nf->eim[i] * si;
        //                    pe.x = nf->p.x + sce * ev[0]; pe.y = nf->p.y + sce * ev[1]; pe.z = nf->p.z + sce * ev[2];
        //                    proj(&pe, c1);
        //                    SetForeground(&c_efield);
        //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                }
        //                else valid = 0;

        //                if (nf->hvalid)
        //                {
        //                    int i;
        //                    for (i = 0; i < 3; i++) hv[i] = nf->hre[i] * co + nf->him[i] * si;
        //                    pe.x = nf->p.x + sch * hv[0]; pe.y = nf->p.y + sch * hv[1]; pe.z = nf->p.z + sch * hv[2];
        //                    proj(&pe, c1);
        //                    SetForeground(&c_hfield);
        //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                }
        //                else valid = 0;

        //                if (valid && show_p)
        //                {
        //                    pv[0] = ev[1] * hv[2] - ev[2] * hv[1];
        //                    pv[1] = ev[2] * hv[0] - ev[0] * hv[2];
        //                    pv[2] = ev[0] * hv[1] - ev[1] * hv[0];
        //                    pe.x = nf->p.x + scp * pv[0]; pe.y = nf->p.y + scp * pv[1]; pe.z = nf->p.z + scp * pv[2];
        //                    proj(&pe, c1);
        //                    SetForeground(&c_poynting);
        //                    DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //                }

        //                nf = nf->next;
        //            }
        //        }

        //        /* ------------------------- gain pattern drawing (non-opaque) ----------------------------------------------------------- */

        private void Draw_polref()
        {
            //            SetForeground(&c_gain_lin); DrawString(winsizex - 5, 5 + 0 * fontheight, "linear", 1, 1);
            //            SetForeground(&c_gain_lhcp); DrawString(winsizex - 5, 5 + 1 * fontheight, "left-hand circular", 1, 1);
            //            SetForeground(&c_gain_rhcp); DrawString(winsizex - 5, 5 + 2 * fontheight, "right-hand circular", 1, 1);
        }


        private void Setpolcolour(double axial)
        {
            //            if (axial > Polthr) SetForeground(&c_gain_lhcp);
            //            else if (axial < -Polthr) SetForeground(&c_gain_rhcp);
            //            else SetForeground(&c_gain_lin);
        }


        //        void draw_gain_theta(int i, Radpattern* rp)
        //        {
        //            int j;
        //            double c0[2], c1[2];
        //            double a0, a1;

        //            proj(&rp->gpo[0][i], c1);
        //            a1 = rp->gain[0][i].axial;
        //            for (j = 1; j < rp->numphi; j++)
        //            {
        //                c0[0] = c1[0]; c0[1] = c1[1];
        //                proj(&rp->gpo[j][i], c1);
        //                if (polarization == POLcolour)
        //                {
        //                    a0 = a1;
        //                    a1 = rp->gain[j][i].axial;
        //                    setpolcolour((a0 + a1) / 2);
        //                }
        //                DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //            }
        //            c0[0] = c1[0]; c0[1] = c1[1];
        //            proj(&rp->gpo[0][i], c1);
        //            if (polarization == POLcolour)
        //            {
        //                a0 = a1;
        //                a1 = rp->gain[0][i].axial;
        //                setpolcolour((a0 + a1) / 2);
        //            }
        //            DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //        }

        //        void draw_gain_phi(int j, Radpattern* rp)
        //        {
        //            int i;
        //            double c0[2], c1[2];
        //            double a0, a1;

        //            proj(&rp->gpo[j][0], c1);
        //            a1 = rp->gain[j][0].axial;
        //            for (i = 1; i < rp->numtheta; i++)
        //            {
        //                c0[0] = c1[0]; c0[1] = c1[1];
        //                proj(&rp->gpo[j][i], c1);
        //                if (polarization == POLcolour)
        //                {
        //                    a0 = a1;
        //                    a1 = rp->gain[j][i].axial;
        //                    setpolcolour((a0 + a1) / 2);
        //                }
        //                DrawLine(c0[0], c0[1], c1[0], c1[1]);
        //            }
        //        }


        private static readonly double R90 = (Math.PI / 2);
        private static readonly double R180 = Math.PI;
        private static readonly double R270 = (1.5 * Math.PI);
        private static readonly double R360 = (2 * Math.PI);

        //# include <sys/timeb.h>

        //        void draw_gain(int ie, Radpattern* rp)
        //        {
        //            int i, j;
        //            int gainplot2 = gainplot;

        //            SetForeground(&c_gain);
        //            SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);

        //            if (gainplot == GPopaque && !quick)
        //            {
        //                Radpattern* p;
        //                int ok = 0;
        //                p = rp;
        //                while (p)
        //                {
        //                    if (!opaque_impossible(p))
        //                    {
        //                        draw_opaque(p);
        //                        ok = 1;
        //                    }
        //                    p = p->next;
        //                }
        //                if (ok) return;
        //                fprintf(stderr, "Opaque gain plot not possible; reason: %s.\n", opaque_impossible(rp));
        //                gainplot2 = GPframe;
        //            }

        //            while (rp)
        //            {
        //                if (gainplot2 == GPframe && !quick)
        //                {

        //                    for (i = 0; i < rp->numtheta; i++) draw_gain_theta(i, rp);
        //                    if (Pending() && ie) return;
        //                    for (j = 0; j < rp->numphi; j++) draw_gain_phi(j, rp);

        //                }
        //                else
        //                {

        //                    if (!scaleplot || theta != 90)      /* suppress if scaleplotting and XY-plane is perpendicular to screen */
        //                        for (i = 0; i < rp->numtheta; i++)
        //                            if (fabs(rp->gtheta[i] - R90) < 0.001)     /* points in XY-plane? */
        //                                draw_gain_theta(i, rp);

        //                    for (j = 0; j < rp->numphi; j++)
        //                        if ((
        //                              ((fabs((rp->gphi[j] + R90) - R90) < 0.001 || fabs((rp->gphi[j]) - R180) < 0.001)  /* points in XZ-plane? */
        //                                 && (!scaleplot || (phi != 0 && theta != 0)))  /* suppress if that plane is perpendicular to screen and we're plotting the scale */
        //                           ) || (
        //                              (fabs((rp->gphi[j]) - R90) < 0.001 || fabs((rp->gphi[j]) - R270) < 0.001)       /* points in YZ-plane? */
        //                                 && (!scaleplot || (phi != 90 && theta != 0))  /* suppress if that plane is perpendicular to screen and we're plotting the scale */
        //                           ))
        //                            draw_gain_phi(j, rp);

        //                }
        //                rp = rp->next;
        //            }
        //        }


        private void Draw_gainscale(int ie)
        {
            double a;
            ////            Point po;
            ////            double c0[2], c1[2], orig[2];
            double Npt = 90;
            double da = R360 / Npt;
            //            static double si[Npt], co[Npt];
            int firsttime = 1;
            //            static double radii[numGS][6]=     /* this array determines the radius of all gain scale circles; note: if you change this, also update the labels in GSnumbers[][] and in labels, here below */
            //           {
            //             {    1,     0.794328,    0.501187,    0.251189,    0.1,       0   },   /* lin.power:    0, -1, -3, -6, -10 dB */
            //             {    1,     0.891251,    0.707946,    0.501187,    0.316228,  0   },   /* lin.voltage:  0, -1, -3, -6, -10 dB */
            //             {    1,     0.839624,    0.558406,    0.311817,    0.174121,  0   },   /* arrl:         0, -3, -10, -20, -30 dB */
            //             {    1,     0.75,        0.5,         0.25,        0   }               /* log:          0, -10, -20, -30 dB */
            //           } ;
            //   static char labels[numGS][6] [4]=
            //           {
            //             { "0dB", "-1", "-3", "-6", "-10" },
            //             { "0dB", "-1", "-3", "-6", "-10" },
            //             { "0dB", "-3", "-10", "-20", "-30" },
            //             { "0dB", "-10", "-20", "-30" }
            //           } ;
            int i, j;
            double r;
            //        double* rp;

            //   if (!scaleplot) return;

            //   SetForeground(&c_scale);
            //        SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);

            //   if (firsttime) {
            //      for (i=Npt-1, a=0; i>=0; i--, a+=da) {
            //         si[i]=sin(a);
            //        co[i]=cos(a);
            //    }
            //    firsttime=0;
            //   }

            //   if (phi!=90 && theta!=0) {
            //      /* draw radial lines (yz plane, 10 degrees spacing) */
            //      po.x=0; po.y=0; po.z=0;
            //      proj(&po, orig);
            //      for (i=0;i<360;i=i+10) {
            //         po.x=0;
            //         po.y=extr* cos(i* M_PI/180);
            //po.z=extr* sin(i* M_PI/180);
            //proj(&po, c1);
            //DrawLine(orig[0], orig[1], c1[0], c1[1]);
            //      }
            //      /* draw circles */
            //      rp=radii[gainscale];
            //      j=0;
            //      while ((r=* rp++* extr)) {
            //         po.x=0; po.y=r; po.z=0;
            //         proj(&po, c0);
            //         for (i=0;i<Npt;i++) {
            //            po.y=r* co[i];
            //po.z=r* si[i];
            //proj(&po, c1);
            //DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //c0[0]=c1[0]; c0[1]=c1[1];
            //         }
            //         po.y=0; po.z=r; proj(&po, c1);
            //DrawStringLL(c1[0], c1[1], labels[gainscale][j++]);
            //      }
            //      if (Pending() && ie) return;
            //   }

            //   if (theta!=90) {
            //      /* draw radial lines (xy plane, 10 degrees spacing) */
            //      po.x=0; po.y=0; po.z=0;
            //      proj(&po, orig);
            //      for (i=0;i<360;i=i+10) {
            //         po.z=0;
            //         po.x=extr* cos(i* M_PI/180);
            //po.y=extr* sin(i* M_PI/180);
            //proj(&po, c1);
            //DrawLine(orig[0], orig[1], c1[0], c1[1]);
            //      }
            //      /* draw circles */
            //      rp=radii[gainscale];
            //      j=0;
            //      while ((r=* rp++* extr)) {
            //         po.x=r; po.y=0; po.z=0;
            //         proj(&po, c0);
            //         for (i=0;i<Npt;i++) {
            //            po.x=r* co[i];
            //po.y=r* si[i];
            //proj(&po, c1);
            //DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //c0[0]=c1[0]; c0[1]=c1[1];
            //         }
            //         po.x=-r; po.y=0; proj(&po, c1);
            //DrawStringLL(c1[0], c1[1], labels[gainscale][j++]);
            //      }
            //      if (Pending() && ie) return;
            //   }

            //   if (phi!=0 && theta!=0) {
            //      /* draw radial lines (xz plane, 10 degrees spacing) */
            //      po.x=0; po.y=0; po.z=0;
            //      proj(&po, orig);
            //      for (i=0;i<360;i=i+10) {
            //         po.y=0;
            //         po.x=extr* cos(i* M_PI/180);
            //po.z=extr* sin(i* M_PI/180);
            //proj(&po, c1);
            //DrawLine(orig[0], orig[1], c1[0], c1[1]);
            //      }
            //      /* draw circles */
            //      rp=radii[gainscale];
            //      j=0;
            //      while ((r=* rp++* extr)) {
            //         po.x=r; po.y=0; po.z=0;
            //         proj(&po, c0);
            //         for (i=0;i<Npt;i++) {
            //            po.x=r* co[i];
            //po.z=r* si[i];
            //proj(&po, c1);
            //DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //c0[0]=c1[0]; c0[1]=c1[1];
            //         }
            //         po.x=0; po.z=r; proj(&po, c1);
            //DrawStringLL(c1[0], c1[1], labels[gainscale][j++]);
            //      }
            //   }
        }


        ///* ------------------------- draw complete picture ------------------------------------------------------------------------- */

        void Draw_axes()        /* draw the axes */
        {
            //    double c0[2], c1[2], c2[2], c3[2];
            //    Point pp;

            //    SetForeground(&c_axis);
            //    SetLineAttributes(0, GDK_LINE_SOLID, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //    pp.x = 0; pp.y = 0; pp.z = 0; proj(&pp, c0);
            //    pp.x = Axislen * extr; pp.y = 0; pp.z = 0; proj(&pp, c1);
            //    pp.x = 0; pp.y = Axislen * extr; pp.z = 0; proj(&pp, c2);
            //    pp.x = 0; pp.y = 0; pp.z = Axislen * extr; proj(&pp, c3);
            //    DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //    DrawLine(c0[0], c0[1], c2[0], c2[1]);
            //    DrawLine(c0[0], c0[1], c3[0], c3[1]);
            //    if (!quick)
            //    {
            //        if (c1[0] >= c0[0] || !ascending(c0, c1)) DrawStringLL(c1[0] + 2, c1[1], "X"); else DrawStringUL(c1[0] + 2, c1[1], "X");
            //        if (c2[0] >= c0[0] || !ascending(c0, c2)) DrawStringLL(c2[0] + 2, c2[1], "Y"); else DrawStringUL(c2[0] + 2, c2[1], "Y");
            //        if (c3[0] >= c0[0] || !ascending(c0, c3)) DrawStringLL(c3[0] + 2, c3[1], "Z"); else DrawStringUL(c3[0] + 2, c3[1], "Z");
            //    }

            //    SetLineAttributes(0, GDK_LINE_ON_OFF_DASH, GDK_CAP_ROUND, GDK_JOIN_ROUND);
            //    pp.x = 0; pp.y = 0; pp.z = 0; proj(&pp, c0);
            //    pp.x = -Axislen * extr; pp.y = 0; pp.z = 0; proj(&pp, c1);
            //    DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //    pp.x = 0; pp.y = -Axislen * extr; pp.z = 0; proj(&pp, c1);
            //    DrawLine(c0[0], c0[1], c1[0], c1[1]);
            //    pp.x = 0; pp.y = 0; pp.z = -Axislen * extr; proj(&pp, c1);
            //    DrawLine(c0[0], c0[1], c1[0], c1[1]);
        }


        public static void Draw_all(int ie)
        {
            //    char s[100];
            //    do
            //    {
            //        if (Pending() && ie) return;
            //        ClearWindow();
            //        if (gainplot == GPslice || gainplot == GPframe || gainplot == GPopaque) draw_gainscale(ie);
            //        draw_axes();
            //        if (Pending() && ie) break;
            //        if (structplot)
            //        {
            //            if (structplot == SPcurrents)
            //            {
            //                if (rp_index >= 0)
            //                {
            //                    draw_currents_noanim(ie, neco[rp_index].cu);
            //                    if (!quick) draw_phasecircle();
            //                }
            //            }
            //            else if (structplot == SPanim)
            //            {
            //                if (rp_index >= 0)
            //                {
            //                    draw_currents_anim(ie, neco[rp_index].cu);
            //                }
            //            }
            //            else draw_antenna(ie);
            //        }
            //        if (Pending() && ie) break;
            //        if (rp_index >= 0)
            //        {
            //            if (gainplot == GPslice || gainplot == GPframe || gainplot == GPopaque)
            //            {
            //                draw_gain(ie, neco[rp_index].rp);
            //                if (polarization == POLcolour) draw_polref();
            //            }
            //            if (gainplot == GPnearfield) draw_nearfield(ie, neco + rp_index);
            //            if (neco[rp_index].rp)
            //            {
            //                double maxgain, vgain;
            //                if (polarization == POLnone || polarization == POLcolour)
            //                {
            //                    maxgain = neco[rp_index].d[neco_maxgain],
            //               vgain = neco[rp_index].d[neco_vgain];
            //                }
            //                else
            //                {
            //                    maxgain = neco[rp_index].d[Neco_polgain + Neco_gsize * polarization];
            //                    vgain = neco[rp_index].d[neco_vgain2];
            //                }
            //                sprintf(s, "f = %g MHz   maxgain = %g dBi   vgain = %g dBi",
            //                   neco[rp_index].f, maxgain, vgain);
            //            }
            //            else
            //            {
            //                sprintf(s, "f = %g MHz", neco[rp_index].f);

            //            }
            //            SetForeground(&c_axis);
            //            DrawString(winsizex / 2, winsizey, s, 0.5, 0);
            //        }
            //    } while (0);
            //   out->Complete();
        }
    }
}
