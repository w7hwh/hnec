/*  XNECVIEW - a program for visualizing NEC2 input and output data
*
*  Copyright (C) 1998-2006,2011, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
*
*  Distributed on the conditions of version 2 of the GPL: see the files
*  README and COPYING, which accompany this source file.
*
*  This module parses NEC2's output.
*  Furthermore, it contains some functions for later processing of this
*  data.
*/
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    /// <summary>
    /// Handles reading the NEC output file
    /// </summary>
    class ReadNecOutput
    {
        double globalmaxdb = -1e30;


        //        void calcswr(ref NECoutput ne)
        //        {
        //            double zr, zi, gamma;

        //            zr = ne.d[neco_zr];
        //            zi = ne.[neco_zi];
        //            gamma = sqrt(((zr - r0) * (zr - r0) + zi * zi) / ((zr + r0) * (zr + r0) + zi * zi));
        //            ne.d[neco_swr] = (1 + gamma) / (1 - gamma);
        //        }


        void calcphiabs(ref Nec.NECoutput ne)
        {
            double zr, zi;
            //zr = ne.d[Nec.NECO.neco_zr];
            //zi = ne.d[Nec.NECO.neco_zi];
            //ne.d[Nec.NECO.neco_zabs] = Math.Sqrt(zr * zr + zi * zi);
            //ne.d[Nec.NECO.neco_zphi] = 180 * (1/Math.PI) * Math.Atan2(zi, zr);
        }


        //        void* mymalloc(size_t n)
        //        {
        //            void* p;
        //            p = malloc(n);
        //            if (p == NULL) { fprintf(stderr, "Out of memory\n"); exit(1); }
        //            return p;
        //        }

        //        void* myrealloc(void* p, size_t n)
        //        {
        //            p = realloc(p, n);
        //            if (p == NULL) { fprintf(stderr, "Out of memory\n"); exit(1); }
        //            return p;
        //        }


        //        double calc_polfactor(int pol, Gain* g)
        //        {
        //            /* note: for efficiency reasons, these formulas are also in update_maxgains() */
        //            double a, b, f;
        //            switch (pol)
        //            {
        //                case POLhor:
        //                    a = g->axial * g->axial;
        //                    b = sin(g->tilt);
        //                    f = (a + (1 - a) * b * b) / (1 + a);
        //                    break;
        //                case POLvert:
        //                    a = g->axial * g->axial;
        //                    b = cos(g->tilt);
        //                    f = (a + (1 - a) * b * b) / (1 + a);
        //                    break;
        //                case POLlhcp:
        //                    a = g->axial * g->axial;
        //                    f = (1 + 2 * g->axial + a) / 2 / (1 + a);
        //                    break;
        //                case POLrhcp:
        //                    a = g->axial * g->axial;
        //                    f = (1 - 2 * g->axial + a) / 2 / (1 + a);
        //                    break;
        //                default: f = 1; /* POLnone, POLcolour */
        //            }
        //            return f;
        //        }

        //        void update_maxgains(Gain* g, NECoutput* ne, double phi, double theta)
        //        {
        //            double f;
        //            double a, b;

        //            if (g->total > ne->d[neco_maxgain])
        //            {
        //                ne->d[neco_maxgain] = g->total;
        //                ne->d[neco_phi] = phi;
        //                ne->d[neco_theta] = theta;
        //            }

        //            a = g->axial * g->axial;

        //            f = (1 + 2 * g->axial + a) / 2 / (1 + a);
        //            f = g->total + 10 * log10(f);
        //            if (f > ne->d[neco_maxgain_lhcp])
        //            {
        //                ne->d[neco_maxgain_lhcp] = f;
        //                ne->d[neco_phi_lhcp] = phi;
        //                ne->d[neco_theta_lhcp] = theta;
        //            }

        //            f = (1 - 2 * g->axial + a) / 2 / (1 + a);
        //            f = g->total + 10 * log10(f);
        //            if (f > ne->d[neco_maxgain_rhcp])
        //            {
        //                ne->d[neco_maxgain_rhcp] = f;
        //                ne->d[neco_phi_rhcp] = phi;
        //                ne->d[neco_theta_rhcp] = theta;
        //            }

        //            b = sin(g->tilt);
        //            f = (a + (1 - a) * b * b) / (1 + a);
        //            f = g->total + 10 * log10(f);
        //            if (f > ne->d[neco_maxgain_hor])
        //            {
        //                ne->d[neco_maxgain_hor] = f;
        //                ne->d[neco_phi_hor] = phi;
        //                ne->d[neco_theta_hor] = theta;
        //            }

        //            b = cos(g->tilt);
        //            f = (a + (1 - a) * b * b) / (1 + a);
        //            f = g->total + 10 * log10(f);
        //            if (f > ne->d[neco_maxgain_vert])
        //            {
        //                ne->d[neco_maxgain_vert] = f;
        //                ne->d[neco_phi_vert] = phi;
        //                ne->d[neco_theta_vert] = theta;
        //            }
        //        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ne"></param>
        /// <param name="r"></param>
        void find_fb(ref Nec.NECoutput ne, ref Nec.Radpattern r)
        {
            int i, j, k;
            double theta, phi;

            //            theta = ne->d[neco_theta] * M_PI / 180;
            //            phi = ne->d[neco_phi] * M_PI / 180;
            //            for (k = POLnone; k <= POLrhcp; k++)
            //            {

            //                for (i = 0; i < r->numtheta; i++)
            //                {
            //                    double d;
            //                    d = fmod(fabs(r->gtheta[i] + theta), 2 * M_PI);
            //                    if (fabs(d - M_PI) < 0.00001) break;
            //                }
            //                for (j = 0; j < r->numphi; j++)
            //                {
            //                    double d;
            //                    d = fmod(fabs(r->gphi[j] - phi), 2 * M_PI);
            //                    if (fabs(d - M_PI) < 0.00001) break;
            //                }

            //                if (k == POLnone)
            //                {
            //                    if (i < r->numtheta && j < r->numphi)
            //                        ne->d[neco_fb] = ne->d[neco_maxgain] - r->gain[j][i].total;
            //                    else
            //                        ne->d[neco_fb] = -DBL_MAX;
            //                }
            //                else
            //                {
            //                    if (i < r->numtheta && j < r->numphi)
            //                    {
            //                        double f;
            //                        Gain* g = &r->gain[j][i];
            //                        ne->d[Neco_gsize * k + Neco_polfb2] = ne->d[Neco_gsize * k + Neco_polgain] - g->total;
            //                        f = calc_polfactor(k, g);
            //                        ne->d[Neco_gsize * k + Neco_polfb1] = ne->d[Neco_gsize * k + Neco_polgain] - r->gain[j][i].total + 10 * log10(f);
            //                    }
            //                    else
            //                    {
            //                        ne->d[Neco_gsize * k + Neco_polfb1] = -DBL_MAX;
            //                        ne->d[Neco_gsize * k + Neco_polfb2] = -DBL_MAX;
            //                    }
            //                }

            //                theta = ne->d[Neco_gsize * (k + 1) + Neco_poltheta] * M_PI / 180;
            //                phi = ne->d[Neco_gsize * (k + 1) + Neco_polphi] * M_PI / 180;
            //            }
        }

        /// <summary>
        /// tries to read (far field) radiation pattern data from f
        /// </summary>
        /// <param name="file"></param>
        /// <param name="ne"></param>
        private void Read_nec_output_rad(System.IO.StreamReader file, ref Nec.NECoutput ne)
        {
            char[] s = new char[200];
            double phi, theta;
            float db, axial, tilt;
            char polsense;
            int i;
            int end = 0;
            int maxthetas, maxphis;
            //            Radpattern* r;

            //            /* Some versions of NEC2 output extraneous data after "RADIATION PATTERNS" before the actual data
            //             * and column labels (e.g. RANGE = and EXP (-JKR) values ).  Discard until
            //             * the true bottom of the column labels (look for DB) */
            //            do fgets(s, 200, f); while (!feof(f) && !strstr(s, "DB"));
            //            if (feof(f)) return;

            //            /* allocate some memory; several of these arrays may need to be reallocated later on, if it turns out that their initial size is too small */
            //            r = mymalloc(sizeof(Radpattern));
            //            maxphis = 128;
            //            r->gain = mymalloc(maxphis * sizeof(Gain*));
            //            r->gphi = mymalloc(maxphis * sizeof(double));

            //            /* read the first block of radiation data, i.e., for one value of phi and the full range of theta */
            //            /* after this, we know how many thetas to expect per phi, which makes the memory allocation simpler */
            //            fgets(s, 200, f);
            //            removecommas(s);
            //            if (sscanf(s, "%lg%lg%*g%*g%g%g%g %c", &theta, &phi, &db, &axial, &tilt, &polsense) != 6) return;
            //            r->gphi[0] = phi;
            //            r->numphi = 1;
            //            r->numtheta = 0;
            //            maxthetas = 16;
            //            r->gtheta = mymalloc(maxthetas * sizeof(double));
            //            r->gain[0] = mymalloc(maxthetas * sizeof(Gain));
            //            do
            //            {
            //                Gain* g;
            //                if (r->numtheta >= maxthetas)
            //                {
            //                    maxthetas *= 2;
            //                    r->gtheta = myrealloc(r->gtheta, maxthetas * sizeof(double));
            //                    r->gain[0] = myrealloc(r->gain[0], maxthetas * sizeof(Gain));
            //                }
            //                g = &r->gain[r->numphi - 1][r->numtheta];
            //                g->total = db;
            //                g->tilt = tilt * M_PI / 180;
            //                if (polsense == 'R') axial = -axial;
            //                g->axial = axial;
            //                update_maxgains(g, ne, phi, theta);
            //                r->gtheta[r->numtheta++] = theta * M_PI / 180;
            //                fgets(s, 200, f);
            //                removecommas(s);
            //                if (sscanf(s, "%lg%lg%*g%*g%g%g%g %c", &theta, &phi, &db, &axial, &tilt, &polsense) != 6) { end = 1; break; }
            //            } while (phi == r->gphi[0]);
            //            r->gtheta = myrealloc(r->gtheta, r->numtheta * sizeof(double));
            //            r->gain[0] = myrealloc(r->gain[0], r->numtheta * sizeof(Gain));
            //            r->gphi[0] *= M_PI / 180;

            //            /* next, read the rest of the data, i.e., for the same values of theta and the entire range of phi */
            //            if (!end)
            //            {
            //                i = 0;
            //                while (1)
            //                {
            //                    Gain* g;
            //                    if (i == 0)
            //                    {
            //                        if (r->numphi >= maxphis)
            //                        {
            //                            maxphis *= 2;
            //                            r->gain = myrealloc(r->gain, maxphis * sizeof(Gain));
            //                            r->gphi = myrealloc(r->gphi, maxphis * sizeof(double));
            //                        }
            //                        r->gain[r->numphi] = mymalloc(r->numtheta * sizeof(Gain));
            //                        r->gphi[r->numphi++] = phi * M_PI / 180;
            //                    }
            //                    g = &r->gain[r->numphi - 1][i];
            //                    g->total = db;
            //                    g->tilt = tilt * M_PI / 180;
            //                    if (polsense == 'R') axial = -axial;
            //                    g->axial = axial;
            //                    update_maxgains(g, ne, phi, theta);
            //                    if (++i == r->numtheta) i = 0;
            //                    fgets(s, 200, f);
            //                    removecommas(s);
            //                    if (sscanf(s, "%lg%lg%*g%*g%g%g%g %c", &theta, &phi, &db, &axial, &tilt, &polsense) != 6) break;
            //                }
            //            }
            //            r->gain = myrealloc(r->gain, r->numphi * sizeof(Gain));
            //            r->gphi = myrealloc(r->gphi, r->numphi * sizeof(double));

            //            /* further processing of maximum gain info */
            //            if (ne->d[neco_maxgain] > globalmaxdb) globalmaxdb = ne->d[neco_maxgain];
            //            find_fb(ne, r);

            //            /* allocate arrays which will later be filled with the grid coordinates in 3D space */
            //            r->gpo = mymalloc(r->numphi * sizeof(Point*));
            //            for (i = 0; i < r->numphi; i++) r->gpo[i] = mymalloc(r->numtheta * sizeof(Point));

            //            /* allocate and fill the arrays of sines and cosines of theta and phi */
            //            r->sin_gphi = mymalloc(r->numphi * sizeof(double));
            //            r->cos_gphi = mymalloc(r->numphi * sizeof(double));
            //            for (i = 0; i < r->numphi; i++)
            //            {
            //                r->sin_gphi[i] = sin(r->gphi[i]);
            //                r->cos_gphi[i] = cos(r->gphi[i]);
            //            }
            //            r->sin_gtheta = mymalloc(r->numtheta * sizeof(double));
            //            r->cos_gtheta = mymalloc(r->numtheta * sizeof(double));
            //            for (i = 0; i < r->numtheta; i++)
            //            {
            //                r->sin_gtheta[i] = sin(r->gtheta[i]);
            //                r->cos_gtheta[i] = cos(r->gtheta[i]);
            //            }

            //            /* finally, attach the set of radiation data just read to the *ne structure */
            //            r->next = NULL;
            //            if (!ne->rp) ne->rp = r;
            //            else
            //            {
            //                Radpattern* p;
            //                p = ne->rp;
            //                while (p->next) p = p->next;
            //                p->next = r;
            //            }
        }

        /// <summary>
        /// tries to read segmentation data from f
        /// </summary>
        /// <param name="f"></param>
        /// <returns>returns NULL if not succesful</returns>
        Nec.Currents Read_nec_output_seg(System.IO.StreamReader f)           
        {
            //            char s[200];
            //            char* p;
            Nec.Currents cu = new Nec.Currents();

            //            /* skip some column labels etc., until we find a line that starts with a digit */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ') p++;
            //            } while (!isdigit(*p) && !feof(f));
            //            if (feof(f)) return NULL;

            //            /* allocate memory */
            //            cu = mymalloc(sizeof(Currents));
            //            cu->numseg = 0;
            //            cu->maxseg = 1;
            //            cu->s = mymalloc(cu->maxseg * sizeof(*(cu->s)));

            //            cu->maxI = 0;
            //            cu->maxQ = 0;

            //            /* read the segmentation geometry */
            //            do
            //            {
            //                double x, y, z, l, alpha, beta;
            //                double dx, dy, dz;
            //                Segcur* se;
            //                if (sscanf(s, "%*d%lg%lg%lg%lg%lg%lg", &x, &y, &z, &l, &alpha, &beta) != 6) break;
            //                EXPAND_IF_NECESSARY(cu->numseg, cu->maxseg, cu->s)
            //               se = cu->s + cu->numseg;
            //                alpha *= M_PI / 180;
            //                beta *= M_PI / 180;
            //                l /= 2;
            //                dx = l * cos(alpha) * cos(beta);
            //                dy = l * cos(alpha) * sin(beta);
            //                dz = l * sin(alpha);
            //                se->c.x = x;
            //                se->c.y = y;
            //                se->c.z = z;
            //                se->p0.x = x - dx;
            //                se->p0.y = y - dy;
            //                se->p0.z = z - dz;
            //                se->p1.x = x + dx;
            //                se->p1.y = y + dy;
            //                se->p1.z = z + dz;
            //                updateextremes(&se->c);
            //                se->re[0] = dx / l;
            //                se->re[1] = dy / l;
            //                se->re[2] = dz / l;
            //                if (cos(alpha) == 0)
            //                {
            //                    se->q0[0] = 1; se->q0[1] = 0; se->q0[2] = 0;
            //                    se->q1[0] = 0; se->q1[1] = 1; se->q1[2] = 0;
            //                }
            //                else
            //                {
            //                    se->q0[0] = sin(alpha) * cos(beta); se->q0[1] = sin(alpha) * sin(beta); se->q0[2] = -cos(alpha);
            //                    se->q1[0] = sin(beta); se->q1[1] = -cos(beta); se->q1[2] = 0;
            //                }
            //                se->qre = 0;
            //                se->qim = 0;
            //                cu->numseg++;
            //                fgets(s, 200, f);
            //                removecommas(s);
            //            } while (!feof(f));

            //            cu->numrealseg = cu->numseg;
            //            cu->numanimseg = cu->numseg;

            return cu;
        }

        /// <summary>
        /// tries to read surface patch data from f
        /// </summary>
        /// <param name="f"></param>
        /// <param name="cu"></param>
        void Read_nec_output_patches(System.IO.StreamReader f, ref Nec.Currents cu)
        {
            //            char s[200];
            //            char* p;

            //            /* skip some column labels etc., until we find a line that starts with a digit */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ') p++;
            //            } while (!isdigit(*p) && !feof(f));
            //            if (feof(f)) return;

            //            /* read the surface patch geometry */
            //            do
            //            {
            //                double x, y, z, area;
            //                if (sscanf(s, "%*d%lg%lg%lg%*g%*g%*g%lg", &x, &y, &z, &area) != 4) break;
            //                EXPAND_IF_NECESSARY(cu->numseg, cu->maxseg, cu->s)
            //               cu->s[cu->numseg].c.x = x;
            //                cu->s[cu->numseg].c.y = y;
            //                cu->s[cu->numseg].c.z = z;
            //                updateextremes(&cu->s[cu->numseg].c);
            //                cu->s[cu->numseg].area = area;
            //                cu->numseg++;
            //                fgets(s, 200, f);
            //                removecommas(s);
            //            } while (!feof(f));
            //            cu->numanimseg = cu->numseg;
        }

        /// <summary>
        /// tries to read segment currents from f
        /// </summary>
        /// <param name="f"></param>
        /// <param name="cug"></param>
        /// <returns></returns>
        Nec.Currents read_nec_output_currents(System.IO.StreamReader f, ref Nec.Currents cug)
        {
            //            char s[200];
            //            char* p;
            Nec.Currents cu = new Nec.Currents();
            //            int i;

            //            /* skip column labels etc., until we find a line that starts with a digit */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ') p++;
            //            } while (!isdigit(*p) && !feof(f));
            //            if (feof(f)) return NULL;

            //            /* allocate memory and copy the geometry info that we collected earlier: */
            //            cu = mymalloc(sizeof(Currents));
            //            memcpy(cu, cug, sizeof(Currents));
            //            cu->maxseg = cu->numseg;
            //            cu->s = mymalloc(cu->maxseg * sizeof(*(cu->s)));
            //            memcpy(cu->s, cug->s, cu->maxseg * sizeof(*(cu->s)));

            //            /* read the amplitude and phase for every segment */
            //            i = 0;
            //            do
            //            {
            //                int j;
            //                if (sscanf(s, "%d%*d%*g%*g%*g%*g%*g%*g%g%g", &j, &cu->s[i].a, &cu->s[i].phase) != 3) break;
            //                if (i != j - 1) break;
            //                for (j = 0; j < 3; j++)
            //                {
            //                    cu->s[i].im[j] = cu->s[i].re[j] * cu->s[i].a * sin(cu->s[i].phase * M_PI / 180);
            //                    cu->s[i].re[j] = cu->s[i].re[j] * cu->s[i].a * cos(cu->s[i].phase * M_PI / 180);
            //                }
            //                if (cu->s[i].a > cu->maxI) cu->maxI = cu->s[i].a;
            //                i++;
            //                fgets(s, 200, f);
            //            } while (i < cu->numrealseg && !feof(f));
            //            if (i != cu->numrealseg) { fprintf(stderr, "Segment currents data in output file incorrect or incomplete\n"); free(cu->s); free(cu); return NULL; }

            return cu;
        }

        /// <summary>
        /// tries to read segment charges from f
        /// </summary>
        /// <param name="f"></param>
        /// <param name="cu"></param>
        void Read_nec_output_charges(System.IO.StreamReader f, ref Nec.Currents cu)
        {
            //            char s[200];
            //            char* p;
            //            int i;

            //            /* skip column labels etc., until we find a line that starts with a digit */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ') p++;
            //            } while (!isdigit(*p) && !feof(f));
            //            if (feof(f)) return;

            //            /* read the amplitude and phase for every segment */
            //            i = 0;
            //            do
            //            {
            //                int j;
            //                double a;
            //                double l;
            //                if (sscanf(s, "%d%*d%*g%*g%*g%*g%g%g%lg", &j, &cu->s[i].qre, &cu->s[i].qim, &a) != 4) break;
            //                if (i != j - 1) break;
            //                /* note: we read charge densities, so we need to multiply them by the segment length to obtain the charges themselves */
            //                l = 0;
            //                l += (cu->s[i].p1.x - cu->s[i].p0.x) * (cu->s[i].p1.x - cu->s[i].p0.x);
            //                l += (cu->s[i].p1.y - cu->s[i].p0.y) * (cu->s[i].p1.y - cu->s[i].p0.y);
            //                l += (cu->s[i].p1.z - cu->s[i].p0.z) * (cu->s[i].p1.z - cu->s[i].p0.z);
            //                l = sqrt(l);
            //                cu->s[i].qre *= l;
            //                cu->s[i].qim *= l;
            //                a *= l;
            //                if (a > cu->maxQ) cu->maxQ = a;
            //                i++;
            //                fgets(s, 200, f);
            //            } while (i < cu->numrealseg && !feof(f));
            //            if (i != cu->numrealseg) { fprintf(stderr, "Segment charges data in output file incorrect or incomplete\n"); return; }

            //            return;
        }

        /// <summary>
        /// tries to read surface patch currents from f
        /// </summary>
        /// <param name="f"></param>
        /// <param name="cu"></param>
        /// <returns></returns>
        Nec.Currents read_nec_output_patchcurrents(System.IO.StreamReader f, ref Nec.Currents cu)
        {
            //            char s[200];
            //            char* p;
            //            int i;
            //            int oldnumseg;

            //            /* skip some column labels etc., until we find a line that starts with a digit */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ') p++;
            //            } while (!isdigit(*p) && !feof(f));
            //            if (feof(f)) return NULL;

            //            /* skip the line containing only the patch number */
            //            fgets(s, 200, f);

            //            /* read the amplitude and phase for every segment */
            //            i = cu->numrealseg;
            //            oldnumseg = cu->numseg;
            //            EXPAND_IF_NECESSARY(2 * cu->numseg - cu->numrealseg, cu->maxseg, cu->s)
            //           do
            //            {
            //                double xr, xi, yr, yi, zr, zi;  /* real and imaginary components of current density in X, Y and Z direction */
            //                double xr2, xi2, yr2, yi2, zr2, zi2;  /* and their squares */
            //                double totl;   /* length of the segment that we're going to use to represent the patch current */
            //                double dxr, dyr, dzr;
            //                double dxi, dyi, dzi;
            //                double h, hh;
            //                double ta, si, co, phi;
            //#if 0
            //      double dx,dy,dz;
            //      double l;
            //#endif

            //                if (sscanf(s, "%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg", &xr, &xi, &yr, &yi, &zr, &zi) != 6) break;
            //                totl = sqrt(cu->s[i].area) / 2;  /* in principle we can choose any length we wish, but this one gives a nice picture */

            //                cu->s[i].re[0] = xr; cu->s[i].im[0] = xi;
            //                cu->s[i].re[1] = yr; cu->s[i].im[1] = yi;
            //                cu->s[i].re[2] = zr; cu->s[i].im[2] = zi;

            //#if 0
            //      /* This (disabled) piece of code just represents the currents on the surface patch by two
            //         segment currents corresponding to the real and imaginary vectors supplied by NEC.
            //         The code is here only for testing purposes.
            //      */
            //      xr2=xr*xr; xi2=xi*xi;
            //      yr2=yr*yr; yi2=yi*yi;
            //      zr2=zr*zr; zi2=zi*zi;

            //      cu->s[cu->numseg].c = cu->s[i].c;

            //      hh=sqrt(xr2+xi2+yr2+yi2+zr2+zi2);
            //      cu->s[i].a=         hh*cu->s[i].area/totl;
            //      cu->s[cu->numseg].a=hh*cu->s[i].area/totl;

            //      h=sqrt(xr2+yr2+zr2);
            //      l=totl*h/hh;
            //      dx=xr*totl/hh/2; dy=yr*totl/hh/2; dz=zr*totl/hh/2;
            //      cu->s[i].phase=0;
            //      cu->s[i].p0.x = cu->s[i].c.x-dx; cu->s[i].p0.y = cu->s[i].c.y-dy; cu->s[i].p0.z = cu->s[i].c.z-dz;
            //      cu->s[i].p1.x = cu->s[i].c.x+dx; cu->s[i].p1.y = cu->s[i].c.y+dy; cu->s[i].p1.z = cu->s[i].c.z+dz;

            //      h=sqrt(xi2+yi2+zi2);
            //      l=totl*h/hh;
            //      dx=xi*totl/hh/2; dy=yi*totl/hh/2; dz=zi*totl/hh/2;
            //      cu->s[cu->numseg].phase=90;
            //      cu->s[cu->numseg].p0.x = cu->s[cu->numseg].c.x-dx; cu->s[cu->numseg].p0.y = cu->s[cu->numseg].c.y-dy; cu->s[cu->numseg].p0.z = cu->s[cu->numseg].c.z-dz;
            //      cu->s[cu->numseg].p1.x = cu->s[cu->numseg].c.x+dx; cu->s[cu->numseg].p1.y = cu->s[cu->numseg].c.y+dy; cu->s[cu->numseg].p1.z = cu->s[cu->numseg].c.z+dz;

            //      cu->numseg++;
            //#endif

            //#if 1
            //      /* In general, the current on a surface patch can be considered as "running around"
            //         an ellipse. In many linearly polarized antennas, this ellipse collapses to a
            //         straight line.
            //         The below code represents the "elliptic" surface current by two segment currents
            //         oriented along the major and minor axes of the ellipse. If the ellipse collapses to
            //         a straight line, the minor axis component becomes zero, leaving one segment aligned
            //         according to the real direction of the current in the surface patch.
            //         Particularly in the latter case, this is much more natural than taking the "raw"
            //         real and imaginary current vectors supplied by NEC, which in this case would have
            //         the same direction and thus be indistinguishable in the picture.
            //      */

            //      xr2=xr*xr; xi2=xi*xi;
            //      yr2=yr*yr; yi2=yi*yi;
            //      zr2=zr*zr; zi2=zi*zi;

            //      cu->s[i].a=sqrt(xr2+yr2+zr2+xi2+yi2+zi2)*cu->s[i].area/totl;
            //      if (cu->s[i].a>cu->maxI) cu->maxI=cu->s[i].a;

            //      hh=xr2+yr2+zr2-xi2-yi2-zi2;
            //      if (hh==0) {
            //         co = sqrt(0.5);
            //         si = sqrt(0.5);
            //         phi = M_PI/4;
            //      } else {
            //         ta = 2*(xr*xi+yr*yi+zr*zi)/(xr2+yr2+zr2-xi2-yi2-zi2);
            //         co = sqrt(0.5+0.5*sqrt(1/(ta*ta+1)));
            //         si = sqrt(0.5-0.5*sqrt(1/(ta*ta+1)));
            //         if (ta<0) si=-si;
            //         phi = atan(ta)/2;
            //      }
            //      dxr =  xr*co + xi*si;  dyr =  yr*co + yi*si;  dzr =  zr*co + zi*si;
            //      dxi = -xr*si + xi*co;  dyi = -yr*si + yi*co;  dzi = -zr*si + zi*co;
            //      h = totl/sqrt(xr2+yr2+zr2+xi2+yi2+zi2)/2;

            //      cu->s[i].phase = phi*180/M_PI;
            //      cu->s[i].p0.x = cu->s[i].c.x-dxr*h;  cu->s[i].p0.y = cu->s[i].c.y-dyr*h;  cu->s[i].p0.z = cu->s[i].c.z-dzr*h;
            //      cu->s[i].p1.x = cu->s[i].c.x+dxr*h;  cu->s[i].p1.y = cu->s[i].c.y+dyr*h;  cu->s[i].p1.z = cu->s[i].c.z+dzr*h;

            //      cu->s[cu->numseg].a = cu->s[i].a;
            //      cu->s[cu->numseg].c = cu->s[i].c;
            //      cu->s[cu->numseg].phase = phi*180/M_PI+90;
            //      cu->s[cu->numseg].p0.x = cu->s[cu->numseg].c.x-dxi*h;  cu->s[cu->numseg].p0.y = cu->s[cu->numseg].c.y-dyi*h;  cu->s[cu->numseg].p0.z = cu->s[cu->numseg].c.z-dzi*h;
            //      cu->s[cu->numseg].p1.x = cu->s[cu->numseg].c.x+dxi*h;  cu->s[cu->numseg].p1.y = cu->s[cu->numseg].c.y+dyi*h;  cu->s[cu->numseg].p1.z = cu->s[cu->numseg].c.z+dzi*h;

            //      cu->numseg++;
            //#endif

            //                i++;
            //                fgets(s, 200, f);  /* again, skip the line containing only the patch number */
            //                fgets(s, 200, f);
            //            } while (i < oldnumseg && !feof(f));
            //            if (i != oldnumseg) { fprintf(stderr, "Surface patch currents data in output file incorrect\n"); free(cu->s); free(cu); return NULL; }

            return cu;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="cu"></param>
        void Normalize_currents(ref Nec.Currents cu)
        {
            int i;
            double max = 0;
            //for (i = 0; i < cu.numseg; i++) if (max < cu.s[i].a) max = cu.s[i].a;
            //for (i = 0; i < cu.numseg; i++) cu.s[i].a /= max;
        }

        /// <summary>
        /// tries to read near electric or magnetic field from f
        /// </summary>
        /// <param name="f"></param>
        /// <param name="ne"></param>
        /// <param name="magnetic"></param>
        void Read_nec_output_nearfield(System.IO.StreamReader f, ref Nec.NECoutput ne, int magnetic)
        {
            //            char s[200];
            //            char* p;
            //            Nearfield* nf,*np;

            //            /* skip some column labels etc., until we find a line that starts with a digit, possibly preceeded by a minus */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                p = s;
            //                while (*p == ' ' || *p == '-') p++;
            //            } while (!isdigit(*p) && !feof(f));

            //            do
            //            {
            //                double mag[3], phase[3]; /* note: "mag" is an abbreviation for magnitude, not for magnetic... */
            //                Point p;
            //                int i;
            //                double l;

            //                /* parse a line of data */
            //                if (sscanf(s, "%g%g%g%lg%lg%lg%lg%lg%lg",
            //                   &p.x, &p.y, &p.z,
            //                   &mag[0], &phase[0],
            //                   &mag[1], &phase[1],
            //                   &mag[2], &phase[2]) != 9) break;

            //                np = ne->nf;
            //                if (np == NULL)
            //                {
            //                    /* no near field data yet? then this is the first point */
            //                    nf = mymalloc(sizeof(Nearfield));
            //                    if (magnetic) nf->evalid = 0; else nf->hvalid = 0;
            //                    ne->nf = nf;
            //                    nf->next = NULL;
            //                    nf->p = p;
            //                }
            //                else
            //                {
            //                    /* otherwise, search whether this point is already in the data set */
            //                    do
            //                    {
            //                        if (np->p.x == p.x && np->p.y == p.y && np->p.z == p.z)
            //                        {
            //                            /* if so, just copy the new data into the existing record */
            //                            nf = np;
            //                            break;
            //                        }
            //                        if (!np->next)
            //                        {
            //                            /* otherwise, append a new record to the set */
            //                            nf = mymalloc(sizeof(Nearfield));
            //                            if (magnetic) nf->evalid = 0; else nf->hvalid = 0;
            //                            np->next = nf;
            //                            nf->next = NULL;
            //                            nf->p = p;
            //                            break;
            //                        }
            //                        np = np->next;
            //                    } while (1);
            //                }
            //                if (magnetic) nf->hvalid = 1; else nf->evalid = 1;
            //                l = 0;
            //                for (i = 0; i < 3; i++)
            //                {
            //                    if (magnetic)
            //                    {
            //                        nf->hre[i] = mag[i] * cos(phase[i] * M_PI / 180); l += nf->hre[i] * nf->hre[i];
            //                        nf->him[i] = mag[i] * sin(phase[i] * M_PI / 180); l += nf->him[i] * nf->him[i];
            //                    }
            //                    else
            //                    {
            //                        nf->ere[i] = mag[i] * cos(phase[i] * M_PI / 180); l += nf->ere[i] * nf->ere[i];
            //                        nf->eim[i] = mag[i] * sin(phase[i] * M_PI / 180); l += nf->eim[i] * nf->eim[i];
            //                    }
            //                }
            //                l = sqrt(l);
            //                if (magnetic)
            //                {
            //                    if (l > ne->maxh) ne->maxh = l;
            //                }
            //                else
            //                {
            //                    if (l > ne->maxe) ne->maxe = l;
            //                }

            //                /* read next line */
            //                fgets(s, 200, f);
            //            } while (1);
        }

        /// <summary>
        /// Tries to read NEC output data from filename, updates gain distribution arrays.
        /// </summary>
        /// <param name="filename">Name of file to read from.</param>
        /// <returns>returns true if succesful</returns>
        public Boolean Read_nec_output(string filename)
        {
            System.IO.StreamReader file = new System.IO.StreamReader(filename);
            return Read_nec_output(file);
        }

        /// <summary>
        /// Tries to read NEC output data from StreamReader file, updates gain distribution arrays.
        /// </summary>
        /// <param name="file">StreamReader to read from.</param>
        /// <returns>true if successful</returns>
        public Boolean Read_nec_output(System.IO.StreamReader file)
        {
            //            char s[200];
            //            char* p;
            double zr, zi;
            int new_rp_index = -1;
            int near_index = -1;
            int curr_index = -1;
            //            NECoutput* ne;
            double freq;
            int i;
            //            Currents* cu = NULL;
            int oldnumneco;

            //oldnumneco = numneco;

            //            /* skip all lines until a line containing "STRUCTURE SPECIFICATION": this
            //               effectively skips all comments (from CM cards) and text resulting from
            //               reading Numerical Green's Function parts. Of course, if a user writes
            //               "STRUCTURE SPECIFICATION" in his comment lines, this still fails...
            //            */
            //            do
            //            {
            //                fgets(s, 200, f);
            //            } while (!feof(f) && !strstr(s, "STRUCTURE SPECIFICATION"));
            //            if (feof(f)) return false;    // failed

            //            /* check whether segmentation and surface patch data is included, and if so, read it */
            //            do
            //            {
            //                fgets(s, 200, f);
            //                if (strstr(s, "SEGMENTATION DATA")) cu = read_nec_output_seg(f);
            //                if (cu && strstr(s, "SURFACE PATCH DATA")) read_nec_output_patches(f, cu);
            //            } while (!feof(f) && !strstr(s, "FREQUENCY"));

            //            printf("#  freq.       Zr       Zi      SWR     gain      f/b      phi    theta\n");
            //            while (!feof(f))
            //            {

            //                EXPAND_IF_NECESSARY(numneco, maxfreqs, neco)

            //      /* search for the frequency and read it */
            //      do fgets(s, 200, f); while (!feof(f) && !strstr(s, "FREQUENCY=") && !strstr(s, "FREQUENCY : "));
            //                if (feof(f)) break;
            //                p = strchr(s, '=');
            //                if (!p) p = strchr(s, ':');
            //                freq = atof(p + 1);

            //                /* find the right place in memory, and prepare it if needed */
            //                if (numneco == 0 || freq > neco[numneco - 1].f)
            //                {
            //                    ne = neco + numneco;
            //                    numneco++;
            //                    ne->f = 0;
            //                }
            //                else
            //                {
            //                    for (i = 0, ne = neco; i < numneco; i++, ne++)
            //                        if (ne->f >= freq)
            //                        {
            //                            if (ne->f == freq) break;
            //                            memmove(ne + 1, ne, (neco + numneco - ne) * sizeof(NECoutput));
            //                            numneco++;
            //                            if (new_rp_index >= i) new_rp_index++;
            //                            break;
            //                        }
            //                }
            //                if (ne->f != freq)
            //                {
            //                    ne->f = freq;
            //                    ne->rpgpovalid = 0;
            //                    ne->rp = NULL;
            //                    ne->cu = NULL;
            //                    ne->nf = NULL;
            //                    ne->maxe = ne->maxh = 0;
            //                    ne->d[neco_maxgain] = -DBL_MAX;
            //                    ne->d[neco_maxgain_hor] = -DBL_MAX;
            //                    ne->d[neco_maxgain_vert] = -DBL_MAX;
            //                    ne->d[neco_maxgain_lhcp] = -DBL_MAX;
            //                    ne->d[neco_maxgain_rhcp] = -DBL_MAX;
            //                    ne->d[neco_zr] = 0;
            //                    ne->d[neco_zi] = 0;
            //                    ne->d[neco_swr] = 0;
            //                    ne->d[neco_zphi] = 0;
            //                    ne->d[neco_zabs] = 0;
            //                }

            //                /* find and read the input impedance, and calculate the SWR */
            //#if 0
            //      do fgets(s,200,f);  while (!feof(f) && !strstr(s,"ANTENNA INPUT PARAMETERS"));
            //      if (feof(f)) break;
            //      do {  /* skip lines until a line starting with a digit is found */
            //         fgets(s,200,f);  
            //         p=s;
            //         while (*p==' ') p++;
            //      } while (!feof(f) && !isdigit(*p));
            //      if ((p=strchr(s,'*'))) *p=' ';
            //      if (sscanf(s,"%*i%*i%*g%*g%*g%*g%lg%lg",&zr,&zi)!=2) break;
            //      ne->d[neco_zr]=zr;
            //      ne->d[neco_zi]=zi;
            //      calcswr(ne);
            //      calcphiabs(ne);
            //#endif

            //                /* check whether radiation pattern data is included for this frequency, and if so, read it */
            //                /* check also whether input impedance data is included for this frequency, and if so, read it */
            //                /* check also whether currents and/or charge data is included for this frequency, and if so, read it */
            //                /* check also whether near field data is included for this frequency, and if so, read it */
            //                do
            //                {
            //                    fgets(s, 200, f);
            //                    if (cu && strstr(s, "CURRENTS AND LOCATION"))
            //                    {
            //                        if (ne->cu)
            //                        {
            //                            fprintf(stderr, "File contains more than one set of currents at a frequency; ignoring earlier data\n");
            //                            free(ne->cu->s);
            //                            free(ne->cu);
            //                        }
            //                        ne->cu = read_nec_output_currents(f, cu);
            //                        if (ne->cu && curr_index < 0) curr_index = ne - neco;
            //                    }
            //                    if (ne->cu && strstr(s, "SURFACE PATCH CURRENTS")) ne->cu = read_nec_output_patchcurrents(f, ne->cu);
            //                    if (cu && strstr(s, "CHARGE DENSITIES"))
            //                    {
            //                        if (!ne->cu) fprintf(stderr, "Charge density information should come after currents information; ignoring it\n");
            //                        else read_nec_output_charges(f, ne->cu);
            //                    }
            //                    if (strstr(s, "NEAR ELECTRIC FIELDS"))
            //                    {
            //                        read_nec_output_nearfield(f, ne, 0);
            //                        if (near_index < 0) near_index = ne - neco;
            //                    }
            //                    if (strstr(s, "NEAR MAGNETIC FIELDS"))
            //                    {
            //                        read_nec_output_nearfield(f, ne, 1);
            //                        if (near_index < 0) near_index = ne - neco;
            //                    }
            //                    if (strstr(s, "RADIATION PATTERNS"))
            //                    {
            //                        read_nec_output_rad(f, ne);
            //                        if (!ne->rp) return 1;
            //                        if (new_rp_index < 0) new_rp_index = ne - neco;
            //                    }

            //                    if (strstr(s, "ANTENNA INPUT PARAMETERS"))
            //                    {
            //                        do
            //                        {  /* skip lines until a line starting with a digit is found */
            //                            fgets(s, 200, f);
            //                            p = s;
            //                            while (*p == ' ') p++;
            //                        } while (!feof(f) && !isdigit(*p));
            //                        if ((p = strchr(s, '*'))) *p = ' ';
            //                        if (sscanf(s, "%*i%*i%*g%*g%*g%*g%lg%lg", &zr, &zi) != 2) continue;
            //                        ne->d[neco_zr] = zr;
            //                        ne->d[neco_zi] = zi;
            //                        calcswr(ne);
            //                        calcphiabs(ne);
            //                    }

            //                } while (!feof(f) && !strstr(s, "FREQUENCY"));

            //                if (ne->cu) normalize_currents(ne->cu);

            //                /* print the data */
            //                printf("%8g %8g %8g %8g ", ne->f, zr, zi, ne->d[neco_swr]);
            //                if (ne->rp)
            //                    if (ne->d[neco_fb] > -DBL_MAX)
            //                        printf("%8g %8g %8g %8g\n", ne->d[neco_maxgain], ne->d[neco_fb], ne->d[neco_phi], ne->d[neco_theta]);
            //                    else
            //                        printf("%8g        - %8g %8g\n", ne->d[neco_maxgain], ne->d[neco_phi], ne->d[neco_theta]);
            //                else
            //                    printf("       -        -        -        -\n");
            //            }

            //            if (cu) { free(cu->s); free(cu); }

            //            if (rp_index < 0 || rp_index >= numneco)
            //            {
            //                if (new_rp_index >= 0) rp_index = new_rp_index;
            //                else if (near_index >= 0) rp_index = near_index;
            //                else if (curr_index >= 0) rp_index = curr_index;
            //            }
            //            if (neco[rp_index].rp == NULL && new_rp_index >= 0) rp_index = new_rp_index;

            //            return (numneco == oldnumneco);
            return true;    // successful
        }


        //#define pow10(x) exp(M_LN10*(x))

        /// <summary>
        /// transform gain distribution into an array of points in 3D space
        /// </summary>
        /// <param name="ne"></param>
        void Process_nec_output(ref Nec.NECoutput ne)
        {
            //            int i, j;
            //            double r = 0;
            //            Radpattern* rp;

            //            if (ne->rpgpovalid) return;
            //            rp = ne->rp;
            //            if (extr_str == 0) extr = 1;
            //            else extr = extr_str * GAINSIZE;
            //            while (rp)
            //            {
            //                for (i = 0; i < rp->numtheta; i++)
            //                    for (j = 0; j < rp->numphi; j++)
            //                    {
            //                        double db;
            //                        double f;
            //                        Gain* g = &rp->gain[j][i];
            //                        f = calc_polfactor(polarization, g);
            //                        if (f < 1e-200) r = 0;
            //                        else
            //                        {
            //                            db = g->total - globalmaxdb;
            //                            switch (gainscale)
            //                            {
            //                                case GSlinpower: r = f * pow10(db / 10); break;        /* linear in power */
            //                                case GSlinvolt: r = sqrt(f) * pow10(db / 20); break;        /* linear in voltage */
            //                                case GSarrl: r = exp(0.116534 * (db + 10 * log10(f)) / 2); break;  /* arrl */
            //                                case GSlog: r = db + 10 * log10(f); if (r < -40) r = 0; else r = r / 40 + 1; break;  /* log */
            //                            }
            //                            r *= extr;
            //                        }
            //                        if (r < 1e-6) r = 1e-6;  /* prevent points from being just about exactly at the origin, since this confuses the edge hiding calculations */
            //                        rp->gpo[j][i].x = rp->cos_gphi[j] * rp->sin_gtheta[i] * r;
            //                        rp->gpo[j][i].y = rp->sin_gphi[j] * rp->sin_gtheta[i] * r;
            //                        rp->gpo[j][i].z = rp->cos_gtheta[i] * r;
            //                    }
            //                rp = rp->next;
            //            }
            //            ne->rpgpovalid = 1;
        }

        /// <summary>
        /// update the vgain records in accordance with the viewing direction
        /// </summary>
        /// <param name=""></param>
        void Calc_vgain()
        {
            //            NECoutput* ne;
            //            int i, j, k;
            //            int i0 = 0, j0 = 0;
            //            Radpattern* rp,*rp0;
            //            double dif, d, d0;
            //            double ph, th;

            //            ph = fmod(phi, 360.) * (M_PI / 180.);
            //            ph += 5 * M_PI;
            //            th = theta * (M_PI / 180.);

            //            for (k = 0, ne = neco; k < numneco; k++, ne++)
            //            {
            //                rp = ne->rp;
            //                ne->d[neco_vgain] = -DBL_MAX;
            //                ne->d[neco_vgain2] = -DBL_MAX;
            //                dif = DBL_MAX;
            //                rp0 = NULL;
            //                while (rp)
            //                {
            //                    d0 = DBL_MAX;
            //                    i0 = 0;
            //                    for (i = 0; i < rp->numtheta; i++)
            //                    {
            //                        d = fabs(rp->gtheta[i] - th);
            //                        if (d < d0) { d0 = d; i0 = i; }
            //                    }
            //                    if (d0 < dif)
            //                    {
            //                        j0 = 0;
            //                        for (j = 0; j < rp->numphi; j++)
            //                        {
            //                            d = d0 + fabs(fmod(ph - rp->gphi[j], 2 * M_PI) - M_PI);
            //                            if (d < dif)
            //                            {
            //                                Gain* g;
            //                                dif = d; j0 = j; rp0 = rp;
            //                                g = &rp->gain[j0][i0];
            //                                ne->d[neco_vgain] = g->total;
            //                                ne->d[neco_vgain2] = g->total + 10 * log10(calc_polfactor(polarization, g));
            //                            }
            //                        }
            //                    }
            //                    rp = rp->next;
            //                }
            //                if (rp0)
            //                {
            //                    double ph0, th0;
            //                    th0 = rp0->gtheta[i0];
            //                    ph0 = rp0->gphi[j0];
            //                    for (i = 0; i < rp0->numtheta; i++)
            //                    {
            //                        double d;
            //                        d = fmod(fabs(rp0->gtheta[i] + th0), 2 * M_PI);
            //                        if (fabs(d - M_PI) < 0.00001) break;
            //                    }
            //                    for (j = 0; j < rp0->numphi; j++)
            //                    {
            //                        double d;
            //                        d = fmod(fabs(rp0->gphi[j] - ph0), 2 * M_PI);
            //                        if (fabs(d - M_PI) < 0.00001) break;
            //                    }
            //                    if (i < rp0->numtheta && j < rp0->numphi)
            //                    {
            //                        Gain* g;
            //                        g = &rp0->gain[j][i];
            //                        ne->d[neco_vfb] = ne->d[neco_vgain2] - g->total + 10 * log10(calc_polfactor(polarization, g));
            //                        ne->d[neco_vfb2] = ne->d[neco_vgain2] - g->total;
            //                    }
            //                    else
            //                    {
            //                        ne->d[neco_vfb] = -DBL_MAX;
            //                        ne->d[neco_vfb2] = -DBL_MAX;
            //                    }
            //                }
            //            }
        }

        /// <summary>
        /// mark gpo fields as invalid; must be called after every change of gain scale
        /// </summary>
        void Mark_gpo_invalid()
        {
            //NECoutput* ne;
            //int k;
            //for (k = 0, ne = neco; k < numneco; k++, ne++) ne->rpgpovalid = 0;
        }
    }
}
