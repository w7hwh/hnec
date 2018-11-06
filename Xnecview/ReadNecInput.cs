/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 1998-2006,2011, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 *
 */
//
// CM CE GA GE GF GM GR GS GW GX SP SM SC 
//

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Text.RegularExpressions;

#if (DEBUG)
using System.Diagnostics;
#endif

namespace Necview
{
    /// <summary>
    /// This module parses the NEC input file (i.e., the file describing the antenna's structure).
    /// </summary>
    public class ReadNecInput
    {
        /// <summary>
        /// some error codes for internal use
        /// </summary>
        private enum ERR
        {
            Err_none = 0,
            Err_misc = 1,
            Err_scan = 2,
            Err_nwires = 3,
            Err_nsurfaces = 4,
            Err_nexcis = 5,
            Err_nloads = 6,
            Err_nnetws = 7
        }

        /// <summary>
        /// next absolute segment number to be assigned
        /// </summary>
        private int nextabs;
        /// <summary>
        /// 
        /// </summary>
        
        /// <summary>
        /// type of previous card
        /// </summary>
        private char[] lastCard = { ' ', ' ' };

        /// <summary>
        /// read the input file, and update the arrays of wires and surfaces accordingly.
        /// </summary>
        /// <param name="filename">Name of file to read from.</param>
        /// <returns>returns true if succesful</returns>
        public Boolean Read_nec_input(string filename)
        {
            Debug.WriteLine(String.Format("Read_nec_input(string = {0}))", filename));
            System.IO.StreamReader file = new System.IO.StreamReader(filename);
            return Read_nec_input(file);
        }

        /// <summary>
        /// read the input file, and update the arrays of wires and surfaces accordingly.
        /// </summary>
        /// <param name="file">StreamReader to read from.</param>
        /// <returns>returns true if succesful</returns>
        public Boolean Read_nec_input(System.IO.StreamReader file)
        {
            string s;
            ERR err = ERR.Err_none;
            bool last_was_NT_or_TL = false;

            Debug.WriteLine(String.Format("Read_nec_input(StreamReader = {0})", file));

            nextabs = 1;

            // TODO stop on read of EN card
            while ((s = file.ReadLine()) != null)
            {
                err = ERR.Err_none;

                switch (s.Substring(0, 2).ToUpper())
                {
                    case "CM": err = Read_nec_CM(s); break;
                    case "EX": err = Read_nec_EX(s); break;
                    case "GA": err = Read_nec_GA(s); break;
                    //TODO implement parsing of GC cards
                    case "GC": break; // just ignore GC cards: they only modify the radius of GW cards, which xnecview doesn't use anyway
                    case "GH": err = Read_nec_GH(s); break;
                    case "GM": err = Read_nec_GM(s); break;
                    case "GR": err = Read_nec_GR(s); break;
                    case "GW": err = Read_nec_GW(s); break;
                    case "GX": err = Read_nec_GX(s); break;
                    case "GS": err = Read_nec_GS(s); break;
                    case "GF": System.Console.WriteLine("Warning: unhandled cardtype in {0}", s); break;
                    case "LD": err = Read_nec_LD(s); break;
                    case "MM": System.Console.WriteLine("{0}", s); break; // TODO what is a MM card?
                    case "ME": System.Console.WriteLine("{0}", s); break; // TODO what is a ME card?
                    case "NT":
                        if (!last_was_NT_or_TL)
                            Nec.netws.Clear();
                        err = Read_nec_NT(s);
                        last_was_NT_or_TL = true;
                        break;
                    case "SC": err = Read_nec_SC(s); break;
                    case "SM": err = Read_nec_SM(s); break;
                    case "SP": err = Read_nec_SP(s); break;
                    // TODO add support for 4nec2 SY cards
                    case "TL":
                        if (!last_was_NT_or_TL)
                            Nec.netws.Clear();
                        err = Read_nec_TL(s);
                        last_was_NT_or_TL = true;
                        break;
                    default:
                        System.Console.WriteLine("Warning: unhandled cardtype in {0}", s);
                        break;
                }
                switch (err)
                {
                    case ERR.Err_scan: System.Console.WriteLine("Error in parameters: {0}", s); break;
                    case ERR.Err_nwires: System.Console.WriteLine("Too many wires."); break;
                    case ERR.Err_nsurfaces: System.Console.WriteLine("Too many surfaces."); break;
                    case ERR.Err_nexcis: System.Console.WriteLine("Too many excitations."); break;
                    case ERR.Err_nloads: System.Console.WriteLine("Too many loads."); break;
                    case ERR.Err_nnetws: System.Console.WriteLine("Too many networks and/or transmission lines."); break;
                }

                lastCard[0] = s[0];
                lastCard[1] = s[1];
            }
            Nec.extr = Nec.extr_str;

#if (DEBUG)
            Nec.DumpWires();
            Nec.DumpSurfaces();
            Nec.DumpExci();
            Nec.DumpNetws();
#endif
            return true;
        }

        /// <summary>
        ///  calculates length of 3D vector
        /// </summary>
        /// <param name="p1"></param>
        /// <returns>Length of vector</returns>
        double Len3d(Nec.Point p1)
        {
            return Math.Sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
        }

        /// <summary>
        /// calculate intermediate point(3D);  f=0..1 -> p1..p2
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        Nec.Point Int3d(Nec.Point p1, Nec.Point p2, double f)
        {
            Nec.Point p = new Nec.Point();
            p.x = (float)(((1 - f) * p1.x) + (f * p2.x));
            p.y = (float)(((1 - f) * p1.y) + (f * p2.y));
            p.z = (float)(((1 - f) * p1.z) + (f * p2.z));
            return p;
        }

        /// <summary>
        /// calculates su->pa, with f determining the length of the arrow
        /// </summary>
        /// <param name="su"></param>
        /// <param name="f"></param>
        void Out3d(ref Nec.Surface su, double f)
        {
            Nec.Point d0 = new Nec.Point();
            Nec.Point d1 = new Nec.Point();
            Nec.Point no = new Nec.Point();

            d0.x = su.p0.x - su.pc.x; d0.y = su.p0.y - su.pc.y; d0.z = su.p0.z - su.pc.z;
            d1.x = su.p1.x - su.pc.x; d1.y = su.p1.y - su.pc.y; d1.z = su.p1.z - su.pc.z;
            no.x = d0.y * d1.z - d0.z * d1.y; no.y = -d0.x * d1.z + d0.z * d1.x; no.z = d0.x * d1.y - d0.y * d1.x;
            f *= (Len3d(d0) + Len3d(d1)) / Len3d(no);
            no.x *= f; no.y *= f; no.z *= f;
            su.pa.x = su.pc.x + no.x; su.pa.y = su.pc.y + no.y; su.pa.z = su.pc.z + no.z;
        }

        /// <summary>
        /// update the value of extr_str with information about this Point
        /// </summary>
        /// <param name="p"></param>
        void Updateextremes(Nec.Point p)
        {
            if (p.x > Nec.extr_str) Nec.extr_str = p.x;
            if (p.y > Nec.extr_str) Nec.extr_str = p.y;
            if (p.z > Nec.extr_str) Nec.extr_str = p.z;
            if (-p.x > Nec.extr_str) Nec.extr_str = -p.x;
            if (-p.y > Nec.extr_str) Nec.extr_str = -p.y;
            if (-p.z > Nec.extr_str) Nec.extr_str = -p.z;
        }

        /// <summary>
        /// perform a reflection of all antenna elements
        /// </summary>
        /// <param name="plane"></param>
        /// <param name="inc"></param>
        /// <returns></returns>
        int Reflect(int plane, int inc)
        {
            //            Wire* w,*w0;
            //            Surface* s,*s0;
            //            int i;

            //            EXPAND_IF_NECESSARY(numwires * 2, maxwires, wires)
            //           for (i = numwires, w0 = wires, w = wires + numwires; i > 0; i--, w0++, w++, numwires++)
            //            {
            //                *w = *w0;
            //                switch (plane)
            //                {
            //                    case 0: w->p0.x = -w->p0.x; w->p1.x = -w->p1.x; break;
            //                    case 1: w->p0.y = -w->p0.y; w->p1.y = -w->p1.y; break;
            //                    case 2: w->p0.z = -w->p0.z; w->p1.z = -w->p1.z; break;
            //                }
            //                if (w->itg > 0) w->itg += inc;
            //                w->abs = nextabs; nextabs += w->ns;
            //            }
            //            EXPAND_IF_NECESSARY(numsurfaces * 2, maxsurfaces, surfaces)
            //           for (i = numsurfaces, s0 = surfaces, s = surfaces + numsurfaces; i > 0; i--, s0++, s++, numsurfaces++)
            //            {
            //                *s = *s0;
            //                switch (plane)
            //                {
            //                    case 0: s->p0.x = -s->p0.x; s->p1.x = -s->p1.x; s->p2.x = -s->p2.x; s->p3.x = -s->p3.x; s->pc.x = -s->pc.x; s->pa.x = -s->pa.x; break;
            //                    case 1: s->p0.y = -s->p0.y; s->p1.y = -s->p1.y; s->p2.y = -s->p2.y; s->p3.y = -s->p3.y; s->pc.y = -s->pc.y; s->pa.y = -s->pa.y; break;
            //                    case 2: s->p0.z = -s->p0.z; s->p1.z = -s->p1.z; s->p2.z = -s->p2.z; s->p3.z = -s->p3.z; s->pc.z = -s->pc.z; s->pa.z = -s->pa.z; break;
            //                }
            //            }
            return 0;
        }

        /// <summary>
        /// change commas in string into spaces
        /// </summary>
        /// <param name="p"></param>
        private void Removecommas(ref string p)
        {
            string s = p.Replace(',', ' ');
            p = s;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tag"></param>
        /// <param name="seg"></param>
        /// <param name="card"></param>
        /// <param name="wirep"></param>
        /// <param name="segp"></param>
        /// <returns></returns>
        ERR Find_segment(int tag, int seg, string card, ref Nec.Wire wirep, ref int segp)
        {
            //  usage: err = find_segment(tag, seg, "EX", &excis[numexcis].wire, &excis[numexcis].seg);

            int i;

            if (tag > 0)
            {  /* location specified as tag,seg pair */
                for (i = 0; i < Nec.wires.Count; i++)
                    if (Nec.wires[i].itg == tag) break;
                if (i == Nec.wires.Count)
                {
                    Console.WriteLine("Unknown tag number on {0} card: {1}", card, tag);
                    return ERR.Err_misc;
                }
                if (Nec.wires[i].ns > 0)
                {
                    if (seg <= 0 || seg > Nec.wires[i].ns)
                    {
                        Console.WriteLine("Incorrect segment number on {0} card: {1}", card, seg);
                        return ERR.Err_misc;
                    }
                }
                else
                {
                    while (Nec.wires[i].ns != -seg && Nec.wires[i].itg == tag) i++;
                    if (Nec.wires[i].itg != tag)
                    {
                        Console.WriteLine("Incorrect segment number on {0} card: {1}", card, seg);
                        return ERR.Err_misc;
                    }
                }
                //*segp = seg;
            }
            else
            {                       /* location specified as absolute segment number */
                for (i = 0; i < Nec.wires.Count; i++)
                    if (seg >= Nec.wires[i].abs && seg < Nec.wires[i].abs + Nec.wires[i].ns) break;
                if (i == Nec.wires.Count)
                {
                    Console.WriteLine("Incorrect absolute segment number on {0} card: {1}", card, seg);
                    return ERR.Err_misc;
                }
                //*segp = seg - wires[i].abs + 1;
            }
            //*wirep = wires + i;
            return ERR.Err_none;
        }

        /// <summary>
        /// CM -> comment card; check whether it contains xnecview options
        /// </summary>
        /// <param name="s"></param>
        /// <returns>ERR enum</returns>
        private ERR Read_nec_CM(string s)
        {
            Debug.WriteLine(String.Format("Read_nec_CM({0})", s));

            Regex r_cm = new Regex("^CM\\s*xnecview:\\s*(?<options>.*)", RegexOptions.CultureInvariant | RegexOptions.Compiled);

            Match m = r_cm.Match(s);
            if (m.Success)
                Program.Process_optionstring(m.Result("${options}"));

            return ERR.Err_none;
        }

        /// <summary>
        /// EX -> excitation
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_EX(string s)
        {
            int type, tag, seg;
            ERR err = ERR.Err_none;

            Debug.WriteLine(String.Format("Read_nec_EX({0})", s));

            Removecommas(ref s);

            Regex r_ex = new Regex("^EX\\s*(?<type>\\d+)\\s*(?<tag>\\d+)\\s*(?<seg>\\d+)", RegexOptions.CultureInvariant | RegexOptions.Compiled);

            Match m = r_ex.Match(s);
            if (!m.Success)
                return ERR.Err_scan;

            type = int.Parse(m.Result("${type}"));
            tag = int.Parse(m.Result("${tag}"));
            seg = int.Parse(m.Result("${seg}"));

            if (type != 0 && type != 5)
            {
                Console.WriteLine("Only voltage sources can be displayed: {0}", s);
                return ERR.Err_misc;
            }

            //  err = find_segment(tag, seg, "EX", &excis[numexcis].wire, &excis[numexcis].seg);

            if (err != ERR.Err_none)
                return err;

            //            numexcis++;
            return ERR.Err_none;
        }

        /// <summary>
        /// TL -> transmission line
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_TL(string s)
        {
            int tag1, seg1, tag2, seg2;
            ERR err = ERR.Err_none;
            double imp;

            Debug.WriteLine(String.Format("Read_nec_TL({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "TL%d%d%d%d%lg", &tag1, &seg1, &tag2, &seg2, &imp) != 5) return ERR.Err_scan;

            //           err = find_segment(tag1, seg1, "TL", &netws[Nec.netws.Count()].wire1, &netws[Nec.netws.Count()].seg1);
            if (err != ERR.Err_none) return err;
            //            err = find_segment(tag2, seg2, "TL", &netws[Nec.netws.Count()].wire2, &netws[Nec.netws.Count()].seg2);
            if (err != ERR.Err_none) return err;
            //            if (imp < 0) netws[Nec.netws.Count()].type = -1; else netws[Nec.netws.Count()].type = 1;
            //            Nec.netws.Count()++;
            return ERR.Err_none;
        }

        /// <summary>
        /// NT -> network
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_NT(string s)
        {
            int tag1, seg1, tag2, seg2;
            ERR err = ERR.Err_none;

            Debug.WriteLine(String.Format("Read_nec_NT({0})", s));
            Removecommas(ref s);
            Debug.WriteLine(String.Format("no commas Read_nec_NT({0})", s));

            //            if (sscanf(s, "NT%d%d%d%d", &tag1, &seg1, &tag2, &seg2) != 4) return Err_scan;

            //           err = find_segment(tag1, seg1, "NT", &netws[Nec.netws.Count()].wire1, &netws[Nec.netws.Count()].seg1);
            //            if (err) return err;
            //            err = find_segment(tag2, seg2, "NT", &netws[Nec.netws.Count()].wire2, &netws[Nec.netws.Count()].seg2);
            if (err != ERR.Err_none) return err;
            //            netws[Nec.netws.Count()].type = 0;
            //            Nec.netws.Count()++;
            return ERR.Err_none;
        }

        /// <summary>
        /// add one load
        /// </summary>
        /// <param name="wi"></param>
        /// <param name="first"></param>
        /// <param name="last"></param>
        /// <returns></returns>
        int Addload(ref Nec.Wire wi, int first, int last)
        {
            //            EXPAND_IF_NECESSARY(numloads, maxloads, loads)
            //           loads[numloads].wire = wi;
            //            loads[numloads].firstseg = first;
            //            loads[numloads].lastseg = last;
            //            numloads++;
            return 0;
        }

        /// <summary>
        /// LD -> loading
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_LD(string s)
        {
            //            int type, tag, segf, segl;
            //            int i;
            //            Wire* wi;

            Debug.WriteLine(String.Format("Read_nec_LD({0})", s));

            Removecommas(ref s);
            //            tag = -2;
            //            segf = segl = 0;
            //            if (sscanf(s, "LD%d%d%d%d", &type, &tag, &segf, &segl) < 1) return Err_scan;
            //            if (tag == -2) return Err_scan;

            //            if (type == -1) { numloads = 0; return 0; }

            //            if (tag > 0)
            //            {                   /* location specified as tag,seg pair */
            //                for (i = 0; i < numwires; i++) if (wires[i].itg == tag) break;
            //                if (i == numwires) { fprintf(stderr, "Unknown tag number on LD card: %i\n", tag); return Err_misc; }
            //                if (wires[i].ns > 0)
            //                {
            //                    if (segf == 0) segf = 1;
            //                    if (segl == 0) segl = wires[i].ns;
            //                    if (segf <= 0 || segl > wires[i].ns || segf > segl) { fprintf(stderr, "Incorrect segment numbers on LD card: %i, %i\n", segf, segl); return Err_misc; }
            //                    addload(wires + i, segf, segl);
            //                    return 0;
            //                }
            //                else
            //                {
            //                    if (segf == 0) segf = -1;
            //                    else segf = -segf;
            //                    if (segl == 0) segl = -INT_MAX;
            //                    else segl = -segl;
            //                    while (wires[i].ns != segf && wires[i].itg == tag && i < numwires) i++;
            //                    if (wires[i].itg != tag || i >= numwires) { fprintf(stderr, "Incorrect first segment number on LD card: %i\n", segf); return Err_misc; }
            //                    while (wires[i].ns >= segl && wires[i].itg == tag && i < numwires)
            //                    {
            //                        addload(wires + i, 1, 1);
            //                        i++;
            //                    }
            //                    if ((wires[i].itg != tag || i >= numwires) && segl != -INT_MAX) { fprintf(stderr, "Incorrect last segment number on LD card: %i\n", -segl); return Err_misc; }
            //                }
            //            }
            //            else
            //            {                       /* location specified as absolute segment number */
            //                if (segf == 0 && segl == 0)
            //                {
            //                    for (i = 0; i < numwires; i++)  /* load all wires */
            //                        if (wires[i].ns > 0) addload(wires + i, 1, wires[i].ns);
            //                        else addload(wires + i, 1, 1);
            //                    return 0;
            //                }
            //                if (segf == 0 || segl == 0 || segf > segl) { fprintf(stderr, "Incorrect absolute segment numbers on LD card: %i %i\n", segf, segl); return Err_misc; }
            //                for (wi = wires; wi < wires + numwires; wi++)
            //                {
            //                    if (wi->ns > 0)
            //                    {
            //                        int f, l;
            //                        f = segf - wi->abs + 1;
            //                        if (f > wi->ns) continue;
            //                        l = segl - wi->abs + 1;
            //                        if (l < 1) continue;
            //                        if (f < 1) f = 1;
            //                        if (l > wi->ns) l = wi->ns;
            //                        addload(wi, f, l);
            //                    }
            //                    else
            //                    {
            //                        if (wi->abs >= segf && wi->abs <= segl) addload(wi, 1, 1);
            //                    }
            //                }
            //            }
            return ERR.Err_none;
        }

        /// <summary>
        /// GA -> wire arc
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GA(string s)
        {
            //            int ns, itg, segnr;
            //            double rada, ang1, ang2, rad;
            //            double phi, incphi;
            //            Point po;
            //            Wire* w;

            Debug.WriteLine(String.Format("Read_nec_GA({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "GA%d%d%lg%lg%lg%lg", &itg, &ns, &rada, &ang1, &ang2, &rad) != 6) return Err_scan;
            //            EXPAND_IF_NECESSARY(numwires + ns, maxwires, wires)
            //           phi = ang1 * M_PI / 180;
            //            incphi = (ang2 - ang1) * M_PI / 180 / ns;
            //            po.x = cos(phi) * rada; po.y = 0; po.z = sin(phi) * rada; updateextremes(&po);
            //            w = wires + numwires; numwires += ns;
            //            for (segnr = 1; segnr <= ns; segnr++)
            //            {
            //                w->p0 = po;
            //                phi += incphi;
            //                po.x = cos(phi) * rada; po.y = 0; po.z = sin(phi) * rada; updateextremes(&po);
            //                w->p1 = po;
            //                w->rad = rad;
            //                w->itg = itg;
            //                w->ns = -segnr;
            //                w->abs = nextabs; nextabs++;
            //                w->dispnr = (segnr - 1 == ns / 2);
            //                w++;
            //            }
            return ERR.Err_none;
        }

        /// <summary>
        /// GH -> helix
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GH(string s)
        {
            //            int ns, itg, segnr;
            //            double sp, hl, a1, b1, a2, b2, rad;
            //            double z, incz;
            //            double phi, incphi;
            //            double rx, incrx;
            //            double ry, incry;
            //            Point po;
            //            Wire* w;

            Debug.WriteLine(String.Format("Read_nec_GH({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "GH%d%d%lg%lg%lg%lg%lg%lg%lg", &itg, &ns, &sp, &hl, &a1, &b1, &a2, &b2, &rad) != 9) return Err_scan;
            //            EXPAND_IF_NECESSARY(numwires + ns, maxwires, wires)
            //           if (hl == 0) { fprintf(stderr, "Helix with length=0 not yet supported.\n"); return Err_misc; }
            //            z = 0; incz = hl / ns;
            //            rx = a1; incrx = (a2 - a1) / ns;
            //            ry = b1; incry = (b2 - b1) / ns;
            //            phi = 0;
            //            incphi = 2 * M_PI * incz / sp;
            //            if (hl < 0)
            //            {
            //                incz = -incz;
            //                phi = M_PI / 2;
            //            }
            //            po.x = cos(phi) * rx; po.y = sin(phi) * ry; po.z = z; updateextremes(&po);
            //            w = wires + numwires; numwires += ns;
            //            for (segnr = 1; segnr <= ns; segnr++)
            //            {
            //                w->p0 = po;
            //                z += incz; rx += incrx; ry += incry; phi += incphi;
            //                po.x = cos(phi) * rx; po.y = sin(phi) * ry; po.z = z; updateextremes(&po);
            //                w->p1 = po;
            //                w->rad = rad;
            //                w->itg = itg;
            //                w->ns = -segnr;
            //                w->abs = nextabs; nextabs++;
            //                w->dispnr = (segnr - 1 == ns / 2);
            //                w++;
            //            }
            return ERR.Err_none;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="p0"></param>
        /// <param name="p1"></param>
        /// <param name="mat"></param>
        /// <param name="xs"></param>
        /// <param name="ys"></param>
        /// <param name="zs"></param>
        void Pointtransf(ref Nec.Point p0, ref Nec.Point p1, double[] mat, double xs, double ys, double zs)
        {
            Nec.Point po = new Nec.Point();

            po.x = xs + mat[0] * p0.x + mat[1] * p0.y + mat[2] * p0.z;
            po.y = ys + mat[3] * p0.x + mat[4] * p0.y + mat[5] * p0.z;
            po.z = zs + mat[6] * p0.x + mat[7] * p0.y + mat[8] * p0.z;
            p1 = po;
            Updateextremes(po);
        }

        /// <summary>
        /// GM -> coordinate transformation
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GM(string s)
        {
            int nrpt, its, i, itsi;
            double rox, roy, roz, xs, ys, zs;
            double[] mat = new double[9];
            double[] mmat = new double[9];
            double co, si;
            //            Wire* w00,*w0,*w1,*wlast;
            //            Surface* s0,*s1;
            double its_double;

            Debug.WriteLine(String.Format("Read_nec_GM({0})", s));

            Removecommas(ref s);
            rox = roy = roz = xs = ys = zs = its = 0;
            //            i = sscanf(s, "GM%d%d%lg%lg%lg%lg%lg%lg%lg", &itsi, &nrpt, &rox, &roy, &roz, &xs, &ys, &zs, &its_double);
            //            if (i < 2) return Err_scan;

            //            its = (int)(its_double + 0.5);   // round to double; the detour via double is needed because some programs write something like 1.1e1 instead of 11 in this field, and this is not forbidden according to the specification

            rox *= Math.PI / 180.0;
            roy *= Math.PI / 180.0;
            roz *= Math.PI / 180.0;

            mat[0] = mat[4] = mat[8] = 1;
            mat[1] = mat[2] = mat[3] = 0;
            mat[5] = mat[6] = mat[7] = 0;

            co = Math.Cos(rox);
            si = Math.Sin(rox);
            //            memcpy(mmat, mat, sizeof(mat));
            mat[3] = co * mmat[3] - si * mmat[6];
            mat[4] = co * mmat[4] - si * mmat[7];
            mat[5] = co * mmat[5] - si * mmat[8];
            mat[6] = si * mmat[3] + co * mmat[6];
            mat[7] = si * mmat[4] + co * mmat[7];
            mat[8] = si * mmat[5] + co * mmat[8];

            co = Math.Cos(roy);
            si = -Math.Sin(roy);
            //            memcpy(mmat, mat, sizeof(mat));
            mat[0] = co * mmat[0] - si * mmat[6];
            mat[1] = co * mmat[1] - si * mmat[7];
            mat[2] = co * mmat[2] - si * mmat[8];
            mat[6] = si * mmat[0] + co * mmat[6];
            mat[7] = si * mmat[1] + co * mmat[7];
            mat[8] = si * mmat[2] + co * mmat[8];

            co = Math.Cos(roz);
            si = Math.Sin(roz);
            //            memcpy(mmat, mat, sizeof(mat));
            mat[0] = co * mmat[0] - si * mmat[3];
            mat[1] = co * mmat[1] - si * mmat[4];
            mat[2] = co * mmat[2] - si * mmat[5];
            mat[3] = si * mmat[0] + co * mmat[3];
            mat[4] = si * mmat[1] + co * mmat[4];
            mat[5] = si * mmat[2] + co * mmat[5];

            //            EXPAND_IF_NECESSARY(numwires * (nrpt + 1), maxwires, wires)
            //           if (nrpt == 0)
            //            {
            //                w00 = w1 = wires;
            //            }
            //            else
            //            {
            //                w00 = wires;
            //                w1 = wires + numwires;
            //            }
            //            wlast = wires + numwires;
            //            while (w00 < wlast)
            //            {
            //                if (w00->itg >= its)
            //                {
            //                    w0 = w00;
            //                    i = 0;
            //                    do
            //                    {
            //                        if (nrpt > 0)
            //                        {
            //                            numwires++;
            //                        }
            //                        else
            //                        {
            //                            w1 = w0;
            //                        }
            //                        pointtransf(&w0->p0, &w1->p0, mat, xs, ys, zs);
            //                        pointtransf(&w0->p1, &w1->p1, mat, xs, ys, zs);
            //                        w1->rad = w0->rad;
            //                        w1->ns = w0->ns;
            //                        if (w0->itg > 0) w1->itg = w0->itg + itsi;
            //                        else w1->itg = w0->itg;
            //                        w1->abs = nextabs; nextabs += w1->ns;
            //                        w1->dispnr = w0->dispnr;
            //                        w0 = w1;
            //                        w1++;
            //                    } while (++i < nrpt);
            //                }
            //                w00++;
            //            }

            //            if (nrpt == 0)
            //            {
            //                s0 = s1 = surfaces;
            //                i = numsurfaces;
            //            }
            //            else
            //            {
            //                int oldnumsurfaces = numsurfaces;
            //                i = numsurfaces * nrpt;
            //                numsurfaces *= (1 + nrpt);
            //                EXPAND_IF_NECESSARY(numsurfaces, maxsurfaces, surfaces)
            //             s0 = surfaces;
            //                s1 = surfaces + oldnumsurfaces;
            //            }
            //            while (i--)
            //            {
            //                s1->type = s0->type;
            //                pointtransf(&s0->p0, &s1->p0, mat, xs, ys, zs);
            //                pointtransf(&s0->p1, &s1->p1, mat, xs, ys, zs);
            //                pointtransf(&s0->p2, &s1->p2, mat, xs, ys, zs);
            //                pointtransf(&s0->p3, &s1->p3, mat, xs, ys, zs);
            //                pointtransf(&s0->pc, &s1->pc, mat, xs, ys, zs);
            //                pointtransf(&s0->pa, &s1->pa, mat, xs, ys, zs);
            //                s0++; s1++;
            //            }

            return ERR.Err_none;
        }

        /// <summary>
        /// GR -> generate cylindrical structure
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GR(string s)
        {
            char[] gm = new char[60];
            int n, inc;

            Debug.WriteLine(String.Format("Read_nec_GR({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "GR%d%d", &inc, &n) != 2) return Err_scan;
            //            sprintf(gm, "GM %i %i 0 0 %g 0 0 0", inc, n - 1, 360./ n);
            //            read_nec_GM(gm);           /* GR card is translated into equivalent GM card */
            return ERR.Err_none;
        }

        /// <summary>
        /// GS -> geometry scale
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GS(string s)
        {
            double a;
            //            Wire* w;
            //            Surface* su;
            //int i;

            Debug.WriteLine(String.Format("Read_nec_GS({0})", s));

            Removecommas(ref s);

            Regex r_gs = new Regex("^GS\\s(?<unused0>[-+]?[0-9]*\\.?[0-9]+)\\s(?<unused1>[-+]?[0-9]*\\.?[0-9]+)\\s(?<f1>[-+]?[0-9]*\\.?[0-9]+)", RegexOptions.CultureInvariant | RegexOptions.Compiled);

            Match m = r_gs.Match(s);
            if (!m.Success)
                return ERR.Err_scan;

            a = Double.Parse(m.Result("${f1}"));
            Console.WriteLine("GS scale = {0}", a);

            foreach (Nec.Wire w in Nec.wires)
            {
                w.p0.x *= a;
                w.p0.y *= a;
                w.p0.z *= a;
                w.p1.x *= a;
                w.p1.y *= a;
                w.p1.z *= a;
                w.rad *= a;
            }
            //            for (i = 0, su = surfaces; i < numsurfaces; i++, su++)
            //            {
            //                su->p0.x *= a; su->p0.y *= a; su->p0.z *= a;
            //                su->p1.x *= a; su->p1.y *= a; su->p1.z *= a;
            //                su->p2.x *= a; su->p2.y *= a; su->p2.z *= a;
            //                su->p3.x *= a; su->p3.y *= a; su->p3.z *= a;
            //                su->pc.x *= a; su->pc.y *= a; su->pc.z *= a;
            //                su->pa.x *= a; su->pa.y *= a; su->pa.z *= a;
            //            }
            Nec.extr_str *= a;
            return ERR.Err_none;
        }

        /// <summary>
        /// GW -> straight wire
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GW(string s)
        {
            Regex r_gw9 = new Regex(    // expression for with the wire radius
              "^GW\\s(?<itg>[-+]?[0-9]*\\.?[0-9]+)\\s(?<ns>[-+]?[0-9]*\\.?[0" +
              "-9]+)\\s(?<p0x>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p0y>[-+]?[0-9]*\\.?" +
              "[0-9]+)\\s(?<p0z>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p1x>[-+]?[0-9]*\\." +
              "?[0-9]+)\\s(?<p1y>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p1z>[-+]?[0-9]*\\." +
              "?[0-9]+([eE][-+]?[0-9]+)?)\\s(?<rad>[-+]?[0-9]*\\.?[0-9]+([eE" +
              "][-+]?[0-9]+)?)",
            RegexOptions.CultureInvariant | RegexOptions.Compiled
            );
            Regex r_gw8 = new Regex(    // expression for without the wire radius
              "^GW\\s(?<itg>[-+]?[0-9]*\\.?[0-9]+)\\s(?<ns>[-+]?[0-9]*\\.?[0" +
              "-9]+)\\s(?<p0x>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p0y>[-+]?[0-9]*\\.?" +
              "[0-9]+)\\s(?<p0z>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p1x>[-+]?[0-9]*\\." +
              "?[0-9]+)\\s(?<p1y>[-+]?[0-9]*\\.?[0-9]+)\\s(?<p1z>[-+]?[0-9]*\\." +
              "?[0-9]+([eE][-+]?[0-9]+)?)",
            RegexOptions.CultureInvariant | RegexOptions.Compiled
            );
            Nec.Wire w = new Nec.Wire();

            Debug.WriteLine(String.Format("Read_nec_GW({0})", s));

            Removecommas(ref s);

            Match m = r_gw9.Match(s);
            if (m.Success)
            {
                w.rad = Double.Parse(m.Result("${rad}"));
                // the rest of the fields are read below
            }
            else
            {
                m = r_gw8.Match(s);
                // omitting the radius seems to be allowed, for tapered wires
                w.rad = 0;
            }

            if (!m.Success)
                return ERR.Err_scan;

            w.itg = int.Parse(m.Result("${itg}"));
            w.ns = int.Parse(m.Result("${ns}"));
            w.p0.x = int.Parse(m.Result("${p0x}"));
            w.p0.y = int.Parse(m.Result("${p0y}"));
            w.p0.z = int.Parse(m.Result("${p0z}"));
            w.p1.x = int.Parse(m.Result("${p1x}"));
            w.p1.y = int.Parse(m.Result("${p1y}"));
            w.p1.z = int.Parse(m.Result("${p1z}"));

            w.abs = nextabs;
            nextabs += w.ns;
            w.dispnr = true;

            Console.WriteLine("GW rad={0} abs={1}", w.rad, w.abs);

            Updateextremes(w.p0);
            Updateextremes(w.p1);
            Nec.wires.Add(w);
            return ERR.Err_none;
        }

        /// <summary>
        /// GX -> reflection in coordinate planes
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_GX(string s)
        {
            //            int i, inc, err;

            Debug.WriteLine(String.Format("Read_nec_GX({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "GX%d%d", &inc, &i) != 2) return Err_scan;
            //            if (i % 10 == 1) { err = reflect(2, inc); if (err) return err; inc *= 2; }
            //            if ((i / 10) % 10 == 1) { err = reflect(1, inc); if (err) return err; inc *= 2; }
            //            if ((i / 100) == 1) return reflect(0, inc);
            return ERR.Err_none;
        }

        /// <summary>
        /// SC -> continuation of SM, SP, or SC
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_SC(string s)
        {
            //            Surface* su;
            //            int ns;
            //            Point p3;
            //            int n;

            //            if (lastCard[0] != 'S') { fprintf(stderr, "SC card not preceded by SM, SP, or SC.\n"); return Err_misc; }

            //            EXPAND_IF_NECESSARY(numsurfaces, maxsurfaces, surfaces)
            //           su = surfaces + numsurfaces++;

            Debug.WriteLine(String.Format("Read_nec_SC({0})", s));

            Removecommas(ref s);

            //            switch (lastCard[1])
            //            {
            //                case 'M':
            //                    if (sscanf(s, "SC%*i%*i%g%g%g", &su->p2.x, &su->p2.y, &su->p2.z) != 3) return Err_scan;
            //                    updateextremes(&su->p2);
            //                    su->pc = int3d(su->p0, su->p2, 0.5);
            //                    p3 = int3d(su->p1, su->pc, 2);
            //                    updateextremes(&p3);
            //                    out3d(su, 0.05);
            //                    su->type = SU_rect;
            //                    break;
            //                case 'P':
            //                    n = sscanf(s, "SC%*i%*i%g%g%g%g%g%g", &su->p2.x, &su->p2.y, &su->p2.z, &su->p3.x, &su->p3.y, &su->p3.z);
            //                    if (n != 6 && n != 3) return Err_scan;
            //                    updateextremes(&su->p2);
            //                    switch (su->type)
            //                    {
            //                        case SU_rect:
            //                            su->pc = int3d(su->p0, su->p2, 0.5);
            //                            su->p3 = int3d(su->p1, su->pc, 2);    // we don't really need su->p3 now, but do need it in case of continuation using more SC cards
            //                            updateextremes(&su->p3);
            //                            out3d(su, 0.05);
            //                            break;
            //                        case SU_tri:
            //                            su->pc = int3d(int3d(su->p0, su->p1, 0.5), su->p2, 1./ 3);
            //                            out3d(su, 0.05);
            //                            break;
            //                        case SU_quad:
            //                            if (n != 6) return Err_scan;
            //                            updateextremes(&su->p3);
            //                            su->pc = int3d(int3d(su->p0, su->p1, 0.5), int3d(su->p2, su->p3, 0.5), 0.5);
            //                            out3d(su, 0.05);
            //                            break;
            //                        default:
            //                            fprintf(stderr, "SC card for surface type which is not rectangle/triangle/quadrilateral.\n"); return Err_misc;
            //                    }
            //                    break;
            //                case 'C':
            //                    n = sscanf(s, "SC%*i%d%g%g%g%g%g%g", &ns, &su->p2.x, &su->p2.y, &su->p2.z, &su->p3.x, &su->p3.y, &su->p3.z);
            //                    if (n != 7 && n != 4) return Err_scan;
            //                    Surface* suLast = surfaces + (numsurfaces - 2);
            //                    su->p0 = suLast->p3;
            //                    su->p1 = suLast->p2;
            //                    updateextremes(&su->p2);
            //                    switch (ns)
            //                    {
            //                        case 1:
            //                            su->type = SU_rect;
            //                            su->pc = int3d(su->p0, su->p2, 0.5);
            //                            su->p3 = int3d(su->p1, su->pc, 2);
            //                            updateextremes(&su->p3);
            //                            out3d(su, 0.05);
            //                            break;
            //                        case 3:
            //                            if (n != 7) return Err_scan;
            //                            su->type = SU_quad;
            //                            updateextremes(&su->p3);
            //                            su->pc = int3d(int3d(su->p0, su->p1, 0.5), int3d(su->p2, su->p3, 0.5), 0.5);
            //                            out3d(su, 0.05);
            //                            break;
            //                        default:
            //                            fprintf(stderr, "SC card specifying unsupported shape NS=%d.\n", ns); return Err_misc;
            //                    }
            //                    break;
            //                default:
            //                    fprintf(stderr, "SC card following card %c%c.\n", lastCard[0], lastCard[1]); return Err_misc;
            //            }

            return ERR.Err_none;
        }

        /// <summary>
        /// SM -> multiple surface patch (shown as one large rectangular piece of surface)
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_SM(string s)
        {
            //            Surface* su;
            //            EXPAND_IF_NECESSARY(numsurfaces, maxsurfaces, surfaces)
            //           su = surfaces + numsurfaces;

            Debug.WriteLine(String.Format("Read_nec_SM({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "SM%*i%*i%g%g%g%g%g%g", &su->p0.x, &su->p0.y, &su->p0.z, &su->p1.x, &su->p1.y, &su->p1.z) != 6) return Err_scan;
            //            updateextremes(&su->p0);
            //            updateextremes(&su->p1);
            return ERR.Err_none;
        }

        /// <summary>
        /// SP -> surface patch
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        ERR Read_nec_SP(string s)
        {
            //            Surface* su;
            //            int ns;
            //            double x1, y1, z1, x2, y2, z2;
            //            double r;

            //            EXPAND_IF_NECESSARY(numsurfaces, maxsurfaces, surfaces)
            //           su = surfaces + numsurfaces;
            //            su->type = SU_arb;

            Debug.WriteLine(String.Format("Read_nec_SP({0})", s));

            Removecommas(ref s);
            //            if (sscanf(s, "SP%*d%d%lg%lg%lg%lg%lg%lg", &ns, &x1, &y1, &z1, &x2, &y2, &z2) != 7) return Err_scan;

            //            if (ns == 0)
            //            {
            //                su->type = SU_arb;
            //                su->pc.x = x1; su->pc.y = y1; su->pc.z = z1;
            //                updateextremes(&su->pc);
            //                r = sqrt(z2 / M_PI);
            //                x2 *= M_PI / 180; y2 *= M_PI / 180;
            //                su->p0.x = x1 - r * sin(y2);
            //                su->p0.y = y1 + r * cos(y2);
            //                su->p0.z = z1;
            //                su->p1.x = x1 - r * sin(x2) * cos(y2);
            //                su->p1.y = y1 - r * sin(x2) * sin(y2);
            //                su->p1.z = z1 + r * cos(x2);
            //                su->p2 = int3d(su->p0, su->pc, 2);
            //                out3d(su, 0.1);
            //                ++numsurfaces;
            //                return 0;
            //            }

            //            switch (ns)
            //            {
            //                case 1:
            //                    su->type = SU_rect;
            //                    break;
            //                case 2:
            //                    su->type = SU_tri;
            //                    break;
            //                case 3:
            //                    su->type = SU_quad;
            //                    break;
            //                default:
            //                    fprintf(stderr, "Unknown patch shape: NS=%i.\n", ns); return Err_misc;
            //            }
            //            su->p0.x = x1; su->p0.y = y1; su->p0.z = z1;
            //            su->p1.x = x2; su->p1.y = y2; su->p1.z = z2;
            //            updateextremes(&su->p0);
            //            updateextremes(&su->p1);
            return ERR.Err_none;
        }
    }
}
