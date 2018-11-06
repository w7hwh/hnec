/*  XNECVIEW - a program for visualizing NEC2 input and output data
 *
 *  Copyright (C) 2002-2005, Pieter-Tjerk de Boer -- pa3fwm@amsat.org
 *
 *  Distributed on the conditions of version 2 of the GPL: see the files
 *  README and COPYING, which accompany this source file.
 *
 *  This module exclusively contains code for drawing the gain pattern as
 *  an opaque surface, i.e., with hidden-line removal.
 *
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    class DrawOpaque
    {
        /*******************************************************************************

        The principle of the algorithm used here can be summarized as follows:
        a) draw parts of the gain surface in such an order, that for any bit of screen
           area, from all parts of the gain surface that appear on that bit of screen,
           the part closest to the observer is drawn first.
        b) keep track of which parts of the screen have already been covered by the
           parts of the gain surface drawn so far, and of any new part of the gain
           surface to be drawn, only draw those parts that are outside the parts of
           the screen already covered. This "covered area" is henceforth called "mask".
        Clearly, these two rules guarantee that only the visible parts of the gain
        surface are drawn.

        The above is deliberately vague about "parts of the gain surface". What are
        these? Well, when drawing a wire frame, we're simply drawing line segments.
        However, (infinitely thin) lines by themselves do not cover any area on the
        screen, so doing rule (b) would be hard on the basis of line segments.
        Actually, areas on the screen are covered by the quadrangles formed by those
        lines.
        So: we're drawing line segments, but we need to keep track of the screen area
        covered by the quadrangles formed by those line segments.

        How do we implement an ordering as described by rule (a) ?
        One way would be to sort the line segments or quadrangles according to their
        distance from the viewer, i.e., their position on an axis perpendicular to the
        screen. However, this may easily lead to rule (b) needing to keep track of
        multiple disjoint areas on the screen.

        Since the gain surface we're drawing is a single-valued function of the angular
        coordinates (phi and theta), a different way to satisfy rule (a) is sorting
        the line segments or quadrangles by the angle a line from the origin to them
        makes with a line from the origin to the viewer. Basically, start with the
        quadrangle that is between the viewer and the origin (i.e., the antenna gain
        in the direction toward the viewer), then continue "outward" from this initial
        quadrangle. One can easily verify that indeed rule (a) is satisfied.
        Furthermore, since we start with the quadrangle around the origin and work
        outward from there, it is also clear that the "covered area" that rule (b)
        is about, will be one single area.

        It would be even nicer (and efficient for the calculations), if the boundary
        of this single "covered area" (mask) would be a single-valued function of an
        angular coordinate on the screen. In other words: starting from the projection
        of the origin on the screen and moving outward on the screen along a straight
        line in any direction, one should cross the boundary of the mask exactly
        once.
        This can be achieved without violating rule (a) by working one's way through
        the quadrangles in the right order; some more details are in the comments
        in the function draw_opaque(), below.

        Note that the above algorithm bears quite a few similarities to the "obvious"
        way of plotting a function of two variables in a 3D Cartesian coordinate
        system: start with the front-most line, and keep track of the screen area
        already covered. See e.g. Section 15.1 in Foley/van Dam/Feiner/Hughes, 1996.
        In fact, the gain pattern is also a function of two variables, only the
        coordinate system is different.

        A completely different approach would be to draw the quadrangles starting
        with the ones farthest away, and drawing them as actually *filled* polygons.
        Thus, the front-most polygon will simply erase the ones farther away if
        they occupy the same screen position. I decided not to use this algorithm
        because of the large number of unnecessary drawing operations involved;
        I hope my algorithm is faster, although I didn't compare them.

        A note about coordinates used in this file:
        x,y  are generally the on-screen coordinates, except for a translation:
           the projection of the origin is at x=0, y=0. The offset with the true
           screen coordinates (in which the origin may be anywhere inside or
           outside the window) is done at the final call to DrawLine().
           These coordinates are used for actually drawing line segments, and for
           calculating intersections of new quadrangle boundaries with the mask.
        a    is defined as ATAN(y,x), with ATAN(y,x) as #define'd below. Basically,
           it is the angle (in radians) of a line (on screen) from the origin to the
           point x,y, with some reference line. With the present definition of ATAN:
            a=0    is y>0, x=0 (i.e. "south" on the screen)
            a=pi/2 is x<0, y=0 (i.e. "east" on the screen)
           etc.
           This coordinate is used for sorting the points that describe the polygon
           which is the mask.
        s    is defined as x*x+y*y . It is simply the square of the distance of x,y
           from the (on-screen) origin. 
           It is only used to speed up decisions about (in)visibility; e.g., if both
           endpoints of a new line segment have  s  values less than the smallest  s
           value presently in the mask, the new line segment surely is completely
           hidden by the mask.

        Finally: many of the complications in the implementation are caused by the
        possibility that points coincide exactly, or have exactly the same a value;
        particularly, this may happens with points on the Z axis.
        Of course, floating point arithmetic is not exact, so checking for equality
        is not a good idea. On the other hand, checking for "approximate equality"
        may take more time, and may also cause problems by itself.
        An additional problem is the fact that the Maskpoints are sorted by their
        a value (after all, the mask is a function of a), but there may be purely radial
        segments in the mask, so more sorting is needed than can be derived from
        the a values.
        The present implementation seems to work correctly in all examples I tried,
        but there might very well still be a case where it doesn't work correctly;
        if you run into it, please let me know.
        Also, perhaps on other architectures than i386 or with a different compiler,
        it might not work correctly, if the tests for equality of floats are handled
        differently. Again, please let me know if you encounter problems.

        *******************************************************************************/

        //#include <stdio.h>
        //#include <math.h>
        //#include <stdlib.h>
        //#include <string.h>

        //#include <sys/time.h>
        //#include <sys/resource.h>
        //#include <unistd.h>

        //#include <assert.h>

        //#include "xnecview.h"


        ///* #defining  DEBUG  prints lots of information about the drawing process,
        //   and regularly pauses the drawing waiting for a keypress: '>' does the
        //   next step, '<' first draws the present version of the mask.
        //*/
        //#undef DEBUG      

        ///* #defining  MEASURE  enables code that draw the picture 100 times and measures the time */
        //#undef MEASURE


        //#ifdef DEBUG
        // #define PRINTF(a)    printf a /**/
        // #define DEBUGPAUSE   \
        //   {   \
        //   double aa = animphase;   \
        //   printf("i=%i j=%i  theta=%f phi=%f \n", i, j, rp->gtheta[i]*180/M_PI, rp->gphi[j]*180/M_PI);   \
        //   out->Complete();   \
        //   do { gtk_main_iteration_do(1);
        //    } while (aa==animphase);   \
        //   if (animphase<aa) { draw_mask(); SetForeground(&c_gain);
        //}   \
        //   }
        //#else
        // #define PRINTF(a)    /**/
        // #define DEBUGPAUSE  /**/
        //#endif

        //#ifdef MEASURE
        // int count;
        //#define DrawLine(a,b,c,d) { if (count==0) out->DrawLine((a),(b),(c),(d)); }
        //#endif


        ///* ------------------------- separate versions of calcproj() and proj() -------------------------------------------------- */



        //double ATAN(double y, double x)
        //{
        //    double a;
        //    a = -atan2(x, y);
        //    if (a == M_PI) return -M_PI;
        //    if (fabs(a) < 1e-7) return 0;
        //    return a;
        //}



        ///* this data type is used to store all relevant information about a point in the gain pattern
        //   as projected onto the 2D-screen: */
        //typedef struct {
        //   double x, y;        /* x,y,a,s: see above description of the coordinates */
        //double a;
        //double s;
        //double c;          /* this variable is the axial ratio of the radiation in this direction, and is used (in the coloured polarization mode) for determining the colour of the line segments */
        //void* maskpoint;   /* pointer to the corresponding entry in the mask polygon, if this point presently is a corner point of the mask */
        //} Hidpoint;


        ///* we have a separate projection matrix, since we sometimes need to change the viewing direction very slightly */
        //double Xx2, Xy2, Yx2, Yy2, Yz2;         /* projection matrix */
        //double Xtr2, Ytr2;                   /* 2D-translation */
        //double phi2, theta2;                 /* in radians! */

        //static void inline calcproj2(void)     /* similar to calcproj(), but uses phi2 and theta2 */
        //{
        //    double scale;
        //    scale = winsize * zoom / (3.0 * extr);
        //    Yz2 = -scale * sin(theta2);
        //    Xx2 = -scale * sin(phi2);
        //    Xy2 = scale * cos(phi2);
        //    Yx2 = scale * cos(phi2) * cos(theta2);
        //    Yy2 = scale * sin(phi2) * cos(theta2);
        //    Xtr2 = 0.5 + trx * zoom * winsize + 0.5 * winsizex;
        //    Ytr2 = 0.5 + try*zoom * winsize + 0.5 * winsizey;
        //    }

        //static void inline proj2(Point * p, Hidpoint * hp)
        //{
        //        hp->x = Xx2 * p->x + Xy2 * p->y;
        //        hp->y = Yx2 * p->x + Yy2 * p->y + Yz2 * p->z;
        //        hp->a = ATAN(hp->y, hp->x);
        //        hp->s = hp->x * hp->x + hp->y * hp->y;
        //    }

        //    /* ------------------------- handling of the mask ----------------------------------------------------- */


        //    /* data structure for one corner of the polygon that is the mask: */
        //    typedef struct _Maskpoint
        //{
        //    double x, y;               /* x,y,a,s: see coordinate description above */
        //    double a;
        //    double s;
        //    struct _Maskpoint *next;  /* pointer to the next point in the mask (sorted on a) */
        //   short int del;            /* flag: if !=0, this maskpoint is to be deleted at the next commit_mask() call */
        //    short int after;          /* flag, only used in the *new list: if set, means this Maskpoint should be inserted after (rather than before) any already existing point(s) with the samen a value; -1 means unknown, let after-flag of other point decide */
        //    Hidpoint* hp;             /* pointer to the corresponding Hidpoint, if any */
        //}
        //Maskpoint;

        //Maskpoint* Q;         /* points to begin of linked list of Maskpoints, ordered by increasing a */
        //Maskpoint* Q2;        /* pointer into that list, pointing at last Maskpoint with a<0; used to speed up searching the list */
        //Maskpoint* new;       /* pointer to linked list of Maskpoints that need to be included into the list at the next commit_mask() call */
        //Maskpoint* spare;     /* pointer to linked list of Maskpoint buffers that are available for use; this prevents us from needing to call malloc() and free() too often */

        //double min_s = 0;       /* roughly: lowest value of s currently in the mask 
        //                         More precisely: min_s is never higher than the lowest s value corresponding to any point on the mask;
        //                         note that this lowest s value may well occur not on a Maskpoint, but somewhere in between;
        //                         in recalc_min_s(), this lowest s value is not caculated exactly, but a lower bound is calculated.
        //                         min_s is used to quickly detect lines that are well within the mask, so they can be skipped without further analysis.
        //                      */
        //Maskpoint* min_s_mp = NULL;  /* Maskpoint that had the lowest s value contributing to min_s;
        //                              if this Maskpoint is deleted from the list, min_s should be recalculated since it probably can increase
        //                           */

        //int left_half = 0, right_half = 0;  /* two flags, indicating whether we're at the moment drawing on the left half (x<0), right half (x>0) or possibly both. This min_s value above is actually calculated only for the active half(es) */

        ///* a handy macro for traversing the linked list as if head and tail were joined (after all, a is an angular coordinate): */
        //#define NEXT(p) ( (p)->next ? (p)->next : Q )


        //static void alloc_sparemaskpoints(void)
        //{
        //    int i;
        //    spare = mymalloc(100 * sizeof(Maskpoint));
        //    for (i = 0; i < 99; i++) spare[i].next = spare + i + 1;
        //    spare[99].next = NULL;
        //}

        //static void inline mpfree(Maskpoint* p)
        //{
        //    p->next = spare;
        //    spare = p;
        //}

        //static void cleanup_mask(void)
        //{
        //    Maskpoint* p;
        //    while (Q)
        //    {
        //        p = Q;
        //        Q = Q->next;
        //        mpfree(p);
        //    }
        //}

        //static Maskpoint* newmaskpoint(double x, double y, double a, double s, Hidpoint* hp)
        //{
        //    Maskpoint* mp;
        //    if (!spare) alloc_sparemaskpoints();
        //    mp = spare;
        //    spare = mp->next;
        //    PRINTF(("newmaskpoint %f %f %f %p\n", x, y, s, hp));
        //    mp->x = x;
        //    mp->y = y;
        //    mp->a = a;
        //    mp->s = s;
        //    mp->next = NULL;
        //    mp->del = 0;
        //    mp->after = 0;
        //    mp->hp = hp;
        //    return mp;
        //}

        //static void init_mask(void)
        //{
        //    if (!spare) alloc_sparemaskpoints();
        //    new= NULL;
        //    min_s = 0;
        //    left_half = right_half = 1;
        //}


        //static void add_to_mask(Maskpoint* q)
        //{
        //    Maskpoint* p;
        //    p = Q;
        //    PRINTF(("============ add_to_mask: %g %g\n", q->a, p->a));
        //    if (q->a >= 0) p = Q2;
        //    else if (q->a <= p->a)
        //    {
        //        if (fabs(p->x - q->x) < 1e-8 && fabs(p->y - q->y) < 1e-8)
        //        {
        //            mpfree(q);
        //            return;
        //        }
        //        if (q->a == p->a && q->after == -1) q->after = !p->after;
        //        if (q->a == p->a && q->after)
        //        {
        //            while (p->next->a == q->a) p = p->next;
        //            q->next = p->next;
        //            p->next = q;
        //            if (q->hp) q->hp->maskpoint = q;
        //            return;
        //        }
        //        q->next = p;
        //        Q = q;
        //        if (q->hp) q->hp->maskpoint = q;
        //        return;
        //    }
        //    while (p->next && p->next->a < q->a) p = p->next;
        //    if (!p->next)
        //    {
        //        p->next = q;
        //        q->next = NULL;
        //        if (q->hp) q->hp->maskpoint = q;
        //        return;
        //    }
        //    if (fabs(p->x - q->x) < 1e-8 && fabs(p->y - q->y) < 1e-8 && !p->del)
        //    {
        //        mpfree(q);
        //        return;
        //    }
        //    if (p->next && fabs(p->next->x - q->x) < 1e-8 && fabs(p->next->y - q->y) < 1e-8 && !p->next->del)
        //    {
        //        mpfree(q);
        //        return;
        //    }
        //    if (q->a == p->next->a && q->after == -1) q->after = !p->next->after;
        //    if (q->a == p->next->a && q->after)
        //    {
        //        while (p->next && p->next->a == q->a) p = p->next;
        //        q->next = p->next;
        //        p->next = q;
        //        if (q->hp) q->hp->maskpoint = q;
        //        return;
        //    }
        //    q->next = p->next;
        //    p->next = q;
        //    if (q->hp) q->hp->maskpoint = q;
        //    if (q->a < 0 && q->a > Q2->a) Q2 = q;
        //}


        //static void add_to_newlist(Maskpoint* q)
        //{
        //    if (new && new->hp && q->hp == new->hp)
        //    {
        //        PRINTF(("============ add_to_newlist: dropped (%p,%p)\n", q->hp, new->hp));
        //        mpfree(q);
        //        return;
        //    }
        //    q->next = new;
        //    new= q;
        //}

        //int find_in_mask(Hidpoint* h, Maskpoint** pp, Hidpoint* h1)
        ///* For a given point *h, search the corresponding entry in the mask.
        //   Return values:
        //    0 = *h lies inside mask; *pp points to maskpoint just before *h
        //    1 = *h is one of the mask's corner points, *pp points to this maskpoint, and *h1 is such that initial segment of *h--*h1 is invisible
        //    2 = as 1, except that initial segment is visible
        //    3 = *h lies outside mask
        //   *pp is as follows:
        //   - if a Maskpoint coincides with *h, *pp is it
        //   - if the last two Maskpoint "before" (i.e., just counter-clockwise from) *h have the same a value, *pp is the first of these
        //   - otherwise, *pp is the last Maskpoint "before" *h
        //*/
        //{
        //    Maskpoint* p,*q;
        //    double a;
        //    double dx, dy;
        //    double x, y;
        //    x = h->x; y = h->y; a = h->a;

        //    if (h->maskpoint)
        //    {
        //        p = h->maskpoint;
        //        q = NEXT(p);
        //        goto found;
        //    }

        //    p = Q;
        //    if (a >= 0) p = Q2;
        //    if (p->a > a)
        //    {
        //        while (p->next) p = p->next;
        //    }
        //    else
        //    {
        //        while (p->next && p->next->a < a) p = p->next;
        //        if (p->next && p->next->a == a && !p->a == a) p = p->next;
        //    }
        //    q = NEXT(p);
        //    /*
        //       if (p->a==a && fabs(p->s-h->s)<1e-4) { 
        //    */
        //    if (p->a == a && p->s == h->s)
        //    {
        //    found:
        //        *pp = p;
        //        if (!h1) return 1;
        //        if (p->a == q->a)
        //        {
        //            if (q->s < p->s) return 2; /* this assumes that *h1 is clockwise from *h */
        //            else return 1;
        //        }
        //        dx = q->x - p->x;
        //        dy = q->y - p->y;
        //        a = (-dy * h1->x + dx * h1->y) / (-dy * p->x + dx * p->y);
        //        if (a <= 1) return 1;
        //        return 2;
        //    }
        //    /*
        //       if (q->a==a && fabs(q->s-h->s)<1e-4) { 
        //    */
        //    if (q->a == a && q->s == h->s)
        //    {
        //        p = q;
        //        q = NEXT(p);
        //        goto found;
        //    }

        //    *pp = p;
        //    if (a == p->a && a == q->a)
        //    {
        //        if (h->s >= p->s && h->s >= q->s) return 3;
        //        if (h->s >= q->s && h1) return 3;  /* this assumes that *h1 is clockwise from *h */
        //        return 0;
        //    }
        //    dx = q->x - p->x;
        //    dy = q->y - p->y;
        //    a = (-dy * x + dx * y) / (-dy * p->x + dx * p->y);
        //    if (a <= 1) return 0;
        //    return 3;
        //}



        //static void recalc_min_s(void)
        //{
        //    Maskpoint* p,*q;
        //    double mm;
        //    Maskpoint* mmp;

        //    if (right_half)
        //    {
        //        p = Q;
        //        q = p;
        //        min_s = p->s; mmp = p;
        //        if (left_half)
        //        {
        //            while ((p = p->next))
        //            {
        //                mm = (mmp = p)->s; if (mm > q->s) mm = (mmp = q)->s; mm *= cos((p->a - q->a) / 2); if (mm < min_s) { min_s = mm; min_s_mp = mmp; }
        //                q = p;
        //            }
        //            p = Q;
        //            mm = (mmp = p)->s; if (mm > q->s) mm = (mmp = q)->s; mm *= cos(drem(p->a - q->a, 2 * M_PI) / 2); if (mm < min_s) { min_s = mm; min_s_mp = mmp; }
        //        }
        //        else
        //        {
        //            while ((p = p->next))
        //            {
        //                mm = (mmp = p)->s; if (mm > q->s) mm = (mmp = q)->s; mm *= cos((p->a - q->a) / 2); if (mm < min_s) { min_s = mm; min_s_mp = mmp; }
        //                if (p->a >= 1e-4) break;
        //                q = p;
        //            }
        //        }
        //    }
        //    else
        //    {
        //        p = Q2->next; min_s = p->s; mmp = p;
        //        q = p;
        //        while ((p = p->next))
        //        {
        //            mm = (mmp = p)->s; if (mm > q->s) mm = (mmp = q)->s; mm *= cos((p->a - q->a) / 2); if (mm < min_s) { min_s = mm; min_s_mp = mmp; }
        //            q = p;
        //        }
        //        if (Q->a < -M_PI + 1e-4)
        //        {
        //            p = Q;
        //            mm = (mmp = p)->s; if (mm > q->s) mm = (mmp = q)->s; mm *= cos(drem(p->a - q->a, 2 * M_PI) / 2); if (mm < min_s) { min_s = mm; min_s_mp = mmp; }
        //        }
        //    }
        //}


        //static void commit_mask(void)
        //{
        //    Maskpoint** pp,*p,*prev;
        //    int need_recalc_min_s = 0;
        //    PRINTF(("commit!\n"));
        //    while (new)
        //    {
        //        Maskpoint* q;
        //        q = new;
        //        new= new->next;
        //        add_to_mask(q);
        //    }
        //    pp = &Q;
        //    prev = NULL;
        //    while (*pp)
        //    {
        //        p = *pp;
        //        if (p->del)
        //        {
        //            need_recalc_min_s |= (p == min_s_mp);
        //            if (p == Q2) { Q2 = prev; if (!Q2) Q2 = Q; }
        //            if (p->hp) p->hp->maskpoint = NULL;
        //            p = *pp;
        //            *pp = p->next;
        //            mpfree(p);
        //        }
        //        else
        //        {
        //            pp = &(p->next);
        //            prev = p;
        //        }
        //    }
        //    if (need_recalc_min_s) recalc_min_s();
        //}


        //#ifdef DEBUG
        //static void draw_mask(void)
        ///* for debugging purposes only: draw the mask (slightly enlarged, so it doesn't completely overwrite the actual lines) */
        //{
        //    Maskpoint* p = Q;
        //    double x0, y0;

        //    SetForeground(&c_hfield);
        //    x0 = Q->x;
        //    y0 = Q->y;
        //    PRINTF(("%8g %8g %8g\n", Q->x, Q->y, Q->a));
        //    p = p->next;
        //    while (p)
        //    {
        //        DrawLine(1.01 * x0 + Xtr2, 1.01 * y0 + Ytr2, 1.01 * p->x + Xtr2, 1.01 * p->y + Ytr2);
        //        PRINTF(("%8g %8g %8g\n", p->x, p->y, p->a));
        //        x0 = p->x;
        //        y0 = p->y;
        //        p = p->next;
        //    }
        //    DrawLine(1.01 * x0 + Xtr2, 1.01 * y0 + Ytr2, 1.01 * Q->x + Xtr2, 1.01 * Q->y + Ytr2);

        //    PRINTF(("left=%i right=%i min_s=%f  Q2->a=%f\n", left_half, right_half, min_s, Q2->a));
        //}
        //#endif


        ///* ------------------------- code for actually drawing the non-masked parts of a line ----------------------------------------------------- */


        //static int intersect(Maskpoint* p, Maskpoint* q, double x0, double y0, double x1, double y1, double* xx, double* yy)
        ///* calculate intersection of line between *p and *q and line between x0y0 and x1y1; coordinates of intersection go into *xx,*yy;
        //   return values:
        //    0: no intersection
        //    1: intersection, x0y0 visible
        //    2: intersection, x0y0 not visible
        //    3: endpoint identical
        //*/
        //{
        //    double t, tnum, tden;
        //    double nx, ny;
        //    double px, py;
        //    int type;

        //    PRINTF(("intersect: (%5f,%5f)-(%5f,%5f) , (%5f,%5f)-(%5f,%5f)", p->x, p->y, q->x, q->y, x0, y0, x1, y1));
        //    if ((p->x == x0 && p->y == y0) || (q->x == x0 && q->y == y0)) { *xx = x0; *yy = y0; type = 3; goto done; }
        //    if ((p->x == x1 && p->y == y1) || (q->x == x1 && q->y == y1)) { *xx = x1; *yy = y1; type = 3; goto done; }

        //    /* part of the Cyrus-Beck algorithm for calculating intersections; see Foley/van Dam/Feiner/Hughes, 1996 */
        //    nx = q->y - p->y;
        //    ny = p->x - q->x;
        //    tnum = nx * (x0 - p->x) + ny * (y0 - p->y);
        //    tden = -nx * (x1 - x0) - ny * (y1 - y0);
        //    if (tden == 0)
        //    {
        //        /* parallel lines */
        //        type = 0;
        //    }
        //    else
        //    {
        //        t = tnum / tden;
        //        if (t <= 0 || t >= 1)
        //        {
        //            type = 0;
        //        }
        //        else
        //        {
        //            /* here we depart from Cyrus-Beck, since we're only trying to intersect with one line segment */
        //            double a;
        //            *xx = x0 + t * (x1 - x0);
        //            *yy = y0 + t * (y1 - y0);
        //            px = q->x - p->x;
        //            py = q->y - p->y;
        //            a = px * *xx + py * *yy;
        //            if (a <= px * p->x + py * p->y || a >= px * q->x + py * q->y)
        //            {
        //                type = 0;
        //            }
        //            else
        //            {
        //                if (tden > 0) type = 1;
        //                else type = 2;
        //            }
        //        }
        //    }
        //done:
        //    PRINTF((" -> %i (%5f,%5f)\n", type, *xx, *yy));
        //    return type;
        //}



        //int cnt0 = 0, cnt1 = 0;


        //static void draw_masked1(Hidpoint* c0, Hidpoint* c1)
        ///* actually draw a (masked) line between two points; the mask is not yet updated, but preparations
        //   are made so the next commit_mask() call will update the mask.
        //*/
        //{
        //    Maskpoint* p0,*p1,*p;
        //    double x, y;
        //    int i0, i1;
        //    int visible;
        //    int skip = 0;
        //    double a0, a1;
        //    double x0, x1, y0, y1;
        //    double s0, s1;
        //    double s;

        //    if (c0->s < min_s && c1->s < min_s) return;

        //    a0 = c0->a;
        //    a1 = c1->a;

        //    if ((a0 > a1) ^ (fabs(a0 - a1) > M_PI))
        //    {
        //        /* if needed, exchange endpoints to make sure that we're moving in a clockwise direction */
        //        double h;
        //        Hidpoint* hh;
        //        h = a0; a0 = a1; a1 = h;
        //        hh = c0; c0 = c1; c1 = hh;
        //    }
        //    x0 = c0->x;
        //    y0 = c0->y;
        //    x1 = c1->x;
        //    y1 = c1->y;

        //    i0 = find_in_mask(c0, &p0, c1);
        //    i1 = find_in_mask(c1, &p1, NULL);
        //    PRINTF(("draw_masked1(%f,%f,%f,%f) %i %i\n", x0, y0, x1, y1, i0, i1));

        //    if (a0 == a1)
        //    {
        //        /* exactly radial line; this needs to be handled as a special case, because letting the normal code handle it causes problems
        //           because the resulting Maskpoints can no longer be sorted sensibly based on their (identical!) a values.
        //           (Note: exactly radial lines are pretty unlikely, except on the z axis.) */
        //        if (i0 != 3 && i1 != 3) return;
        //        if (i0 + i1 < 6)
        //        {
        //            double xx, yy;
        //            int type;
        //            Maskpoint* mp;
        //            if (p0->a != a0)
        //            {
        //                type = intersect(p0, NEXT(p0), x0, y0, x1, y1, &xx, &yy);
        //                if (type == 1 || type == 2)
        //                {
        //                    mp = newmaskpoint(xx, yy, ATAN(yy, xx), xx * xx + yy * yy, NULL);
        //                    mp->after = -1;         /* note: we only add the intersection to the mask; the line's visible end point is also new, but here we wouldn't know the order of the two; that's also the reason we set after to -1 */
        //                    add_to_newlist(mp);
        //                    if (i0 != 3) { x0 = xx; y0 = yy; }
        //                    else { x1 = xx; y1 = yy; }
        //                }
        //            }
        //            else
        //            {
        //                if (NEXT(p0)->a == a0)
        //                {
        //                    if (p0->s > NEXT(p0)->s) { xx = p0->x; yy = p0->y; }
        //                    else { xx = NEXT(p0)->x; yy = NEXT(p0)->y; }
        //                }
        //                else
        //                {
        //                    xx = p0->x; yy = p0->y;
        //                }
        //                if (i0 != 3) { x0 = xx; y0 = yy; }
        //                else { x1 = xx; y1 = yy; }
        //            }
        //        }
        //        DrawLine(x0 + Xtr2, y0 + Ytr2, x1 + Xtr2, y1 + Ytr2);
        //        return;
        //    }

        //    if (i0 == 3)
        //    {
        //        /* if the initial point is outside the mask, it must become a new corner point of the mask */
        //        Maskpoint* mp;
        //        mp = newmaskpoint(x0, y0, a0, c0->s, c0);
        //        mp->after = 1;
        //        add_to_newlist(mp);
        //    }

        //    if (p0->a == a0 && NEXT(p0)->a == a0)
        //    {
        //        /* if the initial point is at the same angle as a radial piece of the mask, there are several cases:
        //           - initial point == first mask corner point (out of the two that form that radial piece of mask), and initial segment visible: remove second mask point
        //           - initial point == first mask corner point, and initial segment not visible: no changes to mask
        //           - initial point is inside the mask: no changes
        //           - initial point is outside the mask: mask is extended, so drop second mask point (because our line segment is clock wise)
        //           - initial point == second mask corner point: it that case, p0==second mask corner point, and this code is not reached because the if-condition is not satisfied
        //        */
        //        if (i0 >= 2) NEXT(p0)->del = 1;
        //        /* In all cases except the last (but in that case we don't get here), we should make sure that p0 points
        //           to the second mask point, so the while()-loop starts with a non-radial segment.
        //           However, the next block also conditionally skips one point, so we need an extra condition here to prevent duplication.
        //        */
        //        if (i0 != 1 && i0 != 2) p0 = NEXT(p0);
        //    }

        //    if (i0 == 1 || i0 == 2)
        //    {
        //        /* if the initial point is a mask point, the following cases are possible: */
        //        if (p0 == p1)
        //        {
        //            /* final point is in next mask segment:
        //               only if final point is outside the mask, does anything need to be drawn; and we don't need the while loop for that
        //            */
        //            if (i1 != 3) return;
        //            skip = 1;
        //        }
        //        else
        //        {
        //            /* final point is in later segment: this is just a normal case, except that we don't need to intersect with first mask
        //               segment, since that intersection would be the initial point
        //            */
        //            p0 = NEXT(p0);
        //        }
        //    }



        //    s0 = c0->s;
        //    s1 = c1->s;

        //    x = x0; y = y0; s = s0;
        //    if (s1 > s) s = s1;
        //    visible = (i0 >= 2);

        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);

        //    p = p0;
        //    if (!skip) while (1)
        //        {
        //            /* this is a while loop over all mask segments that potentially intersect the line we're drawing;
        //               we keep track of whether our line is visible in the variable 'visible', calculating intersections as we go
        //            */
        //            Maskpoint* q;
        //            double xx, yy;
        //            int type;

        //            q = NEXT(p);

        //            cnt0++;
        //            type = intersect(p, q, x, y, x1, y1, &xx, &yy);
        //            switch (type)
        //            {
        //                case 0:     /* no intersection */
        //                case 3:     /* endpoint identical */
        //                    if (visible)
        //                    {
        //                        if (p != p0 || i0 == 1 || i0 == 2) p->del = 1;
        //                        if ((i1 == 0 || i1 == 3) && p != p1) q->del = 1;
        //                        if ((i1 == 1 || i1 == 2) && q != p1) q->del = 1;
        //                        PRINTF(("visible, no intersection: %c%f %f %f %c%f\n", p->del ? '*' : ' ', p->a, p0->a, p1->a, q->del ? '*' : ' ', q->a));
        //                    }
        //                    break;
        //                case 1:     /* intersection, x0y0 visible */
        //                    DrawLine(x + Xtr2, y + Ytr2, xx + Xtr2, yy + Ytr2);
        //                    x = xx; y = yy; s = x * x + y * y; if (s1 > s) s = s1;
        //                    visible = 0;
        //                    if (p == p0 && i0 == 2) p->del = 1;
        //                    add_to_newlist(newmaskpoint(xx, yy, ATAN(yy, xx), xx * xx + yy * yy, NULL));
        //                    break;
        //                case 2:     /* intersection, x0y0 not visible */
        //                    x = xx; y = yy; s = x * x + y * y; if (s1 > s) s = s1;
        //                    visible = 1;
        //                    add_to_newlist(newmaskpoint(xx, yy, ATAN(yy, xx), xx * xx + yy * yy, NULL));
        //                    break;
        //            }

        //            if (p == p1) break;
        //            p = q;
        //            if ((i1 == 1 || i1 == 2 || a1 == p1->a) && p == p1) break;
        //        }
        //    if (visible)
        //    {
        //        /* draw the last segment of the line, if it is visible */
        //        DrawLine(x + Xtr2, y + Ytr2, x1 + Xtr2, y1 + Ytr2);
        //        add_to_newlist(newmaskpoint(x1, y1, a1, c1->s, c1));
        //        PRINTF(("visible, testing for del: %i %g %g %g\n", i1, a1, p1->a, NEXT(p1)->a));
        //        if ((i1 == 0 || i1 == 3) && p1->a == a1 && NEXT(p1)->a == a1) p1->del = 1;
        //    }
        //}



        ///* ------------------------- finally: the code for drawing the gain pattern ----------------------------------------------------- */


        //static void first_quadrangle(Radpattern* rp, Hidpoint** hp, int* i0, int* i1, int* j0, int* j1)
        ///* Draw the first quadrangle: this is the polygon which is intersected by the
        //   line from the origin to the viewer's eye, i.e., the line orthogonal to
        //   the screen or paper.
        //   Furthermore, this function fills the Hidpoint array; this is not really a
        //   logical place to do this, but:
        //   - we can't do it earlier since this function may modify phi2 a little to
        //     prevent the viewing direction from passing through an edge;
        //   - we already in this function need four elements of the Hidpoint array to
        //     actually draw the quadrangle; then it makes sense to calculate them all.
        //*/
        //{
        //    int i, j, k, l;
        //    Hidpoint h0, h1;
        //    Hidpoint* c0,*c1;

        //    int d;
        //    double y, yprev;

        //    /* find the rp->gphi value closest to the phi2 of the viewing direction: */
        //    j = 0;
        //    for (l = 1; l < rp->numphi; l++) if (fabs(drem(rp->gphi[l] - phi2, 2 * M_PI)) < fabs(drem(rp->gphi[j] - phi2, 2 * M_PI))) j = l;
        //    /* if almost identical to phi2, change phi2 a little */
        //    if (fabs(drem(rp->gphi[j] - phi2, 2 * M_PI)) < 1e-6)
        //    {
        //        if (drem(rp->gphi[j] - phi2, 2 * M_PI) < 0) phi2 -= 1e-6;
        //        else phi2 += 1e-6;
        //        phi2 = drem(phi2, 2 * M_PI);
        //        calcproj2();
        //    }
        //    /* choose j and l such that they are consecutive phi indices around the viewing direction */
        //    if (drem(rp->gphi[j] - phi2, 2 * M_PI) < 0)
        //    {
        //        l = (j + 1) % rp->numphi;
        //    }
        //    else
        //    {
        //        l = j;
        //        j = (j - 1 + rp->numphi) % rp->numphi;
        //    }

        //    /* find the rp->gtheta value closest to the theta2 of the viewing direction: */
        //    i = 0;
        //    for (k = 1; k < rp->numtheta; k++) if (fabs(rp->gtheta[k] - theta2) < fabs(rp->gtheta[i] - theta2)) i = k;
        //    /* now calculate the y (screen coordinates) of the intersection of the corresponding
        //       "horizontal" segment with the (on screen) vertical through the origin: */
        //    proj2(&rp->gpo[j][i], &h0);
        //    proj2(&rp->gpo[l][i], &h1);
        //    y = ((h1.y - h0.y) / (h1.x - h0.x)) * (-h0.x) + h0.y;
        //    /* We need to choose subsequent theta indices such that the quadrangle
        //       surrounds the origin (in the projection on the screen).
        //       We test this by calculating the intersections of the "horizontal"
        //       segments with the (on-screen) vertical through the origin. */
        //    if (y > 0) d = -1;
        //    else d = +1;
        //    do
        //    {
        //        k = i;
        //        yprev = y;
        //        i += d;
        //        proj2(&rp->gpo[j][i], &h0);
        //        proj2(&rp->gpo[l][i], &h1);
        //        y = ((h1.y - h0.y) / (h1.x - h0.x)) * (-h0.x) + h0.y;
        //    } while (y * d < 0);
        //    if (i > k)
        //    {
        //        int h;
        //        double hh;
        //        h = i; i = k; k = h;
        //        hh = y; y = yprev; yprev = hh;
        //    }

        //    {
        //        /* fill the Hidpoint array */
        //        int i, j;
        //        for (i = 0; i < rp->numtheta; i++)
        //            for (j = 0; j < rp->numphi; j++)
        //            {
        //                proj2(&rp->gpo[j][i], &hp[j][i]);
        //                hp[j][i].c = rp->gain[j][i].axial;
        //                hp[j][i].maskpoint = NULL;
        //            }
        //    }

        //    /* draw the quadrangle, and set up the mask */
        //# ifdef DEBUG
        //    SetForeground(&c_wire);
        //#endif
        //    c0 = &hp[l][i]; Q2 = Q = newmaskpoint(c0->x, c0->y, c0->a, c0->s, c0);
        //    c1 = &hp[j][i]; add_to_mask(newmaskpoint(c1->x, c1->y, c1->a, c1->s, c1));
        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //    DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //    c0 = &hp[j][k]; add_to_mask(newmaskpoint(c0->x, c0->y, c0->a, c0->s, c0));
        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //    DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //    c1 = &hp[l][k]; add_to_mask(newmaskpoint(c1->x, c1->y, c1->a, c1->s, c1));
        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //    DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //    c0 = &hp[l][i];
        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //    DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //# ifdef DEBUG
        //    SetForeground(&c_gain);
        //#endif
        //    while (Q2->next && Q2->next->a < 0) Q2 = Q2->next;

        //    /* if a "horizontal" edge of the quadrangle passes through the origin (on screen),
        //       add a maskpoint to make sure the mask really surrounds the origin: */
        //    if (fabs(yprev) < 1e-5) add_to_mask(newmaskpoint(0, -0.01, ATAN(-0.01, 0), 0.0001, NULL));
        //    if (fabs(y) < 1e-5) add_to_mask(newmaskpoint(0, 0.01, 0, 0.0001, NULL));

        //    /* caller wants to know where in the *rp array the initial quadrangle is: */
        //    *i0 = i;
        //    *i1 = k;
        //    *j0 = j;
        //    *j1 = l;
        //}


        //static void first_polygon(Radpattern* rp, Hidpoint** hp, int* i0, int* i1, int* j0, int* j1)
        ///* Draw the first polygon, in cases where data is only available for half of the
        //   sphere (namely theta <= 90 degrees), and we're looking from a direction
        //   outside that range (i.e., from "down under").
        //   In this case, the first polygon is simply the set of points at theta=90 degrees.
        //*/
        //{
        //    int i, j, l;
        //    Hidpoint* c0,*c1;

        //    i = rp->numtheta - 1;
        //    *i0 = i;
        //    *i1 = i;

        //    /* find the rp->gphi value closest to the phi2 of the viewing direction: */
        //    j = 0;
        //    for (l = 1; l < rp->numphi; l++) if (fabs(drem(rp->gphi[l] - phi2, 2 * M_PI)) < fabs(drem(rp->gphi[j] - phi2, 2 * M_PI))) j = l;
        //    /* if almost identical to phi2, change phi2 a little */
        //    if (fabs(drem(rp->gphi[j] - phi2, 2 * M_PI)) < 1e-6)
        //    {
        //        if (drem(rp->gphi[j] - phi2, 2 * M_PI) < 0) phi2 -= 1e-6;
        //        else phi2 += 1e-6;
        //        phi2 = drem(phi2, 2 * M_PI);
        //        calcproj2();
        //    }
        //    /* choose *j0 and *j1 such that they are subsequent phi indices around the viewing direction */
        //    if (drem(rp->gphi[j] - phi2, 2 * M_PI) < 0)
        //    {
        //        *j0 = j;
        //        *j1 = (j + 1) % rp->numphi;
        //    }
        //    else
        //    {
        //        *j1 = j;
        //        *j0 = (j - 1 + rp->numphi) % rp->numphi;
        //    }

        //    {
        //        /* fill the Hidpoint array */
        //        int i, j;
        //        for (i = 0; i < rp->numtheta; i++)
        //            for (j = 0; j < rp->numphi; j++)
        //            {
        //                proj2(&rp->gpo[j][i], &hp[j][i]);
        //                hp[j][i].c = rp->gain[j][i].axial;
        //                hp[j][i].maskpoint = NULL;
        //            }
        //    }


        //    c1 = &hp[0][i]; Q = Q2 = newmaskpoint(c1->x, c1->y, c1->a, c1->s, c1);
        //    for (j = 1; j < rp->numphi; j++)
        //    {
        //        Q2 = Q;  /* during initial setup of mask, using Q2 is not reliable otherwise (and we don't really need the speedup here) */
        //        c0 = c1;
        //        c1 = &hp[j][i]; add_to_mask(newmaskpoint(c1->x, c1->y, c1->a, c1->s, c1));
        //        if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //        DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //    }
        //    Q2 = Q;
        //    c0 = c1;
        //    c1 = &hp[0][i]; add_to_mask(newmaskpoint(c1->x, c1->y, c1->a, c1->s, c1));
        //    if (polarization == POLcolour) setpolcolour((c0->c + c1->c) / 2);
        //    DrawLine(c0->x + Xtr2, c0->y + Ytr2, c1->x + Xtr2, c1->y + Ytr2);
        //    while (Q2->next && Q2->next->a < 0) Q2 = Q2->next;
        //}



        //void draw_opaque(Radpattern* rp)
        //{
        //    int i, j;
        //    Hidpoint* c0,*c1,*c2,*c3;
        //    int i0, i1;
        //    int j0, j1;
        //    Hidpoint** hp;
        //    int halfsphere = 0;
        //    int lookfrombottom = 0;
        //    Hidpoint hpnul1 = { 0, 0, -M_PI, 0, 1, NULL };
        //    Hidpoint hpnul2 = { 0, 0, 0, 0, 1, NULL };


        //#ifdef MEASURE
        //   struct rusage ru1,ru2;
        //   getrusage(RUSAGE_SELF,&ru1);
        //   for (count=99;count>=0;count--) {
        //#endif

        //   /* Prepare the viewing direction: in principle just phi,theta, but the
        //      algorithm doesn't work well when looking exactly from the "top" or
        //      "bottom", i.e., when theta is 0 or 180 degrees.
        //      Note that phi2 may be changed a little bit later on in first_polygon()
        //      or first_quadrangle(). */
        //   phi2=phi* (M_PI/180.);
        //   theta2=theta* (M_PI/180.);
        //   if (theta2<1e-7) theta2=1e-7;
        //   else if (theta2>M_PI-1e-7) theta2=M_PI-1e-7;
        //   if (fabs(theta2-M_PI/2)<1e-7) theta2=M_PI/2-1e-7;

        //   calcproj2();

        //   if (rp->gtheta[rp->numtheta - 1]<M_PI-1e-7) halfsphere=1;

        //   hp=mymalloc(rp->numphi*sizeof(Hidpoint*));
        //   for (j=0;j<rp->numphi;j++)
        //      hp[j]=mymalloc(rp->numtheta*sizeof(Hidpoint));

        ///* initialize the mask, drawing the first polygon */
        //init_mask();
        //   if (!halfsphere || theta<=90) {
        //      first_quadrangle(rp, hp,&i0,&i1,&j0,&j1);
        //   } else {
        //      first_polygon(rp, hp,&i0,&i1,&j0,&j1);
        //lookfrombottom=1;
        //   }
        //   recalc_min_s();

        //PRINTF(("%i %i\n", i, j));

        //   /* next, complete the quadrangles "above" and "below" the initial quadrangle, up to the "poles" */
        //   c0=&hp[j0][i1];
        //   c1=&hp[j1][i1];
        //   for (i=i1+1;i<rp->numtheta;i++) {
        //      c2=&hp[j0][i];
        //      c3=&hp[j1][i];
        //      draw_masked1(c0, c2);
        //draw_masked1(c1, c3);
        //draw_masked1(c2, c3);
        //commit_mask();
        //c0=c2; c1=c3;
        //   }
        //   c0=&hp[j0][i0];
        //   c1=&hp[j1][i0];
        //   for (i=i0-1;i>=0;i--) {
        //      c2=&hp[j0][i];
        //      c3=&hp[j1][i];
        //      draw_masked1(c0, c2);
        //draw_masked1(c1, c3);
        //draw_masked1(c2, c3);
        //commit_mask();
        //c0=c2; c1=c3;
        //   }
        //   j=j1;

        //   /* Next, we do all quadrangles "right" of the vertical line through the origin,
        //      "column" by column, starting with the column just to the right of the initial
        //      quadrangle.
        //      Within a column, we need to be a bit careful: depending on the precise location
        //      of the points, we either need to go downward (starting at the "north pole") or
        //      upward (from the "south pole"). In fact, we usually need to change direction
        //      one or more times. This is needed to ensure that the mask is a single-valued
        //      function of the angle a .
        //   */
        //   left_half=0;
        //   recalc_min_s();
        //   while (1) {
        //      int l;
        //int ii;
        //int jn = (j + 1) % rp->numphi;
        //      if (drem(rp->gphi[jn]-phi2,2* M_PI)<0) break;
        //      ii=0;
        //      for (i=1;i<rp->numtheta-1;i++) {
        //         double a1, a2;
        //DEBUGPAUSE
        //         c1 = &hp[j][i];
        //a1 = c1->a;
        //         c2=&hp[jn][i];
        //         a2 = c2->a;
        //         if (drem(a2-a1,2* M_PI)<0 ) {
        //            draw_masked1(c1, c2);
        //c0=&hp[jn][i - 1];
        //            draw_masked1(c0, c2);
        //l=i-1;
        //            while (l>ii) {
        //               c2=c0;
        //               c1=&hp[j][l];
        //               draw_masked1(c1, c2);
        //commit_mask();
        //c0=&hp[jn][l - 1];
        //               draw_masked1(c0, c2);
        //l--;
        //            }
        //            if (l==0) draw_masked1(c0,&hpnul1);
        //commit_mask();
        //ii=i;
        //         }
        //      }
        //      c2=&hp[jn][i];
        //      if (halfsphere && !lookfrombottom) {
        //         c1=&hp[j][i];
        //         draw_masked1(c1, c2);
        //      }
        //      if (!halfsphere) {
        //         draw_masked1(c2,&hpnul2);
        //      }
        //      c0=&hp[jn][i - 1];
        //      draw_masked1(c0, c2);
        //l=i-1;
        //      while (l>ii) {
        //         c2=c0;
        //         c1=&hp[j][l];
        //         draw_masked1(c1, c2);
        //commit_mask();
        //c0=&hp[jn][l - 1];
        //         draw_masked1(c0, c2);
        //l--;
        //      }
        //      commit_mask();
        //j=jn;
        //   }
        //   j1=j;

        //   /* now exactly the same procedure, but for columns to the left of the origin */
        //   j=j0;
        //   left_half=1;
        //   right_half=0;
        //   recalc_min_s();

        //   while (1) {
        //      int l;
        //int ii;
        //int jn = (j - 1 + rp->numphi) % rp->numphi;
        //      if (drem(rp->gphi[jn]-phi2,2* M_PI)>0) break;
        //      ii=0;
        //      for (i=1;i<rp->numtheta-1;i++) {
        //         double a1, a2;
        //DEBUGPAUSE
        //         c1 = &hp[j][i];
        //a1 = c1->a;
        //         c2=&hp[jn][i];
        //         a2 = c2->a;
        //         if (drem(a2-a1,2* M_PI)>0 ) {
        //            draw_masked1(c1, c2);
        //c0=&hp[jn][i - 1];
        //            draw_masked1(c0, c2);
        //l=i-1;
        //            while (l>ii) {
        //               c2=c0;
        //               c1=&hp[j][l];
        //               draw_masked1(c1, c2);
        //commit_mask();
        //c0=&hp[jn][l - 1];
        //               draw_masked1(c0, c2);
        //l--;
        //            }
        //            commit_mask();
        //ii=i;
        //         }
        //      }
        //      c2=&hp[jn][i];
        //      if (halfsphere && !lookfrombottom) {
        //         c1=&hp[j][i];
        //         draw_masked1(c1, c2);
        //      }
        //      c0=&hp[jn][i - 1];
        //      draw_masked1(c0, c2);
        //l=i-1;
        //      while (l>ii) {
        //         c2=c0;
        //         c1=&hp[j][l];
        //         draw_masked1(c1, c2);
        //commit_mask();
        //c0=&hp[jn][l - 1];
        //         draw_masked1(c0, c2);
        //l--;
        //      }
        //      commit_mask();
        //j=jn;
        //   }
        //   j0=j;
        //   right_half=1;
        //   recalc_min_s();

        ///* finally, draw the set of quadrangles that are "behind" the (3D) Z axis */
        //DEBUGPAUSE
        //   if (halfsphere) {
        //      if (!lookfrombottom) { 
        //         i0=i1=rp->numtheta;
        //      } else {
        //         i1=1;
        //         while (rp->gtheta[i1]<M_PI-theta2) i1++;
        //         i0=i1-1;
        //      }
        //   } else {
        //      i0=rp->numtheta-1-i0;
        //      i1=rp->numtheta-i1;
        //   }
        //   for (i=1;i<i1;i++) {
        //      c2=&hp[j0][i];
        //      c3=&hp[j1][i];
        //      draw_masked1(c2, c3);
        //commit_mask();
        //   }
        //   for (i=rp->numtheta-2;i>=i0;i--) {
        //      c2=&hp[j0][i];
        //      c3=&hp[j1][i];
        //      draw_masked1(c2, c3);
        //commit_mask();
        //   }
        //DEBUGPAUSE

        //PRINTF(("\n"));

        //   cleanup_mask();
        //   for (j=0;j<rp->numphi;j++) free(hp[j]);
        //free(hp);

        //#ifdef MEASURE
        //   }
        //   getrusage(RUSAGE_SELF,&ru2);
        //printf("duration = %i\n", ru2.ru_utime.tv_usec-ru1.ru_utime.tv_usec + 1000000*(ru2.ru_utime.tv_sec-ru1.ru_utime.tv_sec));
        //   printf("cnt1/cnt0 = %g\n", (double) cnt1/cnt0);
        //#endif
        //}



        //char* opaque_impossible(Radpattern* rp)
        ///* returns NULL if opaque drawing of this Radpattern is possible;
        //   otherwise, returns a pointer to a string explaining why opaque drawing is not possible.
        //*/
        //{
        //    double a;
        //    if (rp->gtheta[0] * 180 / M_PI > 0.1) return "theta range does not start at 0 degrees";
        //    a = rp->gtheta[rp->numtheta - 1] * 180 / M_PI;
        //    if (a > 180.1) return "theta range continues beyond 180 degrees";
        //    if ((a < 89.9 || a > 90.1) && a < 179.9) return "theta range does not go on until either 90 or 180 degrees";
        //    if (rp->gphi[0] * 180 / M_PI > 0.1) return "phi range does not start at 0 degrees";
        //    if (rp->gphi[rp->numphi - 1] * 180 / M_PI < 271) return "phi range does not go on until almost 360 degrees";
        //    if (rp->gphi[rp->numphi - 1] * 180 / M_PI > 359.9) return "phi range goes on until 360 degrees or beyond";
        //    if (rp->numtheta < 4) return "not enough theta values";
        //    if (rp->numphi < 4) return "not enough phi values";
        //    return NULL;
        //}


    }
}
