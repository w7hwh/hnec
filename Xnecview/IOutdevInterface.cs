using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using System.Threading.Tasks;

namespace Necview
{
    interface IOutdevInterface
    {
        void SetLineAttributes(uint a, int b, int c, int d);
        void DrawLine(double a, double b, double c, double d);
        void SetForeground(Color color);
        void ClearWindow();
        void DrawString(double a, double b, string c, double d, double e);
        void Complete();
        void SetClipRectangle(double a, double b, double c, double d);
        void ClearRectangle(double a, double b, double c, double d);
    }

    //extern Outdev *out;   /* pointer to Outdev struct, determines destination of drawing commands issued through the below macros: */

    //#define SetLineAttributes(a,b,c,d) out->SetLineAttributes((a),(b),(c),(d))
    //#define DrawLine(a,b,c,d) out->DrawLine((a),(b),(c),(d))
    //#define SetForeground(a) out->SetForeground(a)
    //#define ClearWindow() out->ClearWindow()
    //#define DrawStringLL(a,b,c) DrawString(a,b,c,0,0)     /* draw text specifying lower left corner */
    //#define DrawStringUL(a,b,c) DrawString(a,b,c,0,1)     /* draw text specifying upper left corner */
    //#define DrawString(a,b,c,d,e) out->DrawString(a,b,c,d,e) /* draw text specifying a corner defined by d and e:
    //                                                            d = 0...1 for left...right; e = 0...1 for lower...upper*/
    //#define SetClipRectangle(a,b,c,d) out->SetClipRectangle(a,b,c,d)  /* set clip rectangle to (a,b)..(c,d) */
    //#define ClearRectangle(a,b,c,d) out->ClearRectangle(a,b,c,d)  /* clear part of window */

    //#define c_bg    col[0]
    //#define c_axis  col[1]
    //#define c_wire  col[2]
    //#define c_surf  col[3]
    //#define c_back  col[4]
    //#define c_gain  col[5]
    //#define c_scale col[6]
    //#define c_exci  col[7]
    //#define c_load  col[8]
    //#define c_netw  col[9]
    //#define c_inactive col[10]
    //#define c_efield   col[11]
    //#define c_hfield   col[12]
    //#define c_poynting col[13]
    //#define c_qpos     col[14]
    //#define c_qneg     col[15]
    //#define c_gain_lin   col[16]
    //#define c_gain_lhcp  col[17]
    //#define c_gain_rhcp  col[18]
    //#define c_currents (col+19)
    //#define NC NC_phase+19
    //extern GdkColor col[NC];
}
