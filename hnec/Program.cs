using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using Necview;
using System.Windows.Forms.DataVisualization.Charting;
using Microsoft.Win32;      // for screen scaling value
using System.Text.RegularExpressions;


namespace hnec
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
#if (DEBUG)
            //System.Console.WriteLine("OSVersion = {0}", System.Environment.OSVersion);
            //System.Console.WriteLine("ProcessorCount = {0}", System.Environment.ProcessorCount);
            //System.Console.WriteLine("Version = {0}", System.Environment.Version);
            //System.Console.WriteLine("UserName = {0}", System.Environment.UserName);
            //System.Console.WriteLine("MachineName = {0}", System.Environment.MachineName);
            //System.Console.WriteLine("CurrentDirectory = {0}", System.Environment.CurrentDirectory);

            ReadNecInput ni = new ReadNecInput();
            ni.Read_nec_input(@"c:\rcp\test.nec");

            //DateTime now = DateTime.Now;
            //System.Console.WriteLine("Now = {0}{1}{2}  {3}:{4}:{5}", now.Year, now.Month, now.Day, now.Hour, now.Minute, now.Second);
            //System.Console.WriteLine("now.ToString() = {0}", now.ToString());

            //var currentDPI = (int)Registry.GetValue("HKEY_CURRENT_USER\\Control Panel\\Desktop", "LogPixels", 96);
            //System.Console.WriteLine("LogPixels = {0}", currentDPI);

            //string url = "http://www.contoso.com:8080/letters/readme.html";
            //Regex r = new Regex(@"^(?<proto>\w+)://[^/]+?(?<port>:\d+)?/", RegexOptions.None, TimeSpan.FromMilliseconds(150));
            //Match m = r.Match(url);
            //if (m.Success)
            //    Console.WriteLine(m.Result("${proto}:${port}"));

            // replacement for drem(d1, d2)
            //Console.WriteLine(-5.2f % 2.0f); // output: -1.2
            //Console.WriteLine(5.0 % 2.2);    // output: 0.6
            //Console.WriteLine(.41 % .2);     // output: 0.00999999999999995
#endif

            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new Form1());
        }
    }
}
