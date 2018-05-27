using System;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Reflection;

namespace CourseworkAlgo2
{
    public static class Logger
    {
        public static void WriteToFile(ProblemData problemData, string text, string fileName)
        {
            var path = Directory.GetParent(Directory.GetParent(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location)).FullName).FullName;
            var file = new FileInfo($"{path}\\results\\{fileName}");
            file.Directory?.Create();

            using (var writer = new StreamWriter(file.FullName, true))
            {
                writer.WriteLine(text);
                writer.Close();
            }
        }

        public static void WriteIterationToFile(ProblemData problemData, Complex [] values, int iteration, string fileName)
        {
            var path = Directory.GetParent(Directory.GetParent(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location)).FullName).FullName;
            var file = new FileInfo($"{path}\\results\\{fileName}");
            file.Directory?.Create();

            using (var writer = new StreamWriter(file.FullName, true))
            {
                writer.WriteLine($"Iteration #{iteration}");
                var complexC1 = $"[{string.Join(", ", values)}],";
                var complexC2 = $"[{string.Join(", ", values.Select(problemData.GetC2))}],";

                var c1 = $"[{string.Join(", ", values.Select(v => v.Magnitude))}],";
                var c2 = $"[{string.Join(", ", values.Select(v => problemData.GetC2(v).Magnitude))}],";
                writer.WriteLine($"alpha {problemData.Coef1}");
                writer.WriteLine(complexC1);
                writer.WriteLine(complexC2);
                writer.WriteLine();
                writer.WriteLine(c1);
                writer.WriteLine(c2);
                /*foreach (var eigenValue in values)
                {
                    {
                        writer.WriteLine($"C1 {eigenValue.c1}, C2 {eigenValue.c2}");
                    }
                }*/
                writer.WriteLine();
                writer.Close();
            }
        }

        public static void WriteResults(ProblemData problemData, (Complex c1, Complex c2)[] values, DateTime time)
        {
            var path = Directory.GetParent(Directory.GetParent(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location)).FullName).FullName;
            var file = new FileInfo($"{path}\\results\\{problemData.Coef1}_{problemData.Coef2}\\result_{time:yyyy-MM-dd_hh-mm-ss-fff}.txt");
            file.Directory?.Create();
            using (var writer = new StreamWriter(file.FullName))
            {
                foreach (var eigenValue in values)
                {
                    {
                        writer.WriteLine($"{eigenValue.c1.Magnitude}, {eigenValue.c2.Magnitude}");
                    }
                }
                writer.WriteLine();
                writer.Close();
            }
        }
    }
}
