using System;
using System.Numerics;

namespace CourseworkAlgo2
{
    class Runner
    {
        public static void Run()
        {
            //for (int i = 0; i <= 3; i++)
            //{
            //    var runTime = DateTime.Now;
            //    for (double j = 0.1; j <= 20; j += j < 1 ? 0.1 : 1)
            //    {
            //        RunForProblemData(new ProblemData
            //        {
            //            Coef1 = j,
            //            LambdaLimit = (i, i + 1)
            //        }, runTime);
            //    }
            //}

            var runTime = DateTime.Now;
            for (double j = 0.1; j <= 1; j += 0.1)
            {
                RunForProblemData(new ProblemData
                {
                    Coef1 = j,
                    LambdaLimit = (0.0001, 1.5),
                    Ksi1 = new KsiData
                    {
                        Step = 0.2
                    },
                    Ksi2 = new KsiData
                    {
                        Step = 0.2
                    }
                }, runTime);
            }
        }

        private static void RunForProblemData(ProblemData problemData, DateTime runTime)
        {
            try
            {
                Console.WriteLine($"alpha {problemData.Coef1}.");

                var problemCalculator = new ProblemCalculator(problemData);

                int iteration = 0;
                var iterationsFileName = $"Iterations_{runTime:yyyy-MM-dd_hh-mm-ss-fff}.txt";

                var prevEigenValues = problemCalculator.CalculateEigenValues(iterationsFileName);
                var sk = new Complex[prevEigenValues.Length];
                for (int i = 0; i < sk.Length; i++)
                {
                    sk[i] = problemCalculator.CalculateEigenvalueAmount(i + 1);
                }
                //var prevEigenValues = new Complex[] { 0.5, 0.5, 0.5 };
                //var sk = new Complex[] { 0, 0, 0 };
                Console.WriteLine("Sk");
                sk.ConsoleWrite();

                Console.WriteLine("Eigen Values");
                prevEigenValues.ConsoleWrite();
                Logger.WriteIterationToFile(problemData, prevEigenValues, iteration, iterationsFileName);

                var nextEigenValues = problemCalculator.PreciseEigenValues(prevEigenValues, sk);
                nextEigenValues.ConsoleWrite();
                Logger.WriteIterationToFile(problemData, nextEigenValues, ++iteration, iterationsFileName);

                while (!IsSatisfyPrec(prevEigenValues, nextEigenValues, /*problemData.Prec*/ 0.001) && iteration < 0.5e3)
                {
                    Console.WriteLine($"Iteration: {++iteration}");

                    prevEigenValues = nextEigenValues;
                    nextEigenValues = problemCalculator.PreciseEigenValues(prevEigenValues, sk);
                    nextEigenValues.ConsoleWrite();

                    Logger.WriteIterationToFile(problemData, nextEigenValues, iteration, iterationsFileName);
                }

                Console.WriteLine($"alpha {problemData.Coef1} finished.");
            }
            catch
            {
                // ignored
            }
        }

        private static bool IsSatisfyPrec(Complex[] prev, Complex[] next, double prec)
        {
            for (int i = 0; i < prev.Length; i++)
            {
                if (Math.Abs(prev[i].Magnitude - next[i].Magnitude) > prec)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
