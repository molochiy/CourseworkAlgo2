using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace CourseworkAlgo2
{
    public class ProblemData
    {
        private double? _alpha;
        private (double begin, double end)? _lambdaLimit;

        public int M { get; set; } = 5;
        public int N { get; set; } = 5;

        public KsiData Ksi1 { get; set; } = new KsiData();
        public KsiData Ksi2 { get; set; } = new KsiData();

        public double Prec { get; set; } = 1e-6;

        public Func<double, double, Complex> P { get; set; } = (ksi1, ksi2) => 1;

        public double Coef1 { get; set; } = 0.2;
        public double Coef2 { get; set; } = 0;

        public (double begin, double end) LambdaLimit
        {
            get => _lambdaLimit ?? (0.5, 2.5);
            set => _lambdaLimit = value;
        }

        public double Center => (LambdaLimit.begin + LambdaLimit.end) / 2;
        public double Radius => (LambdaLimit.end - LambdaLimit.begin) / 2;

        public Complex LambdaValue { get; set; } = 0.0;

        public Complex GetC2(Complex c1) => Coef1 * c1 + Coef2;
    }
}
