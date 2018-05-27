using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace CourseworkAlgo2
{
    public class ProblemCalculator
    {
        private readonly ProblemData _problemData;

        public ProblemCalculator(ProblemData problemData)
        {
            _problemData = problemData;
        }

        public Complex[] CalculateEigenValues(string filename)
        {
            var s0 = CalculateEigenvalueAmount();
            Console.WriteLine(s0);
            Logger.WriteToFile(this._problemData, s0.ToString(), filename);
            var m = (int)Math.Floor(s0.Real);
            var result = new Complex[m];
            for (int j = 0; j < m; j++)
            {
                var arg = 2 * Math.PI * (j + 1) / m;
                var c1 = _problemData.Center + _problemData.Radius * new Complex(Math.Cos(arg), Math.Sin(arg));
                result[j] = c1;
            }

            return result;
        }

        public Complex[] PreciseEigenValues(Complex[] eigenValues, Complex[] sk)
        {
            var f = GetFVector(eigenValues);
            var jacobian = GetJacobianMatrix(eigenValues);
            var matrix = Matrix.Build.DenseOfRowArrays(jacobian);
            var matrixInverse = matrix.Inverse();
            Console.WriteLine(matrixInverse);
            var inverseJacobian = matrixInverse.ToRowArrays();
            inverseJacobian.ConsoleWrite();

            var fs = f.Subtract(sk);
            var jacobianFs = inverseJacobian.Multiply(fs);

            return eigenValues.Subtract(jacobianFs);
        }

        private Complex[][] GetJacobianMatrix(Complex[] eigenValues)
        {
            var n = eigenValues.Length;
            var jacobianMatrix = new Complex[n][];
            for (int i = 0; i < n; i++)
            {
                jacobianMatrix[i] = new Complex[n];
                for (int j = 0; j < n; j++)
                {
                    jacobianMatrix[i][j] = (i + 1) * Complex.Pow(eigenValues[j], i);
                }
            }

            return jacobianMatrix;
        }

        private Complex[] GetFVector(Complex[] eigenValues)
        {
            var f = new Complex[eigenValues.Length];
            for (int i = 0; i < f.Length; i++)
            {
                f[i] = eigenValues.Select(x => Complex.Pow(x, i + 1)).Aggregate((x, y) => x + y);
            }
            return f;
        }

        public Complex CalculateEigenvalueAmount(int power = 0)
        {
            int n = 4;

            Complex sum = Complex.Zero;
            for (var j = 1; j <= n; j++)
            {
                var arg = 2 * Math.PI * j / n;
                var complexMultiplier = new Complex(Math.Cos(arg), Math.Sin(arg));
                var lambdaJ = _problemData.Center + _problemData.Radius * complexMultiplier;
                var d = GetDMatrix(lambdaJ);
                var b = GetBMatrix(lambdaJ);
                var (u, v) = MatrixDecompositor.DecomposeMatrixDiagonals(d, b);
                Complex subSum = 0.0;
                for (int r = 0; r < u.Length; r++)
                {
                    subSum += v[r] / u[r];
                }

                sum += Complex.Pow(lambdaJ, power) * _problemData.Radius * complexMultiplier * subSum;
                Console.WriteLine(j);
            }

            return sum / n;
        }

        private Complex[][] GetDMatrix(Complex lambda)
        {
            var d = new Complex[_problemData.Ksi1.PartitionsAmount][];
            for (var i = 0; i < _problemData.Ksi1.PartitionsAmount; i++)
            {
                d[i] = new Complex[_problemData.Ksi2.PartitionsAmount];
                for (var j = 0; j < _problemData.Ksi2.PartitionsAmount; j++)
                {
                    var i1 = i;
                    var j1 = j;

                    Complex Func(double ksi1, double ksi2) =>
                        Complex.Sqrt(_problemData.P(ksi1, ksi2) - _problemData.LambdaValue) *
                        GetKFunction(
                            _problemData.Ksi1.GetKsiForPartition(i1),
                            _problemData.Ksi2.GetKsiForPartition(j1),
                            lambda)(ksi1, ksi2);

                    d[i][j] = IntegralCalculator.ComputeIntegralForProblem(Func, _problemData)
                              * Complex.Sqrt(_problemData.P(
                                  _problemData.Ksi1.GetKsiForPartition(i1),
                                  _problemData.Ksi2.GetKsiForPartition(j1))
                                             - _problemData.LambdaValue);
                }
            }

            for (int i = 0; i < d.Length; i++)
            {
                d[i][i] -= 1;
            }

            return d;
        }

        private Complex[][] GetBMatrix(Complex lambda)
        {
            var b = new Complex[_problemData.Ksi1.PartitionsAmount][];
            for (var i = 0; i < _problemData.Ksi1.PartitionsAmount; i++)
            {
                b[i] = new Complex[_problemData.Ksi2.PartitionsAmount];
                for (var j = 0; j < _problemData.Ksi2.PartitionsAmount; j++)
                {
                    var i1 = i;
                    var j1 = j;

                    Complex Func(double ksi1, double ksi2) =>
                        Complex.Sqrt(_problemData.P(ksi1, ksi2) - _problemData.LambdaValue) *
                        GetKFunctionDerivativeC1(
                            _problemData.Ksi1.GetKsiForPartition(i1),
                            _problemData.Ksi2.GetKsiForPartition(j1),
                            lambda)(ksi1, ksi2);

                    b[i][j] = IntegralCalculator.ComputeIntegralForProblem(Func, _problemData)
                              * Complex.Sqrt(_problemData.P(
                                  _problemData.Ksi1.GetKsiForPartition(i1),
                                  _problemData.Ksi2.GetKsiForPartition(j1))
                                             - _problemData.LambdaValue);
                }
            }

            return b;
        }

        private Func<double, double, Complex> GetKFunction(double inKsi1, double inKsi2, Complex lambda)
        {
            return (ksi1, ksi2) =>
            {
                var result = new Complex();
                for (var n = -_problemData.N; n <= _problemData.N; n++)
                {
                    for (var m = -_problemData.M; m <= _problemData.M; m++)
                    {
                        var angle = lambda * n * (inKsi1 - ksi1) + _problemData.GetC2(lambda) * m * (inKsi2 - ksi2);
                        result += Complex.Cos(angle) + Complex.ImaginaryOne * Complex.Sin(angle);
                    }

                }

                return result;
            };
        }

        private Func<double, double, Complex> GetKFunctionDerivativeC1(double inKsi1, double inKsi2, Complex lambda)
        {
            return (ksi1, ksi2) =>
            {
                var result = new Complex();
                for (var n = -_problemData.N; n <= _problemData.N; n++)
                {
                    for (var m = -_problemData.M; m <= _problemData.M; m++)
                    {
                        var angle = lambda * n * (inKsi1 - ksi1) + _problemData.GetC2(lambda) * m * (inKsi2 - ksi2);
                        result += n * Complex.ImaginaryOne * (inKsi1 - ksi1) * (Complex.Cos(angle) + Complex.ImaginaryOne * Complex.Sin(angle));
                    }

                }

                return result;
            };
        }
    }
}
