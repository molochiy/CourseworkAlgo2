using System;
using System.Linq;
using System.Numerics;

namespace CourseworkAlgo2
{
    public static class MatrixDecompositor
    {
        public static (Complex[][], Complex[][], Complex[][], Complex[][]) DecomposeMatrix(Complex[][] d, Complex[][] b)
        {
            var n = d.Length;
            var l = new Complex[n][];
            var u = new Complex[n][];
            var m = new Complex[n][];
            var v = new Complex[n][];
            for (var i = 0; i < n; i++)
            {
                l[i] = new Complex[n];
                l[i][i] = 1;
                u[i] = new Complex[n];
                m[i] = new Complex[n];
                v[i] = new Complex[n];
            }

            for (int r = 0; r < n; r++)
            {
                for (int k = r; k < n; k++)
                {
                    Complex sum = 0.0;
                    for (int j = 0; j < r; j++)
                    {
                        sum += l[r][j] * u[j][k];
                    }

                    u[r][k] = d[r][k] - sum;
                }

                for (int i = r + 1; i < n; i++)
                {
                    Complex sum = 0.0;
                    for (int j = 0; j < r; j++)
                    {
                        sum += l[i][j] * u[j][r];
                    }

                    l[i][r] = (d[i][r] - sum) / u[r][r];
                }

                for (int k = r; k < n; k++)
                {
                    Complex sum = 0.0;
                    for (int j = 0; j < r; j++)
                    {
                        sum += m[r][j] * u[j][k] + l[r][j]*v[j][k];
                    }

                    v[r][k] = b[r][k] - sum;
                }

                for (int i = r + 1; i < n; i++)
                {
                    Complex sum = 0.0;
                    for (int j = 0; j < r; j++)
                    {
                        sum += m[i][j] * u[j][r] + l[i][j] * v[j][r];
                    }

                    m[i][r] = (b[i][r] - sum - l[i][r] * v[r][r]) / u[r][r];
                }
            }

            return (l, u, m, v);
        }

        public static (Complex[], Complex[]) DecomposeMatrixDiagonals(Complex[][] d, Complex[][] b)
        {
            var (l, u, m, v) = DecomposeMatrix(d, b);

            var uDiagonal = u.SelectMany((row, rowIndex) => row.Where((cell, columnIndex) => columnIndex == rowIndex)).ToArray();
            var vDiagonal = v.SelectMany((row, rowIndex) => row.Where((cell, columnIndex) => columnIndex == rowIndex)).ToArray();

            return (uDiagonal, vDiagonal);
        }

        public static void DecomposeMatrixTest()
        {
            var d = new[]
            {
                new Complex[] { 1, 2, 1 },
                new Complex[] { 2, 1, 1 },
                new Complex[] { 1, -1, 2 }
            };
            var b = new[]
            {
                new Complex[] { 3, 2, 6 },
                new Complex[] { 1, 8, 5 },
                new Complex[] { 5, 1, 4 }
            };
            var (l, u, m, v) = DecomposeMatrix(d, b);
            d.ConsoleWrite();
            Console.WriteLine();
            l.ConsoleWrite();
            Console.WriteLine();
            u.ConsoleWrite();
            Console.WriteLine();

            l.Multiply(u).ConsoleWrite();
            Console.WriteLine();

            b.ConsoleWrite();
            Console.WriteLine();
            v.ConsoleWrite();
            Console.WriteLine();
            m.ConsoleWrite();
            Console.WriteLine();
            m.Multiply(u).Add(l.Multiply(v)).ConsoleWrite();
            Console.WriteLine();

            var (uDiag, vDiag) = DecomposeMatrixDiagonals(d, b);
            foreach (var el in uDiag)
            {
                Console.Write($"{el} ");
            }
            Console.WriteLine();

            foreach (var el in vDiag)
            {
                Console.Write($"{el} ");
            }
            Console.WriteLine();
        }
    }
}
