using System;
using System.Numerics;

namespace CourseworkAlgo2
{
    public static class MatrixExtensions
    {
        public static void ConsoleWrite(this Complex[][] matrix)
        {
            foreach (var matrixRow in matrix)
            {
                foreach (var cell in matrixRow)
                {
                    Console.Write($"({cell.Real:0.##}, {cell.Imaginary:0.##})");
                }

                Console.WriteLine();
            }
        }

        public static void ConsoleWrite(this Complex[] vector)
        {
            foreach (var cell in vector)
            {
                Console.Write($"({cell.Real:0.##}, {cell.Imaginary:0.##}) ");
            }
            Console.WriteLine();
        }

        public static Complex[][] Multiply(this Complex[][] matrix1, Complex[][] matrix2)
        {
            var n = matrix1.Length;
            var result = new Complex[n][];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new Complex[n];
                for (int j = 0; j < result[i].Length; j++)
                {
                    result[i][j] = 0;
                    for (int k = 0; k < matrix1.Length; k++)
                    {
                        result[i][j] += matrix1[i][k] * matrix2[k][j];
                    }
                }
            }

            return result;
        }

        public static Complex[][] Add(this Complex[][] matrix1, Complex[][] matrix2)
        {
            var n = matrix1.Length;
            var result = new Complex[n][];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new Complex[n];
                for (int j = 0; j < result[i].Length; j++)
                {
                    result[i][j] = matrix1[i][j] + matrix2[i][j];
                }
            }

            return result;
        }

        public static Complex[][] Subtract(this Complex[][] matrix1, Complex[][] matrix2)
        {
            var n = matrix1.Length;
            var result = new Complex[n][];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new Complex[n];
                for (int j = 0; j < result[i].Length; j++)
                {
                    result[i][j] = matrix1[i][j] - matrix2[i][j];
                }
            }

            return result;
        }

        public static Complex[] Multiply(this Complex[][] matrix, Complex[] vector)
        {
            var n = matrix.Length;
            var result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                result[i] = 0;
                for (int j = 0; j < n; j++)
                {

                    result[i] += matrix[i][j] * vector[j];
                }
            }

            return result;
        }

        public static Complex[] Add(this Complex[] vector1, Complex[] vector2)
        {
            var n = vector1.Length;
            var result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                result[i] = vector1[i] + vector2[i];
            }

            return result;
        }

        public static Complex[] Subtract(this Complex[] vector1, Complex[] vector2)
        {
            var n = vector1.Length;
            var result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                result[i] = vector1[i] - vector2[i];
            }

            return result;
        }
    }
}
