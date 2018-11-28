using System;
using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
    public class DataWithoutAreaLimitation : Data
    {
        readonly int Size = 5;

        readonly int[] V = { 3000, 5000, 6400, 1500, 80 };
        readonly int[] K = { 4, 6, 7, 6, 4 };
        readonly int[] S = { 40, 6, 14, 6, 16 };
        readonly int[] F = { 4, 3, 5, 40, 20 };


        public override int N { get { return Size; } }

        public override double Function(Vector<double> parameters)
        {
            double[] values = parameters.AsArray();
            double result = 0.0;

            for (int i = 0; i < values.Length; i++)
            {
                result += K[i] * V[i] / values[i] + 1.0 / 2.0 * S[i] * values[i];
            }

            return result;
        }

        private bool CheckPositive(Vector<double> vector)
        {
            foreach (var item in vector)
            {
                if (item < 0)
                {
                    return false;
                }
            }

            return true;
        }

        public double GetArea(Vector<double> vector)
        {
            double[] values = vector.AsArray();
            double result = 0.0;

            for (int i = 0; i < values.Length; i++)
            {
                result += F[i] * values[i];
            }

            return result;
        }
    }
}
