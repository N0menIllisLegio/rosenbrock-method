using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
    /// <summary>
    /// Rosenbrock method. Numerical method for finding unconditional extremum.
    /// </summary>
    public sealed class RosenbrockMethod
    {
        private static Data input;

        /// <summary>
        /// Size of vectors.
        /// </summary>
        private static int n;
        private static double beta;
        private static double alpha;
        private static double epsilon;
        private static int N;
        private static int l;

        /// <summary>
        /// Initial step distance in every search direction.
        /// </summary>
        private static Vector<double> delta0;

        /// <summary>
        /// Gets or sets the stretch ratio (α &gt; 1).
        /// Recommended: 3.
        /// </summary>
        /// <value>Stretch ratio.</value>
        public double Alpha
        {
            set
            {
                if (value <= 1)
                {
                    Console.WriteLine("alpha > 1");
                }
                else
                {
                    alpha = value;
                }
            }
            get { return alpha; }
        }

        /// <summary>
        /// Gets or sets the compression ratio (-1 &lt; β &lt; 0).
        /// Recommended: -0.5;
        /// </summary>
        /// <value>Compression ratio.</value>
        public double Beta
        {
            set
            {
                if (value > 0 || value < -1)
                {
                    Console.WriteLine("-1 < beta < 0");
                }
                else
                {
                    beta = value;
                }
            }
            get { return beta; }
        }

        /// <summary>
        /// Gets or sets the number to stop algorithm (ε &gt; 0).
        /// </summary>
        /// <value>Precision.</value>
        public double Epsilon
        {
            set
            {
                if (value <= 0)
                {
                    Console.WriteLine("epsilon > 0");
                }
                else
                {
                    epsilon = value;
                }
            }
            get { return epsilon; }
        }

        /// <summary>
        /// Gets or sets the maxaximum amount of failed steps on one iteration.
        /// </summary>
        /// <value>Maxaximum amount of failed steps.</value>
        public int MaxFailedSteps
        {
            set
            {
                if (value < 1)
                {
                    Console.WriteLine("MaxMistakes > 0");
                }
                else
                {
                    N = value;
                }
            }
            get { return N; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="T:RosenbrockMethod.RosenbrockMethod"/> class.
        /// </summary>
        /// <param name="input">Input data, witch contains yours function for minimization and constrains for it (if they are needed).</param>
        public RosenbrockMethod(Data input)
        {
            n = input.N - 1;
            RosenbrockMethod.input = input;

            delta0 = Vector<double>.Build.Dense(input.N, 1.0);
            beta = -0.5;
            alpha = 3;
            epsilon = 0.00000075;
            N = 1000;
        }

        /// <summary>
        /// Minimizes given function, using Rosenbrock's method.
        /// </summary>
        /// <returns>Vector of minimum points.</returns>
        public Vector<double> MinimizeFunction()
        {
            Matrix<double> d = Matrix<double>.Build.DenseDiagonal(n + 1, n + 1, 1.0);
            List<Vector<double>> x = new List<Vector<double>>{ Vector<double>.Build.Dense(n + 1, 1.0) };
            List<Vector<double>> y = new List<Vector<double>>{ x[0] };

            Vector<double> lambda = Vector<double>.Build.Dense(n + 1);
            Vector<double> delta = delta0.Clone();

            int iteration = 0;
            int currStep = 2;
            int i = 0;

            while (currStep != 0)
            {
                if (currStep == 2)
                {
                    currStep = Step2(ref delta, ref y, ref d, i);
                }
                if (currStep == 3)
                {
                    currStep = Step3(ref delta, ref y, ref x, ref i, iteration);
                }
                if (currStep == 4)
                {
                    currStep = Step4(ref delta, ref y, ref x, ref lambda, ref d, ref i, ref iteration);
                }
            }

            return x[iteration];
        }

        private static int Step2(ref Vector<double> delta, ref List<Vector<double>> y, ref Matrix<double> d, int i)
        {
            Vector<double> temp = y[i].Add(d.Row(i).Multiply(delta[i]));
            double result = input.Function(temp);

            if (result < input.Function(y[i]) && input.CheckConstraints(temp))
            {
                l = 0;
                y.Add(temp);
                delta[i] *= alpha;
            }
            else
            {
                l++;
                y.Add(y[i]);
                delta[i] *= beta;
            }

            return 3;
        }

        private static int Step3(ref Vector<double> delta, ref List<Vector<double>> y, ref List<Vector<double>> x, ref int i, int iteration)
        {
            if (i < n)
            {
                i++;
                return 2;
            }

            if (i == n)
            {
                if (input.Function(y[n + 1]) < input.Function(y[0]) && input.CheckConstraints(y[n + 1]))
                {
                    Vector<double> temp = y[n + 1];
                    y.Clear();
                    y.Add(temp);
                    i = 0;

                    return 2;
                }

                if (input.Function(y[n + 1]).Equals(input.Function(y[0])))
                {
                    if (input.Function(y[n + 1]) < input.Function(x[iteration]) && input.CheckConstraints(y[n + 1]))
                    {
                        return 4;
                    }

                    //additional constraint??
                    if (input.Function(y[n + 1]).Equals(input.Function(x[iteration])) && input.CheckConstraints(y[n + 1]))
                    {
                        if (l <= N)
                        {
                            if (CheckEndCondition(ref delta))
                            {
                                return 0;
                            }
                            else
                            {
                                Vector<double> temp = y[n + 1];
                                y.Clear();
                                y.Add(temp);
                                i = 0;

                                return 2;
                            }
                        }
                    }
                }
            }

            return 4;
        }

        private static int Step4(ref Vector<double> delta, ref List<Vector<double>> y, ref List<Vector<double>> x, ref Vector<double> lambda, ref Matrix<double> d, ref int i, ref int iteration)
        {
            x.Add(y[n + 1]);
            iteration++;

            if (CheckEndCondition(ref x, ref lambda, iteration))
            {
                return 0;
            }
            else
            {
                lambda = d.Solve(lambda);
                ProcedureGrammaShmita(ref lambda, ref d);
                delta = delta0.Clone();

                y.Clear();

                y.Add(x[iteration]);

                i = 0;
                l = 0;

                return 2;
            }
        }

        private static bool CheckEndCondition(ref List<Vector<double>> x, ref Vector<double> lambda, int iteration)
        {
            lambda = x[iteration].Subtract(x[iteration - 1]);

            return lambda.L2Norm() <= epsilon;
        }

        private static bool CheckEndCondition(ref Vector<double> delta)
        {
            foreach (var item in delta)
            {
                if (Math.Abs(item) > epsilon)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Gram–Schmidt process.
        /// </summary>
        private static void ProcedureGrammaShmita(ref Vector<double> lambda, ref Matrix<double> d)
        {
            List<Vector<double>> a = new List<Vector<double>>();
            List<Vector<double>> b = new List<Vector<double>>();

            for (int i = 0; i < n + 1; i++)
            {
                if (lambda[i].Equals(0))
                {
                    a.Add(d.Row(i));
                }
                else
                {
                    Vector<double> temp = Vector<double>.Build.Dense(n + 1, 0.0);

                    for (int j = i; j < n + 1; j++)
                    {
                        temp = temp.Add(d.Row(j).Multiply(lambda[j]));
                    }

                    a.Add(temp);
                }
            }

            for (int i = 0; i < n + 1; i++)
            {
                if (i != 0)
                {
                    Vector<double> temp = Vector<double>.Build.Sparse(n + 1);

                    for (int j = 0; j < i; j++)
                    {
                        temp = temp.Add(d.Row(j).Multiply(a[i].ConjugateDotProduct(d.Row(j))));
                    }

                    temp = a[i].Subtract(temp);
                    b.Add(temp);
                }
                else
                {
                    b.Add(a[i]);
                }

                d.SetRow(i,b[i].Divide(b[i].L2Norm()));
            }
        }
    }
}