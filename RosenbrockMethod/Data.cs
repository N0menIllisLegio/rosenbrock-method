using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
    /// <summary>
    /// Input data template for Rosenbrock's method.
    /// </summary>
    public abstract class Data
    {
        public abstract int N { get; }

        /// <summary>
        /// Check if the specified vector fit yours constraints.
        /// </summary>
        /// <returns><see langword="true"/> - if <paramref name="vector"/> fit constraints and <see langword="false"/> otherwise.</returns>
        /// <param name="vector">Vector to check.</param>
        public virtual bool CheckConstraints(Vector<double> vector)
        {
            return true;
        }

        /// <summary>
        /// Function for minimization.
        /// </summary>
        /// <returns>Result of function.</returns>
        /// <param name="parameters">Input parameters for function.</param>
        public abstract double Function(Vector<double> parameters);
    }
}
