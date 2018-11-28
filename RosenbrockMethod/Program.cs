using System;
using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
    class MainClass
    {
        public static void Main()
        {
            Console.WriteLine("Without area limitaion:");
            DataWithoutAreaLimitation unlimitData = new DataWithoutAreaLimitation();
            RosenbrockMethod method1 = new RosenbrockMethod(unlimitData);
            Vector<double> temp = method1.MinimizeFunction();

            Console.WriteLine("x*:");
            Console.WriteLine(temp);
            Console.WriteLine("F = " + unlimitData.GetArea(temp));
            Console.WriteLine("L = " + unlimitData.Function(temp));

            Console.WriteLine("\nWith area limitaion:");
            DataWithAreaLimitaion limitData = new DataWithAreaLimitaion();
            RosenbrockMethod method2 = new RosenbrockMethod(limitData);
            temp = method2.MinimizeFunction();
            
            Console.WriteLine("x*:");
            Console.WriteLine(temp);
            Console.WriteLine("F = " + limitData.GetArea(temp));
            Console.WriteLine("L = " + limitData.Function(temp));
        }
    }
}
