using MultiPrecision;
using MultiPrecisionComplex;

namespace BesselWorkbench {
    internal class Program {
        static void Main(string[] args) {
            HankelExpansion<Pow2.N4> bessel_asymp = new(1.25);

            Complex<Pow2.N4> c = bessel_asymp.BesselY((-128, 1));

            Console.WriteLine(c);

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
