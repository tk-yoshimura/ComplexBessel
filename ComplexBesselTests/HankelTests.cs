using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class HankelTests {
        [TestMethod()]
        public void Hankel1Test() {
            double nu = 2.5;

            for (double r = 0; r <= 4; r++) {
                for (double i = 0; i <= 4; i++) {
                    Complex<Pow2.N4> z1 = BesselN4.BesselK(nu, -Complex<Pow2.N4>.ImaginaryOne * (r, i));
                    Complex<Pow2.N4> z2 = BesselN4.BesselK(nu, (i, r)).Conj;

                    Console.WriteLine($"{r},{i}");
                    Console.WriteLine($"{z1}");
                    Console.WriteLine($"{z2}");
                    Console.WriteLine($"{z1 - z2}\n");
                }
            }

            for (double r = -4; r <= 4; r++) {
                for (double i = 0; i <= 4; i++) {
                    Complex<Pow2.N4> c = (-MultiPrecision<Pow2.N4>.SinPI(nu / 2), -MultiPrecision<Pow2.N4>.CosPI(nu / 2));
                    Complex<Pow2.N4> z1 = 2 * MultiPrecision<Pow2.N4>.RcpPI * c * BesselN4.BesselK(nu, (i, r)).Conj;
                    Complex<Pow2.N4> z2 = BesselN4.BesselJ(nu, (r, i)) + (0, 1) * BesselN4.BesselY(nu, (r, i));

                    Console.WriteLine($"{r},{i}");
                    Console.WriteLine($"{z1}");
                    Console.WriteLine($"{z2}");
                    Console.WriteLine($"{z1 - z2}\n");
                }
            }
        }

        [TestMethod()]
        public void Hankel2Test() {
            double nu = 2.5;

            for (double r = 0; r <= 4; r++) {
                for (double i = 0; i <= 4; i++) {
                    Complex<Pow2.N4> z1 = BesselN4.BesselK(nu, Complex<Pow2.N4>.ImaginaryOne * (r, i));
                    Complex<Pow2.N4> z2 = BesselN4.BesselK(nu, (-i, r));

                    Console.WriteLine($"{r},{i}");
                    Console.WriteLine($"{z1}");
                    Console.WriteLine($"{z2}");
                    Console.WriteLine($"{z1 - z2}\n");
                }
            }

            for (double r = -4; r <= 4; r++) {
                for (double i = 0; i <= 4; i++) {
                    Complex<Pow2.N4> c = (-MultiPrecision<Pow2.N4>.SinPI(nu / 2), MultiPrecision<Pow2.N4>.CosPI(nu / 2));
                    Complex<Pow2.N4> z1 = 2 * MultiPrecision<Pow2.N4>.RcpPI * c *
                        (double.IsPositive(r) ? BesselN4.BesselK(nu, (-i, r)) : -BesselN4.BesselK(nu, (-i, double.Abs(r))).Conj);
                    Complex<Pow2.N4> z2 = BesselN4.BesselJ(nu, (r, i)) - (0, 1) * BesselN4.BesselY(nu, (r, i));

                    Console.WriteLine($"{r},{i}");
                    Console.WriteLine($"{z1}");
                    Console.WriteLine($"{z2}");
                    Console.WriteLine($"{z1 - z2}\n");
                }
            }
        }
    }
}