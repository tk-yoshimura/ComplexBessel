using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class ConnectionTests {
        
        readonly Complex<Pow2.N4>[] zs = [
            (1, 0), (1, 0.5), (2, 0), (2, 0.5), (2, 1), (0, 1), (0.5, 1), (0, 2), (0.5, 2), (1, 2)
        ];

        [TestMethod()]
        public void BesselItoJTest() {
            for (double nu = -8; nu <= 8; nu += 0.125) {
                Console.WriteLine(nu);

                Complex<Pow2.N4> c = (MultiPrecision<Pow2.N4>.CosPI(nu / 2), MultiPrecision<Pow2.N4>.SinPI(nu / 2));

                foreach (Complex<Pow2.N4> z in zs) {
                    Complex<Pow2.N4> bj = PowerSeries<Pow2.N4>.BesselJ(nu, z);
                    Complex<Pow2.N4> bi = PowerSeries<Pow2.N4>.BesselI(nu, (z.I, z.R));

                    Complex<Pow2.N4> y = c * bi.Conj;

                    Console.WriteLine($"{z}\n {bj}\n {y}");

                    Assert.IsTrue((y - bj).Magnitude / bj.Magnitude < 1e-30);
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselIKtoYTest() {
            for (double nu = -8; nu <= 8; nu += 0.125) {
                Console.WriteLine(nu);

                Complex<Pow2.N4> c = (MultiPrecision<Pow2.N4>.CosPI(nu / 2), MultiPrecision<Pow2.N4>.SinPI(nu / 2));

                foreach (Complex<Pow2.N4> z in zs) {
                    Complex<Pow2.N4> by = PowerSeries<Pow2.N4>.BesselY(nu, z);
                    Complex<Pow2.N4> bi = PowerSeries<Pow2.N4>.BesselI(nu, (z.I, z.R));
                    Complex<Pow2.N4> bk = PowerSeries<Pow2.N4>.BesselK(MultiPrecision<Pow2.N4>.Abs(nu), (z.I, z.R));

                    Complex<Pow2.N4> y = (0, 1) * c * bi.Conj - 2 * MultiPrecision<Pow2.N4>.RcpPI * (c * bk).Conj;

                    Console.WriteLine($"{z}\n {by}\n {y}");

                    if (nu - double.Floor(nu) != 0.5) {
                        Assert.IsTrue((y - by).Magnitude / by.Magnitude < 1e-30);
                    }
                    else { 
                        Assert.IsTrue((y - by).Magnitude / by.Magnitude < 1e-20);
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselIYtoKTest() {
            for (double nu = -8; nu <= 8; nu += 0.125) {
                Console.WriteLine(nu);

                Complex<Pow2.N4> c = (MultiPrecision<Pow2.N4>.CosPI(nu / 2), MultiPrecision<Pow2.N4>.SinPI(nu / 2));

                foreach (Complex<Pow2.N4> z in zs) {
                    Complex<Pow2.N4> by = PowerSeries<Pow2.N4>.BesselY(nu, (z.I, z.R));
                    Complex<Pow2.N4> bi = PowerSeries<Pow2.N4>.BesselI(nu, z);
                    Complex<Pow2.N4> bk = PowerSeries<Pow2.N4>.BesselK(MultiPrecision<Pow2.N4>.Abs(nu), z);

                    Complex<Pow2.N4> y = (Complex<Pow2.N4>.Pow((0, 1), - 2 * nu - 1) * bi - (c * by).Conj) * MultiPrecision<Pow2.N4>.PI / 2;

                    Console.WriteLine($"{z}\n {bk}\n {y}");

                    if (nu - double.Floor(nu) != 0.5) {
                        Assert.IsTrue((y - bk).Magnitude / bk.Magnitude < 1e-30);
                    }
                    else { 
                        Assert.IsTrue((y - bk).Magnitude / bk.Magnitude < 1e-20);
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }
    }
}