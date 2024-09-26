using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class YoshidaPadeTests {
        readonly Complex<Pow2.N4>[] zs = [
            (64,0.0),(64,1),(64,2),(64,4),(64,8),(64,16),(64,32),(64,64),
            (0.0,64),(1,64),(2,64),(4,64),(8,64),(16,64),(32,64),(64,64),
        ];

        readonly Complex<Pow2.N4>[] zs_mini = [
            (8, 0), (8, 2), (0, 8), (2, 8)
        ];

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                YoshidaPade<Pow2.N4> pade = new(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {pade.BesselK(z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKNu1p25Test() {
            YoshidaPade<Pow2.N4> pade = new(1.25);

            Complex<Pow2.N4>[] expecteds = [
                "0.0001606019006214995920875098806847954111049",
                "-0.0000848796993628445379920041646286883097161-0.0001327193785366760545268043058004086379747i",
                "-0.3778461312343345326806188610167906742794-0.2357537784838237984707580183773611438572i",
                "-0.04673236980325523574728826588564566867056-0.03835216190634850921510516245573954493834i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = pade.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}