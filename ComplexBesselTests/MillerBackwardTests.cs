using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class MillerBackwardTests {
        readonly Complex<Pow2.N4>[] zs = [
            (64,0.0),(64,1),(64,2),(64,4),(64,8),(64,16),(64,32),(64,64),
            (0.0,64),(1,64),(2,64),(4,64),(8,64),(16,64),(32,64),(64,64),
        ];

        readonly Complex<Pow2.N4>[] zs_mini = [
            (8, 0), (8, 2), (0, 8), (2, 8)
        ];

        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {MillerBackward<Pow2.N4>.BesselJ(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {MillerBackward<Pow2.N4>.BesselY(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {MillerBackward<Pow2.N4>.BesselI(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1647987313994701804789546421624231786109",
                "0.6913178999244037281043977127043561292605+0.7297627072739941639722102363516970518100i",
                "-147.3852070550464048910463035826592748950+355.8193657654598022461566541656199097089i",
                "380.8071038174136380895412142134691371888+24.4211076404065761574462513884216499481i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2307132623613207291700148154212940832279",
                "-0.7637050333300319276834674919028063222593+0.6731752912845033119892483630569655672571i",
                "-355.8193266390081256008269123991012091764-147.3851125954361195983671793130670840604i",
                "-24.4212063793515572751294953418478851212+380.8070862281761437277593636914272239233i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "385.1361062175201558966471062006241983383",
                "-123.1664080476302902164025935100940386317+361.1654423457972612830495062782884318498i",
                "-0.0630657441813617428751181251909262640121+0.1522541749237956241801270263453949664500i",
                "0.4096569220419722778270335730130202170039+0.9179625558304949894373042261707180088511i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}