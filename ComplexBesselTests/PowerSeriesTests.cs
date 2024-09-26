using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class PowerSeriesTests {
        readonly Complex<Pow2.N4>[] zs = [
            (64,0.0),(64,1),(64,2),(64,4),(64,8),(64,16),(64,32),(64,64),
            (0.0,64),(1,64),(2,64),(4,64),(8,64),(16,64),(32,64),(64,64),
        ];

        readonly Complex<Pow2.N4>[] zs_mini = [
            (2, 0), (2, 1), (0, 2), (1, 2)
        ];

        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {PowerSeries<Pow2.N4>.BesselJ(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {PowerSeries<Pow2.N4>.BesselY(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {PowerSeries<Pow2.N4>.BesselI(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {PowerSeries<Pow2.N4>.BesselK(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5461734240402840405040192805817497539907",
                "0.7406905982130526525282996354298565947518+0.0653525629576418168615275118992889618600i",
                "-0.5128710957721440678874648166285552435091+1.2381803551622607370178499126274901855705i",
                "0.792458571590892412313944059884863669194+1.185514362715563960859494894140685837034i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2609445010948932850970918858558465760891",
                "-0.2243318849705551971234375715163375427610+0.5762614136691997211320426986448840436622i",
                "-1.1999929157095693761426970016461537386749-0.4206784615331551824762210534830703390545i",
                "-1.258414632053969395019982684095835757057+0.846445903226066623816226923398194765228i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.340196758982897224249773547720237183430",
                "0.792011689027532346837707971706087636462+1.185812960098059774207499817142830917950i",
                "-0.2090115205783295335819237869639756567965+0.5045984476724264067489018689079643633987i",
                "-0.2230721251310020601843197486720117690441+0.7093182267190664075595593323427693430941i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1567475478393932155730106536544565311307",
                "0.0345262060273944179327551655365097380708-0.1382475505015834566909735194340404198340i",
                "-0.9494797542342198127922728913694070705297-0.0503750660859061337611722183200748330109i",
                "-0.3341892916384660458080943517369267249452-0.1318737190109499254844835616673640204460i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}