using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class HankelExpansionTests {
        readonly Complex<Pow2.N4>[] zs = [
            (64,-64),(64,-32),(64,-16),(64,-8),(64,-4),(64,-2),(64,-1),(64,-0.0),
            (64,0.0),(64,1),(64,2),(64,4),(64,8),(64,16),(64,32),(64,64),

            (-64,-64),(-64,-32),(-64,-16),(-64,-8),(-64,-4),(-64,-2),(-64,-1),(-64,-0.0),
            (-64,0.0),(-64,1),(-64,2),(-64,4),(-64,8),(-64,16),(-64,32),(-64,64),

            (-64,64),(-32,64),(-16,64),(-8,64),(-4,64),(-2,64),(-1,64),(-0.0,64),
            (0.0,64),(1,64),(2,64),(4,64),(8,64),(16,64),(32,64),(64,64),

            (-64,-64),(-32,-64),(-16,-64),(-8,-64),(-4,-64),(-2,-64),(-1,-64),(-0.0,-64),
            (0.0,-64),(1,-64),(2,-64),(4,-64),(8,-64),(16,-64),(32,-64),(64,-64),
        ];

        readonly Complex<Pow2.N4>[] zs_mini = [
            (64, -8), (64, 8), (-64, -8), (-64, 8), (-8, 64), (8, 64), (-8, -64), (8, -64)
        ];

        [TestMethod()]
        public void ACoefTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                HankelExpansion<Pow2.N4> hankel = new(nu);

                Console.WriteLine(nu);

                for (int k = 0; k <= 16; k++) {
                    Console.WriteLine($"a_{k}: {hankel.ACoef(k)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                HankelExpansion<Pow2.N4> hankel = new(nu);

                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {hankel.BesselJ(z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                HankelExpansion<Pow2.N4> hankel = new(nu);

                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {hankel.BesselY(z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                HankelExpansion<Pow2.N4> hankel = new(nu);

                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {hankel.BesselI(z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                HankelExpansion<Pow2.N4> hankel = new(nu);

                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {hankel.BesselK(z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJMiniTest() {
            HankelExpansion<Pow2.N4> hankel = new(1.25);

            Complex<Pow2.N4>[] expecteds = [
                "9.2145929015511005131094720273851178236-147.6159570273588415961403709192966375974i",
                "9.2145929015511005131094720273851178236+147.6159570273588415961403709192966375974i",
                "-110.8959453519476388434291742350541699970-97.8645430988272226734726977555188560612i",
                "-110.8959453519476388434291742350541699970+97.8645430988272226734726977555188560612i",
                "-2.722560899830383623022272861089442378455e26-1.409865452209600514674146268526138161281e26i",
                "2.922066696281459872639446880884572748286e26+9.28215852645366138489133371390795409174e25i",
                "-2.722560899830383623022272861089442378455e26+1.409865452209600514674146268526138161281e26i",
                "2.922066696281459872639446880884572748286e26-9.28215852645366138489133371390795409174e25i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = hankel.BesselJ(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYMiniTest() {
            HankelExpansion<Pow2.N4> hankel = new(1.25);

            Complex<Pow2.N4>[] expecteds = [
                "-147.6159903358803008197446577402684011815-9.2145949710669354056471653366179386559i",
                "-147.6159903358803008197446577402684011815+9.2145949710669354056471653366179386559i",
                "-97.8645180827771469333720266685710849958+110.8959674412603533326241154327892910921i",
                "-97.8645180827771469333720266685710849958-110.8959674412603533326241154327892910921i",
                "1.409865452209600514674146268526138161281e26-2.722560899830383623022272861089442378455e26i",
                "-9.28215852645366138489133371390795409174e25+2.922066696281459872639446880884572748286e26i",
                "1.409865452209600514674146268526138161281e26+2.722560899830383623022272861089442378455e26i",
                "-9.28215852645366138489133371390795409174e25-2.922066696281459872639446880884572748286e26i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = hankel.BesselY(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselIMiniTest() {
            HankelExpansion<Pow2.N4> hankel = new(1.25);

            Complex<Pow2.N4>[] expecteds = [
                "-2.60666884921141062870887160247972666843e25-3.054850441793332396273499925212507250408e26i",
                "-2.60666884921141062870887160247972666843e25+3.054850441793332396273499925212507250408e26i",
                "-1.975786140944273667434472492382389247086e26-2.344424784861298183212213072859469902627e26i",
                "-1.975786140944273667434472492382389247086e26+2.344424784861298183212213072859469902627e26i",
                "-139.9056334090550484779115612126481548920-47.9770073249204233223395686709723050197i",
                "132.8530893302299147044482156602733629886+65.0033548892541275819643373472106188313i",
                "-139.9056334090550484779115612126481548920+47.9770073249204233223395686709723050197i",
                "132.8530893302299147044482156602733629886-65.0033548892541275819643373472106188313i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = hankel.BesselI(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKMiniTest() {
            HankelExpansion<Pow2.N4> hankel = new(1.25);

            Complex<Pow2.N4>[] expecteds = [
                "-5.25616659269784449308681704247624052439e-30+2.472840504427643132508872150883116217499e-29i",
                "-5.25616659269784449308681704247624052439e-30-2.472840504427643132508872150883116217499e-29i",
                "-9.597095705753467335604038328256260292002e26-8.18909170702392813709399161901292238171e25i",
                "-9.597095705753467335604038328256260292002e26+8.18909170702392813709399161901292238171e25i",
                "-204.2140852044501196181755107589505155153-417.3703365407477352155598375752493218656i",
                "-0.00001701900642617936969720732433114038442166-0.00004958223421184402515747294308650095460414i",
                "-204.2140852044501196181755107589505155153+417.3703365407477352155598375752493218656i",
                "-0.00001701900642617936969720732433114038442166+0.00004958223421184402515747294308650095460414i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = hankel.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}