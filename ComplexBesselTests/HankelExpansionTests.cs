using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class HankelExpansionTests {
        readonly Complex<Pow2.N4>[] zs = [
            (64,0),(64,1),(64,2),(64,4),(64,8),(64,16),(64,32),(64,64),
            (0,64),(1,64),(2,64),(4,64),(8,64),(16,64),(32,64),(64,64),
        ];

        readonly Complex<Pow2.N4>[] zs_mini = [
            (64, 8), (8, 64)
        ];

        [TestMethod()]
        public void ACoefTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Limit<Pow2.N4>.HankelExpansion hankel = new(nu);

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
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {Limit<Pow2.N4>.BesselJ(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {Limit<Pow2.N4>.BesselY(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -4; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {Limit<Pow2.N4>.BesselI(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 4; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs) {
                    Console.WriteLine($"{z}: {Limit<Pow2.N4>.BesselK(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "9.2145929015511005131094720273851178236+147.6159570273588415961403709192966375974i",
                "2.922066696281459872639446880884572748286e26+9.28215852645366138489133371390795409174e25i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselJ(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-147.6159903358803008197446577402684011815+9.2145949710669354056471653366179386559i",
                "-9.28215852645366138489133371390795409174e25+2.922066696281459872639446880884572748286e26i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselY(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-2.60666884921141062870887160247972666843e25+3.054850441793332396273499925212507250408e26i",
                "132.8530893302299147044482156602733629886+65.0033548892541275819643373472106188313i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselI(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-5.25616659269784449308681704247624052439e-30-2.472840504427643132508872150883116217499e-29i",
                "-0.00001701900642617936969720732433114038442166-0.00004958223421184402515747294308650095460414i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselK(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM1p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-147.2999139438063945323394630815032224390+10.9056339300812418071382865800802817172i",
                "-9.13255376986751184684287116150543446463e25+2.889631506025360495813279339858706432795e26i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselJ(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM1p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-10.9056322416899217958437967766265351684-147.2998805626433483749635926306468094384i",
                "-2.889631506025360495813279339858706432795e26-9.13255376986751184684287116150543446463e25i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselY(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM1p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-2.62076152242149906566470469518604540524e25+3.019159107207759513123694167620311528728e26i",
                "131.9139702088759482017693824281792414315+66.4447286321608627589382194728601568628i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselI(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-5.35222404565129229070775331686841841211e-30-2.500563373120195349488734686337018462950e-29i",
                "-0.00001761577059679953277111031956241000418243-0.00004945855291014150741401368246876143340060i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselK(1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu1p375Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-19.3844491315920837390199087861121994352+146.5832933973515097634608469293879648659i",
                "2.678487439866468568876646704464127999900e26+1.475828794352374901455927115895960295955e26i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselJ(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p375Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-146.5833256942998126925823267036391419176-19.3844406821480655811156969072786874160i",
                "-1.475828794352374901455927115895960295955e26+2.678487439866468568876646704464127999900e26i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselY(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p375Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-2.60981095641589821397947573363886292337e25+3.047007460358397106095076131919830300429e26i",
                "132.6489770521213388555665512843130872017+65.3197340654040189559932767827654910729i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselI(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p375Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-5.27705248199702819773776225908853739551e-30-2.478879067656174354854949126828274145020e-29i",
                "-0.00001714958853999699091968918500505305032212-0.00004955578201524993327661980544066492052557i",
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = Limit<Pow2.N4>.BesselK(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}