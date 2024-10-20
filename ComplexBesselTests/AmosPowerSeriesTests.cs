using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBesselTests {
    [TestClass()]
    public class AmosPowerSeriesTests {
        readonly Complex<Pow2.N4>[] zs_mini = [
            (0.5, 0), (1, 0), (2, 0), (1, 0.5), (0, 0.5), (0, 1), (0, 2), (0.5, 1), (1, 1)
        ];

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = double.ScaleB(1, -16); nu <= 0.5; nu *= 2) {
                Console.WriteLine(nu);

                foreach (Complex<Pow2.N4> z in zs_mini) {
                    Console.WriteLine($"{z}: {AmosPowerSeries<Pow2.N4>.BesselK(nu, z)}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.9244190712276658617819241675302169895388",
                "0.4210244382407083333356273792126090361362",
                "0.1138938727495334356527195749324818329983",
                "0.3078189429740561651197317824086681489296-0.2601600399489966290917754536027510089646i",
                "0.698248393783854194778283242924061814301-1.474144926021783584241760245878163763477i",
                "-0.1386337152040539996810990857549921961745-1.2019697153172064991366624462957556118924i",
                "-0.8016962318836942154259743686714061965573-0.3516868134783004458924008931401674893990i",
                "0.0583659790931038640803753116433600481447-0.6764549973133448353518414219607300433577i",
                "0.0801977269465178187269687365642791668341-0.3572774592853302506059456932500239816608i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNuEmpM16Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.9244190713592872307278672497149446101236",
                "0.4210244382765422549804078760543075399232",
                "0.1138938727550307535000208267324656385440",
                "0.3078189429893138478777387812218825937231-0.2601600399768814560126660771504269865347i",
                "0.698248393618497894307920796172138094816-1.474144926234965897798253760219646380486i",
                "-0.1386337153201636711587706602792928565800-1.2019697153427340640931717951727782375048i",
                "-0.8016962319114101644671424016520487836021-0.3516868134396126483566443287350537285567i",
                "0.0583659790525546228635048135259772095935-0.6764549973529138102384061321482056170998i",
                "0.0801977269364759333169807107950869586719-0.3572774593075529610098245586561226352440i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(double.ScaleB(1, -16), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNuEmpM12Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.9244191049227368194856062156874016012845",
                "0.4210244474141923610715662926755202063049",
                "0.1138938741568468124412096192582511725723",
                "0.3078189468800229628567362320232196348486-0.2601600470875123917423856433952335802843i",
                "0.698248351452640074331513202258879834122-1.474144980596455768556735275291613181210i",
                "-0.1386337449281301053035586886436081113812-1.2019697218522628640412810076416602956350i",
                "-0.8016962389789771227369868335623161873275-0.3516868035742242135409957827670405354598i",
                "0.0583659687124979763639023170313392323692-0.6764550074430024018166872161093204597663i",
                "0.0801977243757951133936497769150671959632-0.3572774649743441388871341221913889966829i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(double.ScaleB(1, -12), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNuEmpM8Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.9244276971990933975611621662079432779763",
                "0.4210267866582995867945764824954606041765",
                "0.1138942330222742763285847715820100974010",
                "0.3078199429023223179910629893544170799532-0.2601618674136785474727241034348215958645i",
                "0.698237556914393820950298669756362164962-1.474158897132216846568753249569255467338i",
                "-0.1386413245811251636548943593860766293547-1.2019713882743368812971653289309773678323i",
                "-0.8016980482730233474039985058995355057436-0.3516842780306432392251952750166702371274i",
                "0.0583633216490707519394969321306519397991-0.6764575905053872096789019636612270751243i",
                "0.0801970688388591429873302134287768054156-0.3572789156745167235709901162998508110720i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(double.ScaleB(1, -8), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNuEmpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.9266295010464901548631112787035406549135",
                "0.4216260055367596629426245602700832855719",
                "0.1139861364205282075443885899447483979437",
                "0.3080749747899242825789422388745447367529-0.2606281755387502471906200077399876012798i",
                "0.695468999204232780158599094861129323283-1.477721157451174529201299256666978859677i",
                "-0.1405826060154866072328721726987680536852-1.2023968580534103282801554507580350944070i",
                "-0.8021610246220980929889348602083177331681-0.3510374674690429240970540576876128372207i",
                "0.0576850883023240000099045377609235978662-0.6771188350127757295871097542568355976287i",
                "0.0800290779085800668629703072579014929372-0.3576504018107018752237188503907149267022i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void Gamma12Test() {
            for (double x = 1 / 32d; x <= 0.5; x += 1 / 32d) {
                MultiPrecision<Pow2.N4> expected_g1 = ((1 / MultiPrecision<Pow2.N8>.Gamma(1 - x) - 1 / MultiPrecision<Pow2.N8>.Gamma(1 + x)) / (2 * x)).Convert<Pow2.N4>();
                MultiPrecision<Pow2.N4> expected_g2 = ((1 / MultiPrecision<Pow2.N8>.Gamma(1 - x) + 1 / MultiPrecision<Pow2.N8>.Gamma(1 + x)) / 2).Convert<Pow2.N4>();

                (MultiPrecision<Pow2.N4> actual_g1, MultiPrecision<Pow2.N4> actual_g2) = AmosPowerSeries<Pow2.N4>.Gamma12(x);

                Console.WriteLine(x);
                Console.WriteLine(expected_g1);
                Console.WriteLine(actual_g1);
                Console.WriteLine(expected_g2);
                Console.WriteLine(actual_g2);

                Assert.IsTrue(MultiPrecision<Pow2.N4>.NearlyEqualBits(expected_g1, actual_g1, 1));
                Assert.IsTrue(MultiPrecision<Pow2.N4>.NearlyEqualBits(expected_g2, actual_g2, 1));
            }
        }

        [TestMethod()]
        public void SinhcTest() {
            for (double x = 0.5d; x >= double.ScaleB(1, -64); x /= 2) {
                Complex<Pow2.N4> expected = Complex<Pow2.N4>.Sinh(x) / x;

                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.Sinhc(x);

                Console.WriteLine(x);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }

            for (double x = 0.5d; x >= double.ScaleB(1, -64); x /= 2) {
                Complex<Pow2.N4> expected = Complex<Pow2.N4>.Sinh((0, x)) / (0, x);

                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.Sinhc((0, x));

                Console.WriteLine((0, x));
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }

            for (double x = 0.5d; x >= double.ScaleB(1, -64); x /= 2) {
                Complex<Pow2.N4> expected = Complex<Pow2.N4>.Sinh((x, x)) / (x, x);

                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.Sinhc((x, x));

                Console.WriteLine((x, x));
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void G1CoefPlot() {
            Console.WriteLine($"new ddouble[]{{");
            for (int i = 0; i < AmosPowerSeries<Pow2.N4>.G1Coef.Count; i++) {
                Console.WriteLine($"    {ToFP128(AmosPowerSeries<Pow2.N4>.G1Coef[i])},");
            }
            Console.WriteLine("},");
        }

        public static string ToFP128<N>(MultiPrecision<N> x) where N : struct, IConstant {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}