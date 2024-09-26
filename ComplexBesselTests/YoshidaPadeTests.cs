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
        public void BesselKNu0Test() {
            YoshidaPade<Pow2.N4> pade = new(0);

            Complex<Pow2.N4>[] expecteds = [
                "0.0001464707052228153870965844086986779219673",
                "-0.0000752672003871676039474869759917141832085-0.0001232098841922133928827997847868472206244i",
                "-0.3511067344897134855257515184958057642762-0.2696284573430488985964559109599364607694i",
                "-0.04196897453996195607319429422296981610021-0.04119361428317056483224333622558133825219i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = pade.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu0p5Test() {
            YoshidaPade<Pow2.N4> pade = new(0.5);

            Complex<Pow2.N4>[] expecteds = [
                "0.0001486480066651728298787091323662822494168",
                "-0.0000767389254864942166147125308804301885050-0.0001246901320077561362108244312234056496959i",
                "-0.3555834816787445634443971339041685998763-0.2644048570026352726652935856384220247056i",
                "-0.04273763947298712951180590896644054718980-0.04077236697162196416456803343181723693125i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = pade.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1Test() {
            YoshidaPade<Pow2.N4> pade = new(1);

            Complex<Pow2.N4>[] expecteds = [
                "0.0001553692118050011339168624506224746211171",
                "-0.0000813038765376162163181515712117942236129-0.0001292246118233234030480928708481418493529i",
                "-0.3685659117707023895587438860176189398935-0.2482807926989488900376426504422463573417i",
                "-0.04502955711915833929891439287722383386132-0.03943245719952738186290546219448643638179i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = pade.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
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

        [TestMethod()]
        public void BesselKNu2Test() {
            YoshidaPade<Pow2.N4> pade = new(2);

            Complex<Pow2.N4>[] expecteds = [
                "0.0001853130081740656705758000213542965772466",
                "-0.0001019989720326845609075281027973211681382-0.0001488330942366650044048127207974809367302i",
                "-0.4131769326644507080351621811063673536116-0.1774869794003733012067699394555317257961i",
                "-0.05389599724097771294087289667327390900541-0.03291798067863491922443438861885140301316i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = pade.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu2p75Test() {
            YoshidaPade<Pow2.N4> pade = new(2.75);

            Complex<Pow2.N4>[] expecteds = [
                "0.0002281520973477027067949682707390443338159",
                "-0.0001324312768976771167409628920725297231692-0.0001754274209677220047745542198008392981400i",
                "-0.4488503118537490779632456563902180546862-0.0825214332107734056661935355351039671171i",
                "-0.06343958146532723270082014302166117683223-0.02276924163375224299848779442436235999847i"
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