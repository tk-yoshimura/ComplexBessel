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
        public void BesselJNuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2438916009944455930509024006822094590044",
                "0.8595350248736214788314033339982778177531+0.4095305874830826991499563334758218582224i",
                "-99.2093378490127544059778361466223235514+239.5125289491410127442475837441638935363i",
                "263.1863222110156427169044319492286989283-8.6304727131861880681028896099643316695i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.1578864406955335255015814633461565067443",
                "-0.4383779385535059719943648529327205271505+0.8277695321264137670785830819869276383458i",
                "-239.5124733658510438511085950787753606038-99.2092036590802701877532774570760269600i",
                "8.6303372703820257216095566129115487650+263.1862870585406456034047796762068619319i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "259.2464932068564478111058220152804824509",
                "-108.6905622308658007172640590065555120097+239.8497174068626787608849610725518054656i",
                "-0.0933332749935713814195454669945824550087+0.2253264583101776746266894266920051249616i",
                "0.0494271141562927878056095328083569037162+0.9508273878338364036046455424570177729892i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.1129917204240752499955502449543583481812",
                "-0.2920321596896477435835730063569339248031+0.9255247080088065395168025748585358079659i",
                "-327.5958315261647606226260197816271964032",
                "118.1993372362088692954958992751177392489+305.9649315869144812446885700585963491194i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2630366048203780940912122170454912468567",
                "-0.9598359655038815632385364472100267815535-0.2710759223232147404864880325988818725306i",
                "0.0001179739250805254279489835832574471-327.5958315261647606226260197816271964032i",
                "-305.9649965214768384196294150273788377167+118.1994319862994430039537174469856505201i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "327.5958315261647606226260197816271964032",
                "-118.1993372362088692954958992751177392489+305.9649315869144812446885700585963491194i",
                "0.1129917204240752499955502449543583481812",
                "0.2920321596896477435835730063569339248031+0.9255247080088065395168025748585358079659i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2796692128288667274428623784399225778619",
                "-1.0288565828861527766808513568245201346070-0.0400133455759930485396999121947086674212i",
                "-147.3851793884671005673716123311225297439-355.8192989724288212771522211475644640253i",
                "-286.5396860688984519914924100709499972379+252.0029421792155869311268357655892558287i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.04660861182185516786099534902358450609043",
                "0.0511854329017749640224751891540278299565-0.9920269723649964009426949306983220909872i",
                "355.8193380988804979224819629140831645579-147.3852738480773858600507366007147205786i",
                "-252.0028847977071317412842311959296137343-286.5396038124517804738619566641598059086i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "385.1360339212643882075688494665499538842",
                "-123.1663698383411667769737547067553404568+361.1655020905082446833876126249957812383i",
                "0.1070247742921935169996122049078168973888+0.2583806616061329528564769948405731583932i",
                "0.4306938795552514466252833659793766480125+0.9352270943925872815596475167233643156748i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2346363468539146243812766515904546115488",
                "-0.9210727287952035715695186441612476276637-0.4500097492594776163232664756529305323059i",
                "-399.8731367825600982190830861458227548896i",
                "-375.2248883656372164073488879012350645410+124.0494248516583719404915270407171445874i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(-1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1580604617312474942555552661874835503273",
                "0.4751132311857351322507571130069904404809-0.8924060223921919876891883801254571132552i",
                "399.8731367825600982190830861458227548896-0.0000989111122522303538817640329167230i",
                "-124.0493425847154086872860238678454883630-375.2248366059818424281756969924508622243i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(-1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "399.8731367825600982190830861458227548896",
                "-124.0494248516583719404915270407171445874+375.2248883656372164073488879012350645410i",
                "0.2346363468539146243812766515904546115488i",
                "0.4500097492594776163232664756529305323059+0.9210727287952035715695186441612476276637i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(-1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.04104480174033306261896375195155558184215",
                "-0.2727735624946555689624092912803813034216-0.9709512394713786735307219342622057525196i",
                "297.3075777403793113791647071017302300085-297.3075777403793113791647071017302300085i",
                "-367.5842730035107596086193916940176872302-190.7538044519412771761723355694098982582i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2790928085709920614516248870672736661963",
                "1.0085439501536119966106488207754848437719-0.2718888790453185955251840608598293608961i",
                "297.3075108252664289364888200644860096655+297.3075108252664289364888200644860096655i",
                "190.7538951268747822868705129281805200074-367.5842945891375648641857461164027952741i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "420.4564086367377224027943364288013287513",
                "-125.0380234332093153675803374381876821568+394.8046407634099467588291362801105027019i",
                "-0.02902305764304691593054948933637538297441+0.02902305764304691593054948933637538297441i",
                "0.4936861698633118036315660339348718604974+0.8794462414000787273877855304183676208941i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1716508071375539060908694078519720010684",
                "0.5352263164214296202071813035509292802714-0.8738208098768826045918880302438020137565i",
                "427.5641157218047851773967913180828851256",
                "-125.3153849445027676104733828199197106494-401.5501656053972599252701629201056081047i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2235214893875662205273234004986203592748",
                "0.9005390888950140019847811488628206257495+0.5090016470734796277222239522170476897953i",
                "-0.0000932461470174678394785425750282663+427.5641157218047851773967913180828851256i",
                "401.5502135219852371490405216060640693202-125.3154633823511954934224467692757147100i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "427.5641157218047851773967913180828851256",
                "-125.3153849445027676104733828199197106494+401.5501656053972599252701629201056081047i",
                "0.1716508071375539060908694078519720010684",
                "0.5352263164214296202071813035509292802714+0.8738208098768826045918880302438020137565i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2790928085709920614516248870672736661963",
                "1.0085439501536119966106488207754848437719-0.2718888790453185955251840608598293608961i",
                "297.3075108252664289364888200644860096655+297.3075108252664289364888200644860096655i",
                "190.7538951268747822868705129281805200074-367.5842945891375648641857461164027952741i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.04104480174033306261896375195155558184215",
                "0.2727735624946555689624092912803813034216+0.9709512394713786735307219342622057525196i",
                "-297.3075777403793113791647071017302300085+297.3075777403793113791647071017302300085i",
                "367.5842730035107596086193916940176872302+190.7538044519412771761723355694098982582i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "420.4563140044775563257531029648905588402",
                "-125.0379745796920404228457803623330811536+394.8047201436134020202089940928665959301i",
                "0.1973484175209474790666099872207749384389+0.1973484175209474790666099872207749384389i",
                "0.5208937961761326290976649941806354211097+0.9054027363804401940962223897362469952305i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2346363468539146243812766515904546115488",
                "0.9210727287952035715695186441612476276637+0.4500097492594776163232664756529305323059i",
                "399.8731367825600982190830861458227548896i",
                "375.2248883656372164073488879012350645410-124.0494248516583719404915270407171445874i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.1580604617312474942555552661874835503273",
                "-0.4751132311857351322507571130069904404809+0.8924060223921919876891883801254571132552i",
                "-399.8731367825600982190830861458227548896+0.0000989111122522303538817640329167230i",
                "124.0493425847154086872860238678454883630+375.2248366059818424281756969924508622243i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "399.8731367825600982190830861458227548896",
                "-124.0494248516583719404915270407171445874+375.2248883656372164073488879012350645410i",
                "0.2346363468539146243812766515904546115488i",
                "0.4500097492594776163232664756529305323059+0.9210727287952035715695186441612476276637i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
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

        [TestMethod()]
        public void BesselJNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.1129917204240752499955502449543583481812",
                "-0.2920321596896477435835730063569339248031+0.9255247080088065395168025748585358079659i",
                "-327.5958315261647606226260197816271964032",
                "118.1993372362088692954958992751177392489+305.9649315869144812446885700585963491194i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2630366048203780940912122170454912468567",
                "-0.9598359655038815632385364472100267815535-0.2710759223232147404864880325988818725306i",
                "0.0001179739250805254279489835832574471-327.5958315261647606226260197816271964032i",
                "-305.9649965214768384196294150273788377167+118.1994319862994430039537174469856505201i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "327.5958315261647606226260197816271964032",
                "-118.1993372362088692954958992751177392489+305.9649315869144812446885700585963491194i",
                "0.1129917204240752499955502449543583481812",
                "0.2920321596896477435835730063569339248031+0.9255247080088065395168025748585358079659i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2840999778108356335127490876541630028614",
                "-0.9177630578292492700452390087895244832607+0.2957395939136044642915387946975069906045i",
                "-99.2092985456914967434350511345434417219-239.5124340625297861885658100666964787744i",
                "-179.9982631431427018260236710899959316364+192.2034740747427997948568611516414999797i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselJ(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.06081483206439675962631471644838838250142",
                "-0.2978030316817218234246120347257738098577-0.8749033049388012975197579062382920970133i",
                "239.5124896458197550817047987320850117069-99.2094327356239809616596098240897383133i",
                "-192.2034031587709596276249574615153547044-179.9981425140639697597553313248743302966i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselY(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "259.2463905023285440105301620419008985411",
                "-108.6905026158461204681000083585844502836+239.8497963769473780615563866751557904423i",
                "0.1087203546434964230675715646752914150022+0.2624741546863417656896731171066431523443i",
                "0.0779849592979988073945703062980994529966+0.9610771477085318169711289133567321192352i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = MillerBackward<Pow2.N4>.BesselI(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}