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
        public void BesselJNuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5565691569806435311057880380939130630575",
                "0.3290298221768464123213215078946398671009-0.4523808740528944171447563711333104663790i",
                "-0.1979043570785148036862449091704256882579+0.4777833829116782736393753645645450514344i",
                "0.1896628776686639242904197490555560680669+0.0685530092007337641378144060412755662610i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.7995200433890299558413376618184818312221",
                "0.5230096510842060982284813045982633139195-0.1523418017528301804148943093992476988243i",
                "-0.3565348133837406668199006866347372867088+0.0948155838941693455523153737777185081147i",
                "-0.3096801318157763734666315501614217847293+0.2523731707250513061278947550254016052414i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5171490070929148318526967220221328202034",
                "-0.0092461189258685963501503177179227210710+0.2014597516151629711787962879669869787286i",
                "-0.2129897953418971305253262788923010675083+0.5142028525614779174104736905985266529528i",
                "-0.5438596921381461056604579132288904866683+0.1308652526761356599726468629745794593138i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(-2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.3528340286156377191506207876191884610951",
                "0.4126719082931705311258787227976070133621+0.2659739227983885388560494647939763326385i",
                "-0.6889484476987382040549500158118671053314",
                "-0.2611297504430097979157677148626194170926+0.7623197829675602274187133336866540367098i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.6174081041906826664849773650046345099556",
                "-0.5737407339596305445767047913047697553493+0.4104130982552780846484455191480970642099i",
                "0.1615484771878957938043278315661031569039-0.6889484476987382040549500158118671053314i",
                "-0.7512454872676915796452160761587209184193-0.1239501069691137006484431890901290084408i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.6889484476987382040549500158118671053314",
                "0.2611297504430097979157677148626194170926+0.7623197829675602274187133336866540367098i",
                "-0.3528340286156377191506207876191884610951",
                "-0.4126719082931705311258787227976070133621+0.2659739227983885388560494647939763326385i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(-2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.5707185580803000794035039721588725359653",
                "-0.6823739418566101337626293708897207126341+0.3612671129063680504785518444927812148370i",
                "-0.4858684983789953056194193905464810147296-1.1729903183164206138746515755640795098954i",
                "-1.450186349651061018773712838579265363548-0.239757606991499747024910985096081570513i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2016873056052210453059326134908999881165",
                "-0.3651207476585298877502164566719517136282-0.4536895937769066570884411363069767003523i",
                "1.2111777577691119747498044865454159567909-0.5780611326179841910306631536919659191842i",
                "0.329480690088410684155416805197902026103-1.436812883148947493913490113969087775686i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM1p25Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.269635571564186098912998845758437679230",
                "0.776469435716799839858808256038737480114+1.248046222790724323797250073554682378759i",
                "0.2184045367206240740667380001971245370696+0.5272751946347432960019028290857091703068i",
                "-0.0726340891574538000094265768000748961408+0.7686822571680042886816102564735492179252i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(-1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.2347857104062484691740346837934112099540",
                "-0.4644716341321899286815541648754870514760-0.4762103386692476645735662030782459997737i",
                "1.500898928317347896363305470899432470243-1.500898928317347896363305470899432470243i",
                "0.066515696518825564688658993687305865069-1.955436413661095548394436703184904215721i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5130161365618277516656918486272844223548",
                "0.6686912967460061875290077268565560774343-0.4259785835635051158727973500392477374496i",
                "1.446907961804159844412741270488718888503+1.446907961804159844412741270488718888503i",
                "1.986613473041921807184818958223952427961+0.001380872835431496519696154959868778199i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "2.122591620177637193816120295731577941119",
                "1.429736048342671509717252365702859667107+1.335668648215055530783472968393760270490i",
                "-0.1660185679539592503648086084118448496451+0.1660185679539592503648086084118448496451i",
                "0.0083005175804988171088098763545968877486+0.6651626019078359820136260671446847456382i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(-0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2238907791412356680518274546499486258252",
                "0.1878537280824617161908894244937405040159-0.6461694351539807163814284884027101536900i",
                "2.279585302336067267437204440811533353286",
                "1.586259450202371272377863210976327103297-1.391602452327335923396791402078602041710i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5103756726497451195966065927271578732681",
                "0.8004511204099939791665822433811426234450+0.0756385502863937918135993035331327647987i",
                "-0.072507091343870251624958589644410047299+2.279585302336067267437204440811533353286i",
                "1.367418716811797870340976513222951111228+1.521506576945447807137181115295457565707i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "2.279585302336067267437204440811533353286",
                "1.586259450202371272377863210976327103297+1.391602452327335923396791402078602041710i",
                "0.2238907791412356680518274546499486258252",
                "0.1878537280824617161908894244937405040159+0.6461694351539807163814284884027101536900i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu0Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1138938727495334356527195749324818329983",
                "0.03798772291598645925495274330008729825901-0.10171357546139087330069486667139851472545i",
                "-0.8016962318836942154259743686714061965573-0.3516868134783004458924008931401674893990i",
                "-0.2423451044918719926095078380947558624556-0.1762671890926997383148292291783157322613i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(0, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5130161365618277516656918486272844223548",
                "0.6686912967460061875290077268565560774343-0.4259785835635051158727973500392477374496i",
                "1.446907961804159844412741270488718888503+1.446907961804159844412741270488718888503i",
                "1.986613473041921807184818958223952427961+0.001380872835431496519696154959868778199i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2347857104062484691740346837934112099540",
                "0.4644716341321899286815541648754870514760+0.4762103386692476645735662030782459997737i",
                "-1.500898928317347896363305470899432470243+1.500898928317347896363305470899432470243i",
                "-0.066515696518825564688658993687305865069+1.955436413661095548394436703184904215721i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "2.046236863089055036605183612020732319268",
                "1.405724282930391339128511248769814108893+1.403771433838611525810476905829853342329i",
                "0.3627571890209923192100892172181200144976+0.3627571890209923192100892172181200144976i",
                "0.1716238053715320540926810159452733768212+0.7740484955275217776394726747009263284040i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu0p5Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1199377719680614473680365016367935162195",
                "0.03771759290947043839144489431972342692531-0.10697560550198205009609217449211541812340i",
                "-0.8305990167544448610670548592703028877569-0.3090363033107886090210885121743461951240i",
                "-0.2565476205422207937049886251375604711910-0.1710375617375823177637363016820303105362i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(0.5, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.5767248077568733872024482422691370869203",
                "0.7906233925534283360794576590957144279037-0.0799326941677760538669949499630600623626i",
                "1.590636854637329063382254424999666247954i",
                "1.291847519239456433209125816448801848102+1.010488365079473626473056461917733683704i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.1070324315409375468883707722774766366875",
                "-0.0163154378204725036411449592642420464086+0.5994068417668535937569835487194162630565i",
                "-1.590636854637329063382254424999666247954+0.089041385844025542179369241921693109605i",
                "-1.089469855204280961140857707673213460862+1.314951464532273343940129400166894471442i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.590636854637329063382254424999666247954",
                "1.010488365079473626473056461917733683704+1.291847519239456433209125816448801848102i",
                "0.5767248077568733872024482422691370869203i",
                "-0.0799326941677760538669949499630600623626+0.7906233925534283360794576590957144279037i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1398658818165224272845988070354110238872",
                "0.03629159240042704557109597463199955990743-0.12406383457283476223609599758866618265197i",
                "-0.9059172095959896177775537048007528118840-0.1681261503124309352277046804553284122626i",
                "-0.3003622555979414051258616290696074264864-0.1511862121880112240333746913117548598634i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
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

        [TestMethod()]
        public void BesselJNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.3528340286156377191506207876191884610951",
                "0.4126719082931705311258787227976070133621+0.2659739227983885388560494647939763326385i",
                "-0.6889484476987382040549500158118671053314",
                "-0.2611297504430097979157677148626194170926+0.7623197829675602274187133336866540367098i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.6174081041906826664849773650046345099556",
                "-0.5737407339596305445767047913047697553493+0.4104130982552780846484455191480970642099i",
                "0.1615484771878957938043278315661031569039-0.6889484476987382040549500158118671053314i",
                "-0.7512454872676915796452160761587209184193-0.1239501069691137006484431890901290084408i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.6889484476987382040549500158118671053314",
                "0.2611297504430097979157677148626194170926+0.7623197829675602274187133336866540367098i",
                "-0.3528340286156377191506207876191884610951",
                "-0.4126719082931705311258787227976070133621+0.2659739227983885388560494647939763326385i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2537597545660558629373183819678928568856",
                "0.0173954630071941908173911239702204731242-0.2154812800798295013180100545951312848100i",
                "-0.9698223821961251506536790491267346088198+0.5542303961176891718851528116605853224851i",
                "-0.4834389764814575338865522427720027209409+0.0035481305104488961725101975526682649824i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.1717922192746526626093285041524474842609",
                "0.1371644524338003115773501242626902895670+0.2121596626442962346883776679923519561052i",
                "-0.1121686713561415324767021066612574724344-0.2707991276613673956103578841255690708853i",
                "-0.3530888281445462507798597484218917583485+0.1299804827326514264057628281377628851272i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselJ(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "-0.9588998694752389125719107550253147367972",
                "-0.6024828893815042678354883775511311202788+0.4276035047995020196735373521301689504990i",
                "0.3920476971893050024298325620553768356110-0.4048886123288256817152623896094016688070i",
                "0.0848650142668125060285480272777539359506-0.2269290780858166843140195780534000126947i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselY(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselINu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.2931108636266483420273007767328058373126",
                "0.0150349370614888906750899530276686514847+0.3759529187537437911766495506975738271829i",
                "-0.0657420361256402136549843911904251195158-0.1587153152325425679150590247230866901759i",
                "-0.2485005333974095081333174739613941018161-0.0455334422815672596624936499071866197056i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselI(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu2p75Test() {
            Complex<Pow2.N4>[] expecteds = [
                "0.4976876225514758230877801935227240069764",
                "-0.0539389446833493112606483330612126576354-0.3876263575526084173560816776177924676567i",
                "-0.327102278551979591598151520305248252672+1.494848323234421118728815288349863674038i",
                "-0.6561230834990135580985217731908020542535+0.3918593760704902663093507954799350779213i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = PowerSeries<Pow2.N4>.BesselK(2.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }
    }
}