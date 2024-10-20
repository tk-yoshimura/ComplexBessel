using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class AmosPowerSeriesTests {
        readonly Complex<Pow2.N4>[] zs_mini = [
            (0.5, 0), (1, 0), (2, 0), (1, 0.5), (0, 0.5), (0, 1), (0, 2), (0.5, 1), (1, 1)
        ];

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 0.5; nu = nu < double.ScaleB(1, -16) ? double.ScaleB(1, -16) : nu * 2) {
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
        public void BesselKNuExpM16Test() {
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
        public void BesselKNuExpM12Test() {
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
        public void BesselKNuExpM8Test() {
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
        public void BesselKNuExpM4Test() {
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
        public void BesselKNu1Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.656441120003300893696445403174091511534",
                "0.6019072301972345747375400015356173392616",
                "0.1398658818165224272845988070354110238872",
                "0.3763244754275179194561865147382986231027-0.4018549385212971723099386790277225506824i",
                "-0.380554403413956781197401076552306136589-2.311383429386515572834192626985747679135i",
                "-0.6912298436920842628831416638469705187187-1.2271262301435714892432807817072674867337i",
                "-0.9059172095959896177775537048007528118840-0.1681261503124309352277046804553284122626i",
                "-0.1568703721986486110345570611228852554130-0.8422945505079033891740771412106772613429i",
                "0.0245683055237403486124773461852012359952-0.4597194738011893647760430085185719981090i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1PlusExpM16Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.656469331455259481715793954532103579293",
                "0.6019136545995693803408905332664963917406",
                "0.1398667507669075088571957344438708120123",
                "0.3763266451106571033096328811711664698928-0.4018599931435846322139776333006270816745i",
                "-0.380599391480143012869523553122251987247-2.311404738206777656626698438750831335017i",
                "-0.6912481844734710144742849417745625044042-1.2271241146474678209069595860748549612352i",
                "-0.9059198927573463370738678247000992771522-0.1681200337937813490662888434374414974309i",
                "-0.1568782735716401849841627014574382519441-0.8422993917639777759698735828064896309310i",
                "0.0245661915432033688595253494090968953375-0.4597228115038011712346730833089195296245i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 + double.ScaleB(1, -16), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1PlusExpM12Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.656892604225782814459136688172978519933",
                "0.6020100396546593001166488865482326340895",
                "0.1398797872056888185609511425140495739339",
                "0.3763591948102353541998357259875407908110-0.4019358275567665759544514336910885190088i",
                "-0.381274388475230335965807611589731240301-2.311724363376317952517210175035243317263i",
                "-0.6915233391656589889437974136132893831539-1.2270923549623681434515177340014455739151i",
                "-0.9059601411138248194533639241633647947167-0.1680282711806665885483668284262008757820i",
                "-0.1569968189770427062189923767873378660424-0.8423720171076278053085933688504441288329i",
                "0.0245344746480383994749262871174829762889-0.4597728846797539333233020625226682909109i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 + double.ScaleB(1, -12), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1PlusExpM8Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.663690784376086234119981647495951109495",
                "0.6035570592365777825039112612898706928421",
                "0.1400889275506483188661172503071105516325",
                "0.3768811256389344611169632723209984701746-0.4031530300842314385056382263770696772121i",
                "-0.392119323043925907494493490044783243367-2.316836467476535006848738100068994627064i",
                "-0.6959367780328729602254071666655578214180-1.2265772250871609365409584989171654993276i",
                "-0.9066043488973520126922130846670843127321-0.1665562835467506185285946952317698380834i",
                "-0.1588998836112527351579239868522528506067-0.8435356767673385629820031660172276014753i",
                "0.0240251677485420064615610101159373277842-0.4605760043298218808554790724418483849029i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 + double.ScaleB(1, -8), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1PlusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.779338224468948852146131065283833018615",
                "0.6295893924357833257752022581164982439835",
                "0.1435806094626021159441643194370303510366",
                "0.3855236897644546089905058953340739773531-0.4236430226960778824193519366325911384766i",
                "-0.577579995023674385484887359210048852806-2.397801247429434323148267595910543160111i",
                "-0.769382529478720517868813430650828633944-1.216432424936543777205460175928132386446i",
                "-0.9169521096524708049986980096988312713945-0.1420216478562205349301200545194847918920i",
                "-0.1910195312870529366068036951608402863521-0.8625605105415607927191126652200254836318i",
                "0.0153898264526934755519609711427045869056-0.4739308488364714625254706585285193540195i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 + double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1MinusExpM16Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.656412909392798587719777849353287969881",
                "0.6019008059533817417222693049987945390935",
                "0.1398650128843276286271923771924863444580",
                "0.3763223057814822773057297706794026351980-0.4018498840246522216089983477667461888942i",
                "-0.380509416814235117001738015173019831096-2.311362120506971840104829058819424344166i",
                "-0.6912115032687791585763328796470242726465-1.2271283454127058630111819606677589923893i",
                "-0.9059145264424439992343848436770167720599-0.1681322667074762863551753659960604556948i",
                "-0.1568624710323866049915160546379708394194-0.8422897093060254075243661924773978321025i",
                "0.0245704194443931130768668539844121577227-0.4597161361622141956961243954843734601422i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 - double.ScaleB(1, -16), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1MinusExpM12Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.655989851193631704283857775852971291806",
                "0.6018044613111953390028444108193871652388",
                "0.1398519810840684946321324878218478524563",
                "0.3762897655433071249939819980928823521998-0.4017740816503105119414941080038202208065i",
                "-0.379834793767616609629326424797454274303-2.311042480220602483694620806673662379925i",
                "-0.6909364398874108137172147568712768537127-1.2271600472206338319794979872673848356882i",
                "-0.9058742800777959561199798919898376852110-0.1682239978015110653849725381023635506324i",
                "-0.1567439783430244431667913080148403461817-0.8422170977824582743799713977859574102930i",
                "0.0246021210690829671081034073189664532910-0.4596660792136039969824131992660505037545i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 - double.ScaleB(1, -12), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1MinusExpM8Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.649246601607350399783088488648951995777",
                "0.6002677874658098169820786325203749068720",
                "0.1396440282029050599742087160431183246896",
                "0.3757702568314626920903713914607671809120-0.4005650810913261903402116288888694888631i",
                "-0.369085590276235019266856720520263626689-2.305926505649747603396389683354051468237i",
                "-0.6865463765619797033402730625332845607319-1.2276603604354509768820123610998692229945i",
                "-0.9052305821862948125389896656240504446183-0.1696879165478171159793915643954398986527i",
                "-0.1548544090507474111972417774922784686991-0.8410569760367025357521347415549696922418i",
                "0.0251075187127437938857159316115274606188-0.4588671137653165063767356336538459332760i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 - double.ScaleB(1, -8), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu1MinusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "1.547680849207326313430353245607947854887",
                "0.5768861417436883679073741880827378332871",
                "0.1364564759363191029726400325654835761652",
                "0.3677476010623996930921426718336188837418-0.3821764564027066435713428239148651415110i",
                "-0.208149705660880753184562545043304137887-2.223933997628376128108617822195260829291i",
                "-0.6190829222220442268337939993060742471432-1.2340052506884796031095691975154783931566i",
                "-0.8950122679356556222426321310933554690682-0.1921567118951016657419284832825046502150i",
                "-0.1261919022008915636485879466232029064825-0.8229360599606896062387667237310913439637i",
                "0.0327411591965760885745077550885054265160-0.4465758813540163411450525861754808277921i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(1 - double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.550183551240869436567705780226583035675",
                "1.624838898635177482810707382283843714659",
                "0.2537597545660558629373183819678928568856",
                "0.5884541528410470984016792627677679053480-1.2041875219250864403526265518377459885387i",
                "-8.547285323762208096558487265018928902238+0.048072687634043540547844060331060782881i",
                "-2.592886175491196978167660649169527169642+0.180489972066962026629620881398185425545i",
                "-0.9698223821961251506536790491267346088198+0.5542303961176891718851528116605853224851i",
                "-1.414801599478460447425793763192031774334-1.099298042201829769035811837132655443771i",
                "-0.3549534413309311974365969257690915952798-0.8415652386102599639944660479537972157651i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2PlusExpM16Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.550397500044793862393393885345466599632",
                "1.624870116422659762286344853553281356619",
                "0.2537627577191850718589996604764836043468",
                "0.5884578627962372257805630990076037463480-1.2042117499612484807003896420986829806847i",
                "-8.547511636157966173217954178218318033295+0.048275869628933799127403289400289440329i",
                "-2.592919393471137138528747784076527036099+0.180547748824681398419072230372860887591i",
                "-0.9698188310133243532146608942722819900403+0.5542469025582847956864781631376636531377i",
                "-1.414838145595104139475146077076976313384-1.099295724988492983787276530897676436277i",
                "-0.3549655329871949982970948665454995278892-0.8415738520183002837660727685719039225172i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 + double.ScaleB(1, -16), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2PlusExpM12Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.553607581009682318351898219051359605351",
                "1.625338476938186087916249148236120353700",
                "0.2538078117012088948738797439848285906405",
                "0.5885135087501393360262804536852351773342-1.2045752408644029326657258810975253663606i",
                "-8.550906645376241421237832585943651289734+0.051324911574460917431461780063004876907i",
                "-2.593417621041812861026638843040054652558+0.181414614359507095774509033165784334629i",
                "-0.9697655326867874287141984781087383035576+0.5544945192146772232999299978231219557400i",
                "-1.415386428851729797870048040792540148262-1.099260909668014742361781649075106225680i",
                "-0.3551469450333695697549819496071604002297-0.8417030634035348584853620392997141289382i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 + double.ScaleB(1, -12), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2PlusExpM8Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.605186003209314681449213469433871952042",
                "1.632856194656740915728210811804029104649",
                "0.2545303829461673151902727923981327996259",
                "0.5894029736734624221103445449609367881602-1.2104090760152445346569666326357916612455i",
                "-8.605308882796067065926665263333194384342+0.100445259466048125713822171391827244306i",
                "-2.601378409326440481553537595493197983473+0.195339173869790859155409997889712945484i",
                "-0.9689049423023784605362205241595544447571+0.5584614941045890547665967670123948737013i",
                "-1.424182328183754167774053355225235465714-1.098689242858581287226407852487381851272i",
                "-0.3580590444478163559159390913280948417699-0.8437730610813117198228623866801278741993i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 + double.ScaleB(1, -8), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2PlusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "8.488816955039522776484168306159830984029",
                "1.759503464462799230214929358767642054037",
                "0.2665405339745429557350631793465931459201",
                "0.603368678097830917806353114804768030548-1.308516450421869064945448311049355417440i",
                "-9.495186302370863093221538187758679107191+0.976993821399441609109472019975728764750i",
                "-2.725501509005642133794475046546049374883+0.432541017088870772191073089374975752724i",
                "-0.9530590254693324113521874181352703245533+0.6232241490367073062140625776173953886359i",
                "-1.571034381212324343728370133999833967707-1.085562065785112411166789237920428771917i",
                "-0.4071707583741840442963837355895266971213-0.8775536193054396219309899569163903639352i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 + double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2MinusExpM16Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.549969609509742275052552391155845922753",
                "1.624807681628436555632535528126329027945",
                "0.2537567514686295563272158407600167949169",
                "0.5884504428577623188654865146454127683596-1.2041632944751798902109377319754034840452i",
                "-8.547059014061213991755920290847574697843+0.047869516570264716991107862361155293763i",
                "-2.592852957160424862944592492321577572346+0.180432197093772349398648739118839673503i",
                "-0.9698259331240941614502725531802443091752+0.5542138898441637110462092285521364288674i",
                "-1.414765054124303364469718837459276120668-1.099300358939284205329608880471050209435i",
                "-0.3549413499846423016905425280690482334982-0.8415566252877624638851599702457354505940i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 - double.ScaleB(1, -16), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2MinusExpM12Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.546761332108207355659822614240055351233",
                "1.624339520201958431978132631691824132893",
                "0.2537117116908459544135198091094393942052",
                "0.5883947897397229252200749414187396856910-1.2037999530671775811273842814307103069764i",
                "-8.543664692007679428333986996304252437162+0.044823262058221981865229284811487652064i",
                "-2.592354640127590819149283423926782380006+0.179565786614019988581763147332889858639i",
                "-0.9698791664685249285845118686036195906087+0.5539663157906621244256345047569793182850i",
                "-1.414216965301763422452499275583458112926-1.099335052909759519795195431205405061314i",
                "-0.3547600169820695392761008573679672734332-0.8414274357159464626055583174955014085729i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 - double.ScaleB(1, -12), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2MinusExpM8Test() {
            Complex<Pow2.N4>[] expecteds = [
                "7.495644625478380146696905363168611386149",
                "1.616872769500342581563561258219644988961",
                "0.2529927767400117384122692348280849911348",
                "0.5875034907328154255827362527241969005939-1.1980043888155514960869510895467561910294i",
                "-8.489438364033819286332191474857436653793-0.003583498375342941677371003746330079748i",
                "-2.584370948886125156662028360014722347039+0.165757721282732058951034912834550468001i",
                "-0.9707231213968255527740955719965557172922+0.5500102471939863592023139622260362604416i",
                "-1.405470841121738257713269853661740411956-1.099875653860005856260780252959524739140i",
                "-0.3518681527955394573455214429024388607077-0.8393630222444596007855414801384063401469i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 - double.ScaleB(1, -8), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu2MinusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "6.730432685573963830226935949733345110740",
                "1.503287521306175352768951162925216722985",
                "0.2419140826108273665812386204748892505985",
                "0.5730640340814938395386491286888242062323-1.1097035609396099823667412474874994763527i",
                "-7.644283491902177700248717738371098786558-0.697159761222871704759189712754588342601i",
                "-2.454342451056385863063314418040290040854-0.041616378887077402966791702059145881014i",
                "-0.9823079420237559046219928132856658427446+0.4880365337206342217554135652124079150307i",
                "-1.271362928289379082084686507803115597941-1.105033026681955588793302877120349745878i",
                "-0.3079409741140201699219154718861379457591-0.8070101273268822780858066703257020401160i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(2 - double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu3Test() {
            Complex<Pow2.N4>[] expecteds = [
                "62.05790952993025638623809164498675579694",
                "7.101262824737944505980369530670992197899",
                "0.6473853909486341531592355709711967376584",
                "0.332677729438730329777357672654762338554-5.196781653227249138881030465336938362563i",
                "0.00402709765839154318535140609618012646+66.06689916071114919963370549316568353877i",
                "0.030730044575763843635341861745771183462+9.144418471821216423427361814970841191834i",
                "0.202543582639388725992751918520417833086+1.771518614079819366079653417798140805377i",
                "-5.938306666410040587830424961054633514416+1.926193700300242412131163961591575706493i",
                "-2.368469054358641974249648601260576386094-1.432943068359846897891781252887983239080i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(3, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu3PlusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "71.81207810354501175814051959110243863685",
                "7.887541183344830150411785863033021716886",
                "0.6933204607850969621477321268393787144968",
                "0.217588184291212680591481460958372039385-5.737305627949666811119814002523331266434i",
                "7.48261903152171888966825680558971345638+75.93748574713018619592942245309855947421i",
                "1.014849166012871417419363063020946346043+10.026261299711730024696749391074321284946i",
                "0.368447697735738014067806056637046717667+1.823662592174277563483766495384510252499i",
                "-6.365581077378259060609018901397981280396+2.530675538913674063159306534410719136138i",
                "-2.634354202511530336041997269650749351523-1.444095499507186091396845990015175666823i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(3 + double.ScaleB(1, -4), z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-37);
            }
        }

        [TestMethod()]
        public void BesselKNu3MinusExpM4Test() {
            Complex<Pow2.N4>[] expecteds = [
                "53.70853416240554599768910685604137246312",
                "6.402125286805117859887059944417952634854",
                "0.6051650109947971257237898597355814991997",
                "0.424205587258635122993506037163349734715-4.710506748141813040193146840593791037864i",
                "-5.61113785513813646506828281889136379304+57.01926306461350104881894465018075476653i",
                "-0.780346390409469163330111844785264536071+8.276571747155015616260774172390645515152i",
                "0.050558516148073182408481651505684866304+1.711064925775925399463182592458472920103i",
                "-5.522406823763491466139090952791116295513+1.405487826379354385594141990922024903544i",
                "-2.127476599845172404440453895446934546117-1.413522365703936675712591783152136260608i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = AmosPowerSeries<Pow2.N4>.BesselK(3 - double.ScaleB(1, -4), z);

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