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
        public void BesselJNu1p25Test() {
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
        public void BesselYNu1p25Test() {
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
        public void BesselINu1p25Test() {
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
        public void BesselKNu1p25Test() {
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

        [TestMethod()]
        public void BesselJNuM1p75Test() {
            HankelExpansion<Pow2.N4> hankel = new(-1.75);

            Complex<Pow2.N4>[] expecteds = [
                "-147.2999139438063945323394630815032224390-10.9056339300812418071382865800802817172i",
                "-147.2999139438063945323394630815032224390+10.9056339300812418071382865800802817172i",
                "-96.4453203127618438412463729175939687224+111.8682157229589337900597395166247943914i",
                "-96.4453203127618438412463729175939687224-111.8682157229589337900597395166247943914i",
                "1.397508963018419728866125463198502317718e26-2.689047103063237008011487557828466479775e26i",
                "-9.13255376986751184684287116150543446463e25+2.889631506025360495813279339858706432795e26i",
                "1.397508963018419728866125463198502317718e26+2.689047103063237008011487557828466479775e26i",
                "-9.13255376986751184684287116150543446463e25-2.889631506025360495813279339858706432795e26i"
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
        public void BesselYNuM1p75Test() {
            HankelExpansion<Pow2.N4> hankel = new(-1.75);

            Complex<Pow2.N4>[] expecteds = [
                "-10.9056322416899217958437967766265351684+147.2998805626433483749635926306468094384i",
                "-10.9056322416899217958437967766265351684-147.2998805626433483749635926306468094384i",
                "111.8681933127851316348824917693777891420+96.4453451106815493494087232567175383028i",
                "111.8681933127851316348824917693777891420-96.4453451106815493494087232567175383028i",
                "2.689047103063237008011487557828466479775e26+1.397508963018419728866125463198502317718e26i",
                "-2.889631506025360495813279339858706432795e26-9.13255376986751184684287116150543446463e25i",
                "2.689047103063237008011487557828466479775e26-1.397508963018419728866125463198502317718e26i",
                "-2.889631506025360495813279339858706432795e26+9.13255376986751184684287116150543446463e25i"
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
        public void BesselINuM1p75Test() {
            HankelExpansion<Pow2.N4> hankel = new(-1.75);

            Complex<Pow2.N4>[] expecteds = [
                "-2.62076152242149906566470469518604540524e25-3.019159107207759513123694167620311528728e26i",
                "-2.62076152242149906566470469518604540524e25+3.019159107207759513123694167620311528728e26i",
                "1.949552053750027174411954077068820097967e26+2.320183702625431595125833192684946828683e26i",
                "1.949552053750027174411954077068820097967e26-2.320183702625431595125833192684946828683e26i",
                "140.2607810578372943130995226449327728200+46.2937446780354935621573347839373412787i",
                "131.9139702088759482017693824281792414315+66.4447286321608627589382194728601568628i",
                "140.2607810578372943130995226449327728200-46.2937446780354935621573347839373412787i",
                "131.9139702088759482017693824281792414315-66.4447286321608627589382194728601568628i"
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
        public void BesselKNuM1p75Test() {
            HankelExpansion<Pow2.N4> hankel = new(-1.75);

            Complex<Pow2.N4>[] expecteds = [
                "-5.35222404565129229070775331686841841211e-30+2.500563373120195349488734686337018462950e-29i",
                "-5.35222404565129229070775331686841841211e-30-2.500563373120195349488734686337018462950e-29i",
                "-9.484968071222616258294302283203169493183e26-8.23336514565018365948586434406526860170e25i",
                "-9.484968071222616258294302283203169493183e26+8.23336514565018365948586434406526860170e25i",
                "-208.7422488243166520251860480998194073164-414.4199122853585266887857775703335245101i",
                "-0.00001761577059679953277111031956241000418243-0.00004945855291014150741401368246876143340060i",
                "-208.7422488243166520251860480998194073164+414.4199122853585266887857775703335245101i",
                "-0.00001761577059679953277111031956241000418243+0.00004945855291014150741401368246876143340060i"
            ];

            foreach ((Complex<Pow2.N4> z, Complex<Pow2.N4> expected) in zs_mini.Zip(expecteds)) {
                Complex<Pow2.N4> actual = hankel.BesselK(z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 1e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu1p375Test() {
            HankelExpansion<Pow2.N4> hankel = new(1.375);

            Complex<Pow2.N4>[] expecteds = [
                "-19.3844491315920837390199087861121994352-146.5832933973515097634608469293879648659i",
                "-19.3844491315920837390199087861121994352+146.5832933973515097634608469293879648659i",
                "-128.0071970497257571256627952069408871441-74.0038936463615894478416609665561282765i",
                "-128.0071970497257571256627952069408871441+74.0038936463615894478416609665561282765i",
                "-2.388500783627850265566881132161440929619e26-1.909824495175187236085662122110732590600e26i",
                "2.678487439866468568876646704464127999900e26+1.475828794352374901455927115895960295955e26i",
                "-2.388500783627850265566881132161440929619e26+1.909824495175187236085662122110732590600e26i",
                "2.678487439866468568876646704464127999900e26-1.475828794352374901455927115895960295955e26i"
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
        public void BesselYNu1p375Test() {
            HankelExpansion<Pow2.N4> hankel = new(1.375);

            Complex<Pow2.N4>[] expecteds = [
                "-146.5833256942998126925823267036391419176+19.3844406821480655811156969072786874160i",
                "-146.5833256942998126925823267036391419176-19.3844406821480655811156969072786874160i",
                "-74.0038734805861684890514483972788984304+128.0072236547530183316730182054670692015i",
                "-74.0038734805861684890514483972788984304-128.0072236547530183316730182054670692015i",
                "1.909824495175187236085662122110732590600e26-2.388500783627850265566881132161440929619e26i",
                "-1.475828794352374901455927115895960295955e26+2.678487439866468568876646704464127999900e26i",
                "1.909824495175187236085662122110732590600e26+2.388500783627850265566881132161440929619e26i",
                "-1.475828794352374901455927115895960295955e26-2.678487439866468568876646704464127999900e26i"
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
        public void BesselINu1p375Test() {
            HankelExpansion<Pow2.N4> hankel = new(1.375);

            Complex<Pow2.N4>[] expecteds = [
                "-2.60981095641589821397947573363886292337e25-3.047007460358397106095076131919830300429e26i",
                "-2.60981095641589821397947573363886292337e25+3.047007460358397106095076131919830300429e26i",
                "-2.715194686571793667532226846517769785438e26-1.407154366007622034466960613285609780914e26i",
                "-2.715194686571793667532226846517769785438e26+1.407154366007622034466960613285609780914e26i",
                "-111.1101312101308589510940681955478744057-97.5548948736905780898205433405748677084i",
                "132.6489770521213388555665512843130872017+65.3197340654040189559932767827654910729i",
                "-111.1101312101308589510940681955478744057+97.5548948736905780898205433405748677084i",
                "132.6489770521213388555665512843130872017-65.3197340654040189559932767827654910729i"
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
        public void BesselKNu1p375Test() {
            HankelExpansion<Pow2.N4> hankel = new(1.375);

            Complex<Pow2.N4>[] expecteds = [
                "-5.27705248199702819773776225908853739551e-30+2.478879067656174354854949126828274145020e-29i",
                "-5.27705248199702819773776225908853739551e-30-2.478879067656174354854949126828274145020e-29i",
                "-9.572456252895233493029366205546736747497e26-8.19896292793433789719356411048917496803e25i",
                "-9.572456252895233493029366205546736747497e26+8.19896292793433789719356411048917496803e25i",
                "-205.2080358950215402640576240225591685438-416.7290866214760641928695100914379240115i",
                "-0.00001714958853999699091968918500505305032212-0.00004955578201524993327661980544066492052557i",
                "-205.2080358950215402640576240225591685438+416.7290866214760641928695100914379240115i",
                "-0.00001714958853999699091968918500505305032212+0.00004955578201524993327661980544066492052557i"
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