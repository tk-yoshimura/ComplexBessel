using Algebra;
using Regression;

namespace ComplexBesselMillerBackwardIteration {
    internal class Program {
        static void Main() {
            List<double> reals = [], imags = [], nus = [], ms = [];

            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamReader sr = new($"../../../../results_ddouble/besseli_nu{nu}_millerbackward_iteration.csv");

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrWhiteSpace(line)) {
                        break;
                    }

                    string[] items = line.Split(",");

                    double r = double.Parse(items[0]), i = double.Parse(items[1]);
                    int m = int.Parse(items[2]);

                    reals.Add(r);
                    imags.Add(i);
                    nus.Add(nu);
                    ms.Add(m);
                }
            }

            {
                Vector r = reals.ToArray();
                Vector i = imags.ToArray();
                Vector nu = nus.ToArray();
                Vector m = ms.ToArray();

                RobustRegressor fitter = new([r, i, r * i], m);

                Vector param = fitter.Fit();

                Vector m_fit = param[0] + r * param[1] + i * param[2] + r * i * param[3];

                for (double bias = 0; bias <= 20; bias += 0.125) {
                    Vector c = m_fit + bias - m;

                    double rate = c.Count(v => v.val > 0) / (double)c.Dim;

                    Console.WriteLine(rate);

                    if (rate > 0.999) {
                        param[0] += bias;
                        break;
                    }
                }

                Console.WriteLine($"{param:e4}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
