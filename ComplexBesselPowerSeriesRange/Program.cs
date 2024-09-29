
using Algebra;
using Regression;

namespace ComplexBesselPowerSeriesRange {
    internal class Program {
        static void Main() {
            List<double> xs = [], ys = [], nus = [], ws = [];

            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamReader sr = new($"../../../../results/besselj_nu{nu}_powerseries_range.csv");

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrWhiteSpace(line)) {
                        break;
                    }

                    string[] items = line.Split(",");

                    double r = double.Parse(items[0]), i = double.Parse(items[1]);

                    xs.Add(i);
                    ys.Add(r);
                    nus.Add(nu);
                    ws.Add(1 / (i + 1));
                }
            }

            {
                Vector x = xs.ToArray();
                Vector nu = nus.ToArray();
                Vector w = ws.ToArray();
                Vector y = ys.ToArray();

                Regressor fitter = new([x, nu, x * nu], y - 6, intercept: false);

                Vector param = fitter.Fit(w);

                Console.WriteLine($"{param:e2}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
