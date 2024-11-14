
using Algebra;
using Regression;

namespace ComplexBesselPowerSeriesRange {
    internal class Program {
        static void Main() {
            List<double> xs = [], ys = [], nus = [];

            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamReader sr = new($"../../../../results_ddouble/besselj_nu{nu}_powerseries_range.csv");

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
                }
            }

            {
                Vector x = xs.ToArray();
                Vector nu = nus.ToArray();
                Vector y = ys.ToArray();

                RobustRegressor fitter = new([nu, nu * nu, x, x * nu], y, intercept: 7.5);

                Vector param = fitter.Fit();

                Console.WriteLine($"{param:e2}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
