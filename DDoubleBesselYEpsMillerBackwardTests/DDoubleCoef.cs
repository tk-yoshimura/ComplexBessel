using DoubleDouble;
using MultiPrecision;
using System.Collections.ObjectModel;

namespace DDoubleBesselYEpsMillerBackwardTests {
    public static class DDoubleCoef {

        public static (ddouble eta0, ddouble xi1) BesselYEta0Xi1(ddouble alpha, ddouble x) {
            ddouble s = ddouble.Pow(ddouble.Ldexp(1 / x, 1), alpha), sqs = s * s;

            ddouble rcot = 1d / ddouble.TanPI(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            ddouble r = ddouble.Ldexp(ddouble.RcpPI * sqs, 1);
            ddouble p = sqs * rsqgamma * ddouble.RcpPI;

            ddouble eta0 = rcot - p / alpha;
            ddouble xi1 = rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d));

            return (eta0, xi1);
        }

        public static (ddouble eta0, ddouble xi1) BesselYEta0Xi1Eps(ddouble alpha, ddouble x) {
            ddouble lnx = ddouble.Log(x), g = ddouble.EulerGamma;

            ddouble eta0, xi1;

            {
                ddouble r0 = eta0_coef[0] + lnx * eta0_coef[1]
                           + g * eta0_coef[2];
                ddouble r1 = eta0_coef[3] + lnx * (eta0_coef[4] + lnx * eta0_coef[5])
                           + g * (eta0_coef[6] + lnx * eta0_coef[7]
                           + g * eta0_coef[8]);
                ddouble r2 = eta0_coef[9] + lnx * (eta0_coef[10] + lnx * (eta0_coef[11] + lnx * eta0_coef[12]))
                           + g * (eta0_coef[13] + lnx * (eta0_coef[14] + lnx * eta0_coef[15])
                           + g * (eta0_coef[16] + lnx * eta0_coef[17]
                           + g * eta0_coef[18]));
                ddouble r3 = eta0_coef[19] + lnx * (eta0_coef[20] + lnx * (eta0_coef[21] + lnx * (eta0_coef[22] + lnx * eta0_coef[23])))
                           + g * (eta0_coef[24] + lnx * (eta0_coef[25] + lnx * (eta0_coef[26] + lnx * eta0_coef[27]))
                           + g * (eta0_coef[28] + lnx * (eta0_coef[29] + lnx * eta0_coef[30])
                           + g * (eta0_coef[31] + lnx * eta0_coef[32]
                           + g * eta0_coef[33])));
                ddouble r4 = eta0_coef[34] + lnx * (eta0_coef[35] + lnx * (eta0_coef[36] + lnx * (eta0_coef[37] + lnx * (eta0_coef[38] + lnx * eta0_coef[39]))))
                           + g * (eta0_coef[40] + lnx * (eta0_coef[41] + lnx * (eta0_coef[42] + lnx * (eta0_coef[43] + lnx * eta0_coef[44])))
                           + g * (eta0_coef[45] + lnx * (eta0_coef[46] + lnx * (eta0_coef[47] + lnx * eta0_coef[48]))
                           + g * (eta0_coef[49] + lnx * (eta0_coef[50] + lnx * eta0_coef[51])
                           + g * (eta0_coef[52] + lnx * eta0_coef[53]
                           + g * eta0_coef[54]))));
                ddouble r5 = eta0_coef[55] + lnx * (eta0_coef[56] + lnx * (eta0_coef[57] + lnx * (eta0_coef[58] + lnx * (eta0_coef[59] + lnx * (eta0_coef[60] + lnx * eta0_coef[61])))))
                           + g * (eta0_coef[62] + lnx * (eta0_coef[63] + lnx * (eta0_coef[64] + lnx * (eta0_coef[65] + lnx * (eta0_coef[66] + lnx * eta0_coef[67]))))
                           + g * (eta0_coef[68] + lnx * (eta0_coef[69] + lnx * (eta0_coef[70] + lnx * (eta0_coef[71] + lnx * eta0_coef[72])))
                           + g * (eta0_coef[73] + lnx * (eta0_coef[74] + lnx * (eta0_coef[75] + lnx * eta0_coef[76]))
                           + g * (eta0_coef[77] + lnx * (eta0_coef[78] + lnx * eta0_coef[79])
                           + g * (eta0_coef[80] + lnx * eta0_coef[81]
                           + g * eta0_coef[82])))));

                eta0 = r0 + alpha * (r1 + alpha * (r2 + alpha * (r3 + alpha * (r4 + alpha * r5))));
            }

            {
                ddouble r0 = xi1_coef[0] + lnx * xi1_coef[1]
                           + g * xi1_coef[2];
                ddouble r1 = xi1_coef[3] + lnx * (xi1_coef[4] + lnx * xi1_coef[5])
                           + g * (xi1_coef[6] + lnx * xi1_coef[7]
                           + g * xi1_coef[8]);
                ddouble r2 = xi1_coef[9] + lnx * (xi1_coef[10] + lnx * (xi1_coef[11] + lnx * xi1_coef[12]))
                           + g * (xi1_coef[13] + lnx * (xi1_coef[14] + lnx * xi1_coef[15])
                           + g * (xi1_coef[16] + lnx * xi1_coef[17]
                           + g * xi1_coef[18]));
                ddouble r3 = xi1_coef[19] + lnx * (xi1_coef[20] + lnx * (xi1_coef[21] + lnx * (xi1_coef[22] + lnx * xi1_coef[23])))
                           + g * (xi1_coef[24] + lnx * (xi1_coef[25] + lnx * (xi1_coef[26] + lnx * xi1_coef[27]))
                           + g * (xi1_coef[28] + lnx * (xi1_coef[29] + lnx * xi1_coef[30])
                           + g * (xi1_coef[31] + lnx * xi1_coef[32]
                           + g * xi1_coef[33])));
                ddouble r4 = xi1_coef[34] + lnx * (xi1_coef[35] + lnx * (xi1_coef[36] + lnx * (xi1_coef[37] + lnx * (xi1_coef[38] + lnx * xi1_coef[39]))))
                           + g * (xi1_coef[40] + lnx * (xi1_coef[41] + lnx * (xi1_coef[42] + lnx * (xi1_coef[43] + lnx * xi1_coef[44])))
                           + g * (xi1_coef[45] + lnx * (xi1_coef[46] + lnx * (xi1_coef[47] + lnx * xi1_coef[48]))
                           + g * (xi1_coef[49] + lnx * (xi1_coef[50] + lnx * xi1_coef[51])
                           + g * (xi1_coef[52] + lnx * xi1_coef[53]
                           + g * xi1_coef[54]))));
                ddouble r5 = xi1_coef[55] + lnx * (xi1_coef[56] + lnx * (xi1_coef[57] + lnx * (xi1_coef[58] + lnx * (xi1_coef[59] + lnx * (xi1_coef[60] + lnx * xi1_coef[61])))))
                           + g * (xi1_coef[62] + lnx * (xi1_coef[63] + lnx * (xi1_coef[64] + lnx * (xi1_coef[65] + lnx * (xi1_coef[66] + lnx * xi1_coef[67]))))
                           + g * (xi1_coef[68] + lnx * (xi1_coef[69] + lnx * (xi1_coef[70] + lnx * (xi1_coef[71] + lnx * xi1_coef[72])))
                           + g * (xi1_coef[73] + lnx * (xi1_coef[74] + lnx * (xi1_coef[75] + lnx * xi1_coef[76]))
                           + g * (xi1_coef[77] + lnx * (xi1_coef[78] + lnx * xi1_coef[79])
                           + g * (xi1_coef[80] + lnx * xi1_coef[81]
                           + g * xi1_coef[82])))));

                xi1 = r0 + alpha * (r1 + alpha * (r2 + alpha * (r3 + alpha * (r4 + alpha * r5))));
            }

            return (eta0, xi1);
        }

        private static readonly ReadOnlyCollection<ddouble> eta0_coef = new([
            (-1, -2, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 0, 0xF03677ABD9FF2936uL, 0x3565C32E09B19596uL),
            (+1, -1, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -1, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (-1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, -1, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (+1, 0, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (-1, -1, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 0, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (-1, 0, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (+1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, -1, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (+1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 0, 0xBEC615D3618A833AuL, 0x44E78617E13792F3uL),
            (+1, 0, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (-1, 0, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (+1, -1, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 0, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (-1, 1, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (+1, 0, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (-1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 0, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (+1, 0, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (-1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -1, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, -2, 0x91FDA4152F0A06F0uL, 0x67C83E946B398D5CuL),
            (+1, 0, 0xCD2801D8B434667FuL, 0xE1831DD92E39839BuL),
            (-1, 0, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (+1, 0, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (-1, -2, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (+1, 0, 0xCD2801D8B434667FuL, 0xE1831DD92E39839BuL),
            (-1, 1, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (+1, 1, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (-1, 0, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 0, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (+1, 1, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (-1, 0, 0xE1EE4C7BF4B4DDD9uL, 0x6C10B05CABFDE7BBuL),
            (+1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 0, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (-1, 0, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (+1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, -2, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 0, 0x94C825A9B12DD9ADuL, 0x5D12B0E09466609DuL),
            (+1, -1, 0x91FDA4152F0A06F0uL, 0x67C83E946B398D5CuL),
            (-1, 0, 0xCD2801D8B434667FuL, 0xE1831DD92E39839BuL),
            (+1, -1, 0xD0EFA8D243C36A9CuL, 0xE25A4A80C14CAE10uL),
            (-1, -1, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (+1, -4, 0xF0FE2F732727530AuL, 0x0CDE99FC7330F72DuL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (+1, -1, 0x91FDA4152F0A06F0uL, 0x67C83E946B398D5CuL),
            (-1, 1, 0xCD2801D8B434667FuL, 0xE1831DD92E39839BuL),
            (+1, 1, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (-1, 1, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (+1, -1, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -3, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 0, 0xCD2801D8B434667FuL, 0xE1831DD92E39839BuL),
            (+1, 1, 0x9CB3BE9DB2D28FF5uL, 0xA9C3B7E090F9828CuL),
            (-1, 1, 0xD457CBD4DCC7FA26uL, 0x0F18584CBCE85A15uL),
            (+1, 0, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -1, 0xD0EFA8D243C36A9CuL, 0xE25A4A80C14CAE10uL),
            (-1, 1, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (+1, 0, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -1, 0x90DDCA29292012CFuL, 0xFC94BFD712D9A8ABuL),
            (-1, -1, 0x8D8FDD38932FFC19uL, 0x5F659033289AE6B9uL),
            (+1, -1, 0x969EDDA7F87893E6uL, 0x480B203DC7FE9A7CuL),
            (-1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -4, 0xF0FE2F732727530AuL, 0x0CDE99FC7330F72DuL),
            (-1, -3, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
        ]);

        private static readonly ReadOnlyCollection<ddouble> xi1_coef = new([
            (-1, 0, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0xEDB4603E47865401uL, 0x3445CE9CEBCC17E2uL),
            (+1, 1, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, 1, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (-1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 2, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (+1, 2, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (-1, 1, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 2, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (-1, 2, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (+1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (+1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 2, 0xEF6AD8880470F02AuL, 0xE4DE7223A19CEC2FuL),
            (+1, 3, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (-1, 2, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (+1, 0, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 3, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (-1, 3, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (+1, 2, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (-1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 2, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (+1, 2, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (-1, 0, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, 0, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 3, 0x87DCD2CF31332DB5uL, 0x03DD297CA6D818C0uL),
            (+1, 3, 0xD95E534E4294DC2CuL, 0x4FD4F458CF1637E5uL),
            (-1, 3, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (+1, 1, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (-1, -1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (+1, 3, 0xD95E534E4294DC2CuL, 0x4FD4F458CF1637E5uL),
            (-1, 4, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (+1, 3, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (-1, 1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 3, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (+1, 3, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (-1, 2, 0x89F854D6244F420BuL, 0x5917D8002599E8CFuL),
            (+1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 1, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (-1, 1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (+1, -1, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, -1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 3, 0xABD844C29C147F8AuL, 0x30D5F36AEFA24978uL),
            (+1, 4, 0x87DCD2CF31332DB5uL, 0x03DD297CA6D818C0uL),
            (-1, 3, 0xD95E534E4294DC2CuL, 0x4FD4F458CF1637E5uL),
            (+1, 2, 0xC21FD7B122CEA821uL, 0x260C137554C3BA7DuL),
            (-1, 0, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (+1, -2, 0x932B0528AF4379A5uL, 0xB45DB3335B4ED632uL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (+1, 4, 0x87DCD2CF31332DB5uL, 0x03DD297CA6D818C0uL),
            (-1, 4, 0xD95E534E4294DC2CuL, 0x4FD4F458CF1637E5uL),
            (+1, 4, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (-1, 2, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (+1, 0, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -3, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 3, 0xD95E534E4294DC2CuL, 0x4FD4F458CF1637E5uL),
            (+1, 4, 0x9197E1C4DA1AFE18uL, 0xDC890E97FF92CBDEuL),
            (-1, 3, 0xAAAF175D91B8BDEFuL, 0x9D590319162D639CuL),
            (+1, 1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, 2, 0xC21FD7B122CEA821uL, 0x260C137554C3BA7DuL),
            (-1, 2, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (+1, 1, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -1, 0x90DDCA29292012CFuL, 0xFC94BFD712D9A8ABuL),
            (-1, 0, 0xE3941F276CF65294uL, 0xD1CC042172E72F7BuL),
            (+1, 0, 0xB7F5C672DB14580FuL, 0x2175200032228BBEuL),
            (-1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -2, 0x932B0528AF4379A5uL, 0xB45DB3335B4ED632uL),
            (-1, -3, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
        ]);
    }
}
