using DoubleDouble;
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
            ddouble lnx = ddouble.Log(x);

            ddouble eta0, xi1;

            {
                ddouble r0 = eta0_coef[0] + lnx * eta0_coef[1];
                ddouble r1 = eta0_coef[2] + lnx * (eta0_coef[3] + lnx * eta0_coef[4]);
                ddouble r2 = eta0_coef[5] + lnx * (eta0_coef[6] + lnx * (eta0_coef[7] + lnx * eta0_coef[8]));
                ddouble r3 = eta0_coef[9] + lnx * (eta0_coef[10] + lnx * (eta0_coef[11] + lnx * (eta0_coef[12]
                           + lnx * eta0_coef[13])));
                ddouble r4 = eta0_coef[14] + lnx * (eta0_coef[15] + lnx * (eta0_coef[16] + lnx * (eta0_coef[17]
                           + lnx * (eta0_coef[18] + lnx * eta0_coef[19]))));
                ddouble r5 = eta0_coef[20] + lnx * (eta0_coef[21] + lnx * (eta0_coef[22] + lnx * (eta0_coef[23]
                           + lnx * (eta0_coef[24] + lnx * (eta0_coef[25] + lnx * eta0_coef[26])))));
                ddouble r6 = eta0_coef[27] + lnx * (eta0_coef[28] + lnx * (eta0_coef[29] + lnx * (eta0_coef[30]
                           + lnx * (eta0_coef[31] + lnx * (eta0_coef[32] + lnx * (eta0_coef[33] + lnx * eta0_coef[34]))))));
                ddouble r7 = eta0_coef[35] + lnx * (eta0_coef[36] + lnx * (eta0_coef[37] + lnx * (eta0_coef[38]
                           + lnx * (eta0_coef[39] + lnx * (eta0_coef[40] + lnx * (eta0_coef[41] + lnx * (eta0_coef[42]
                           + lnx * eta0_coef[43])))))));

                eta0 = r0 + alpha * (r1 + alpha * (r2 + alpha * (r3 + alpha * (r4 + alpha * (r5 + alpha * (r6 + alpha * r7))))));
            }

            {
                ddouble r0 = xi1_coef[0] + lnx * xi1_coef[1];
                ddouble r1 = xi1_coef[2] + lnx * (xi1_coef[3] + lnx * xi1_coef[4]);
                ddouble r2 = xi1_coef[5] + lnx * (xi1_coef[6] + lnx * (xi1_coef[7] + lnx * xi1_coef[8]));
                ddouble r3 = xi1_coef[9] + lnx * (xi1_coef[10] + lnx * (xi1_coef[11] + lnx * (xi1_coef[12]
                           + lnx * xi1_coef[13])));
                ddouble r4 = xi1_coef[14] + lnx * (xi1_coef[15] + lnx * (xi1_coef[16] + lnx * (xi1_coef[17]
                           + lnx * (xi1_coef[18] + lnx * xi1_coef[19]))));
                ddouble r5 = xi1_coef[20] + lnx * (xi1_coef[21] + lnx * (xi1_coef[22] + lnx * (xi1_coef[23]
                           + lnx * (xi1_coef[24] + lnx * (xi1_coef[25] + lnx * xi1_coef[26])))));
                ddouble r6 = xi1_coef[27] + lnx * (xi1_coef[28] + lnx * (xi1_coef[29] + lnx * (xi1_coef[30]
                           + lnx * (xi1_coef[31] + lnx * (xi1_coef[32] + lnx * (xi1_coef[33] + lnx * xi1_coef[34]))))));
                ddouble r7 = xi1_coef[35] + lnx * (xi1_coef[36] + lnx * (xi1_coef[37] + lnx * (xi1_coef[38]
                           + lnx * (xi1_coef[39] + lnx * (xi1_coef[40] + lnx * (xi1_coef[41] + lnx * (xi1_coef[42]
                           + lnx * xi1_coef[43])))))));

                xi1 = r0 + alpha * (r1 + alpha * (r2 + alpha * (r3 + alpha * (r4 + alpha * (r5 + alpha * (r6 + alpha * r7))))));
            }

            return (eta0, xi1);
        }

        public static (ddouble eta0, ddouble xi1) BesselYEta0Xi1EpsR2(ddouble alpha, ddouble x) {
            const int N = 7;

            ReadOnlyCollection<ddouble> eta0_coef = new(DDoubleCoef.eta0_coef.Reverse().ToArray());
            ReadOnlyCollection<ddouble> xi1_coef = new(DDoubleCoef.xi1_coef.Reverse().ToArray());

            ddouble lnx = ddouble.Log(x);

            ddouble eta0 = 0d, xi1 = 0d;
            for (int i = N, k = 0; i >= 0; i--) {
                ddouble s = eta0_coef[k], t = xi1_coef[k];
                k++;

                for (int j = i; j >= 0; j--, k++) {
                    s = eta0_coef[k] + lnx * s;
                    t = xi1_coef[k] + lnx * t;
                }

                eta0 = s + alpha * eta0;
                xi1 = t + alpha * xi1;
            }

            return (eta0, xi1);
        }

        private static readonly ReadOnlyCollection<ddouble> eta0_coef = new([
            (-1, -4, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 0, 0xCA28399BC43B5702uL, 0xF0B44645EF33148BuL),
            (+1, -3, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -3, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (+1, 0, 0x883B4FB4B14055BFuL, 0x85B55E7C87EB57FFuL),
            (-1, -3, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 0, 0x9F9A4CA27EFFD464uL, 0xD17DF0A9838B4A3CuL),
            (-1, -2, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (-1, 0, 0x883B4FB4B14055BFuL, 0x85B55E7C87EB57FFuL),
            (+1, -4, 0xC988F11328ABF05AuL, 0xB14D0ED3FA8AE911uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -2, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (+1, 0, 0x8ED06F76EF1F08D4uL, 0xFAAFF2FC72E0F22EuL),
            (+1, -2, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (+1, -1, 0xB5A46A4641AB1CFFuL, 0x5CF1D350B539CAA9uL),
            (-1, -5, 0xC988F11328ABF05AuL, 0xB14D0ED3FA8AE911uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 0, 0xA1EFBF6811197CC5uL, 0xEA0BD6E32998A1A1uL),
            (-1, -1, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (-1, 0, 0x8ED06F76EF1F08D4uL, 0xFAAFF2FC72E0F22EuL),
            (-1, -3, 0xB59DD3B1A81E56E5uL, 0x59E099DAB4D8985AuL),
            (-1, -2, 0xB5A46A4641AB1CFFuL, 0x5CF1D350B539CAA9uL),
            (+1, -6, 0xA13A5A75BA2326AEuL, 0xF43DA5766208BA74uL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (+1, -1, 0xA2EEA6FBB819E07EuL, 0x95CCA4D4A4C628CEuL),
            (+1, 0, 0x9E1267360C8533E2uL, 0xF1588EFC26ADAEDBuL),
            (+1, -1, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (+1, -1, 0xBE6B3F493ED40BC6uL, 0xA39543FB43D6983EuL),
            (+1, -4, 0xB59DD3B1A81E56E5uL, 0x59E099DAB4D8985AuL),
            (+1, -3, 0x915055050155B0CCuL, 0x4A5B0F73C42E3BBAuL),
            (-1, -8, 0xD6F8789CF82EDE3EuL, 0x9AFCDC9DD80BA345uL),
            (+1, -7, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
            (-1, 0, 0xB0EDF45B2A0853AFuL, 0x2A42C979CA5D0596uL),
            (-1, 0, 0xA2EEA6FBB819E07EuL, 0x95CCA4D4A4C628CEuL),
            (-1, 0, 0x9E1267360C8533E2uL, 0xF1588EFC26ADAEDBuL),
            (-1, -1, 0x8E8A6C2954F709B2uL, 0xDFC03FBAF929AA92uL),
            (-1, -2, 0xBE6B3F493ED40BC6uL, 0xA39543FB43D6983EuL),
            (-1, -5, 0x914B0FC1534B78B7uL, 0x7B1A14AEF713AD15uL),
            (-1, -5, 0xC1C0715C01C79665uL, 0xB87969EFB03DA4F9uL),
            (+1, -10, 0xF5AE40B364C7D96CuL, 0x1ED7D78FD2567173uL),
            (-1, -9, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
        ]);

        private static readonly ReadOnlyCollection<ddouble> xi1_coef = new([
            (-1, -1, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0xABA41964255F42B5uL, 0x773880C3A34BE05AuL),
            (+1, 0, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (+1, 1, 0xD13DA106DF235947uL, 0xC0976A7F9B5A5829uL),
            (-1, 0, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 1, 0xF03C087CD2E2F132uL, 0xABA01CBF4A537060uL),
            (+1, 2, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (-1, 1, 0xD13DA106DF235947uL, 0xC0976A7F9B5A5829uL),
            (+1, -1, 0xF27DCD6022C59A43uL, 0x5108C19D1B97DA23uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 1, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (+1, 2, 0xC422FE094F2AC935uL, 0x818D2129A54607CEuL),
            (-1, 2, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (+1, 1, 0x8B7E6B59EA1790DAuL, 0x8064F1AA6791901BuL),
            (-1, -2, 0xF27DCD6022C59A43uL, 0x5108C19D1B97DA23uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 2, 0x8F45E37E7DC47E26uL, 0x8FE6337E8ACC854DuL),
            (+1, 2, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (-1, 2, 0xC422FE094F2AC935uL, 0x818D2129A54607CEuL),
            (+1, 1, 0xB3D9F038F945D22EuL, 0xD78DC5578966D089uL),
            (-1, 0, 0x8B7E6B59EA1790DAuL, 0x8064F1AA6791901BuL),
            (+1, -3, 0xC1FE3DE68237AE9CuL, 0x40D3CE174946481CuL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (-1, 1, 0xD9277CDD4B4A0DEAuL, 0x740B8DE15C3A3708uL),
            (+1, 2, 0xF5188116761D8AE2uL, 0xE71C9F4A8A782581uL),
            (-1, 2, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (+1, 2, 0x82C1FEB0DF71DB79uL, 0x0108C0C66E2EAFDEuL),
            (-1, 0, 0xB3D9F038F945D22EuL, 0xD78DC5578966D089uL),
            (+1, -2, 0xDF30ABC31025B490uL, 0xCD6E4F770C1C19C6uL),
            (-1, -4, 0x8154294456CFC9BDuL, 0x808D340F862EDABDuL),
            (+1, -7, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
            (-1, 2, 0x9833B38CBAB4864DuL, 0x9007F69410F14DABuL),
            (+1, 2, 0xD9277CDD4B4A0DEAuL, 0x740B8DE15C3A3708uL),
            (-1, 2, 0xF5188116761D8AE2uL, 0xE71C9F4A8A782581uL),
            (+1, 2, 0x83112A774F6C8904uL, 0x4EB45C7FA1DD88D1uL),
            (-1, 1, 0x82C1FEB0DF71DB79uL, 0x0108C0C66E2EAFDEuL),
            (+1, -1, 0x8FE18CFA6104A825uL, 0x793E37793AB8A6D4uL),
            (-1, -3, 0x94CB1D2CB56E7860uL, 0x88F434FA0812BBD9uL),
            (+1, -6, 0x93CDE604F57FC1FDuL, 0x2533A93650358C46uL),
            (-1, -9, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
        ]);
    }
}
