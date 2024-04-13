// Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
// Copyright (c) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mborgerding/kissfft

using System;
using System.Diagnostics;

namespace KissFft;

internal static partial class KissFft<T>
{
    // TODO: Use T.FusedMultiplyAdd

    internal static class Fft
    {
#pragma warning disable IDE1006 // Naming Styles

        private static void kf_bfly2(Span<kiss_fft_cpx> fOut, int fStride, ReadOnlySpan<kiss_fft_cpx> twiddles, int m)
        {
            kiss_fft_cpx temp0, temp1;

            int tw1 = 0;
            int k = m;
            do
            {
                temp1 = fOut[m] * twiddles[tw1];
                tw1 += fStride;

                temp0 = fOut[0];

                fOut[0] = temp0 + temp1;

                fOut[m] = temp0 - temp1;

                fOut = fOut.Slice(1);
            } while (--k != 0);
        }

        private static void kf_bfly4(Span<kiss_fft_cpx> fOut, int fStride, ReadOnlySpan<kiss_fft_cpx> twiddles, int m)
        {
            kiss_fft_cpx temp0, temp1, temp2, temp3, tempA, tempB, tempC, tempD;

            int tw1 = 0;
            int k = m;
            do
            {
                temp1 = fOut[m] * twiddles[tw1];
                temp2 = fOut[2 * m] * twiddles[2 * tw1];
                temp3 = fOut[3 * m] * twiddles[3 * tw1];
                tw1 += fStride;

                temp0 = fOut[0];

                tempA = temp0 + temp2;
                tempB = temp0 - temp2;
                tempC = temp1 + temp3;
                tempD = temp1 - temp3;

                fOut[0] = tempA + tempC;

                fOut[2 * m] = tempA - tempC;

                fOut[m].r = tempB.r + tempD.i;
                fOut[m].i = tempB.i - tempD.r;

                fOut[3 * m].r = tempB.r - tempD.i;
                fOut[3 * m].i = tempB.i + tempD.r;

                fOut = fOut.Slice(1);
            } while (--k != 0);
        }

        private static void kf_bfly3(Span<kiss_fft_cpx> fOut, int fStride, ReadOnlySpan<kiss_fft_cpx> twiddles, int m)
        {
            kiss_fft_cpx temp0, temp1, temp2, tempA, tempB, tempC;

            T epi3 = twiddles[fStride * m].i;

            int tw1 = 0;
            int k = m;
            do
            {
                temp1 = fOut[m] * twiddles[tw1];
                temp2 = fOut[2 * m] * twiddles[tw1 * 2];
                tw1 += fStride;

                tempA = temp1 + temp2;
                tempB = temp1 - temp2;

                tempB *= epi3;

                temp0 = fOut[0];

                tempC.r = temp0.r - ToT(0.5) * tempA.r;
                tempC.i = temp0.i - ToT(0.5) * tempA.i;

                fOut[0] = temp0 + tempA;

                fOut[m].r = tempC.r - tempB.i;
                fOut[m].i = tempC.i + tempB.r;

                fOut[2 * m].r = tempC.r + tempB.i;
                fOut[2 * m].i = tempC.i - tempB.r;

                fOut = fOut.Slice(1);
            } while (--k != 0);
        }

        private static void kf_bfly5(Span<kiss_fft_cpx> fOut, int fStride, ReadOnlySpan<kiss_fft_cpx> twiddles, int m)
        {
            kiss_fft_cpx temp0, temp1, temp2, temp3, temp4, tempA, tempB, tempC,
                         tempD, tempE, tempF, tempG, tempH;

            var ya = twiddles[fStride * m];
            var yb = twiddles[2 * fStride * m];

            int tw1 = 0;
            int k = m;
            do
            {
                temp1 = fOut[m] * twiddles[tw1];
                temp2 = fOut[2 * m] * twiddles[2 * tw1];
                temp3 = fOut[3 * m] * twiddles[3 * tw1];
                temp4 = fOut[4 * m] * twiddles[4 * tw1];
                tw1 += fStride;

                tempA = temp1 + temp4;
                tempB = temp1 - temp4;
                tempC = temp2 + temp3;
                tempD = temp2 - temp3;

                temp0 = fOut[0];

                tempE = temp0 + tempA * ya.r + tempC * yb.r;
                tempF = temp0 + tempA * yb.r + tempC * ya.r;

                fOut[0] = temp0 + tempA + tempC;

                tempG.r = tempB.i * ya.i + tempD.i * yb.i;
                tempG.i = -tempB.r * ya.i - tempD.r * yb.i;

                fOut[m] = tempE - tempG;
                fOut[4 * m] = tempE + tempG;

                tempH.r = -tempB.i * yb.i + tempD.i * ya.i;
                tempH.i = tempB.r * yb.i - tempD.r * ya.i;

                fOut[2 * m] = tempF + tempH;
                fOut[3 * m] = tempF - tempH;

                fOut = fOut.Slice(1);
            } while (--k != 0);
        }

        private static void kf_bfly_generic(Span<kiss_fft_cpx> fOut, int fStride, ReadOnlySpan<kiss_fft_cpx> twiddles, int m, int p, Span<kiss_fft_cpx> scratch)
        {
            scratch = scratch.Slice(0, p);

            int k;

            for (int u = 0; u < m; ++u)
            {
                k = u;
                for (int q1 = 0; q1 < p; ++q1)
                {
                    scratch[q1] = fOut[k];

                    k += m;
                }

                k = u;
                for (int q1 = 0; q1 < p; ++q1)
                {
                    int twidx = 0;
                    var temp = scratch[0];
                    for (int q = 1; q < p; ++q)
                    {
                        twidx += fStride * k;
                        if (twidx >= twiddles.Length)
                            twidx -= twiddles.Length;
                        temp += scratch[q] * twiddles[twidx];
                    }
                    fOut[k] = temp;

                    k += m;
                }
            }
        }

        private static void kf_work(
            Span<kiss_fft_cpx> fOut,
            ReadOnlySpan<kiss_fft_cpx> fIn,
            int fStride,
            int fInStride,
            int factorsIndex,
            ref readonly kiss_fft_state state)
        {
            int p = state.factors[factorsIndex++]; /* the radix  */
            int m = state.factors[factorsIndex++]; /* stage's fft length/p */

            int fInIndex = 0;
            var fOutSlice = fOut.Slice(0, p * m);

            if (m == 1)
            {
                for (int fOutIndex = 0; fOutIndex < fOutSlice.Length;)
                {
                    fOut[fOutIndex] = fIn[fInIndex];
                    fInIndex += fStride * fInStride;
                    fOutIndex += 1;
                }
            }
            else
            {
                for (int fOutIndex = 0; fOutIndex < fOutSlice.Length;)
                {
                    // recursive call:
                    // DFT of size m*p performed by doing
                    // p instances of smaller DFTs of size m,
                    // each one takes a decimated version of the input
                    kf_work(fOut.Slice(fOutIndex), fIn.Slice(fInIndex), fStride * p, fInStride, factorsIndex, in state);
                    fInIndex += fStride * fInStride;
                    fOutIndex += m;
                }
            }

            // recombine the p smaller DFTs
            switch (p)
            {
                case 2:
                    kf_bfly2(fOut, fStride, state.twiddles, m);
                    break;
                case 3:
                    kf_bfly3(fOut, fStride, state.twiddles, m);
                    break;
                case 4:
                    kf_bfly4(fOut, fStride, state.twiddles, m);
                    break;
                case 5:
                    kf_bfly5(fOut, fStride, state.twiddles, m);
                    break;
                default:
                    kf_bfly_generic(fOut, fStride, state.twiddles, m, p, state.scratch);
                    break;
            }
        }

        private static void kiss_fft_stride(ref readonly kiss_fft_state state, ReadOnlySpan<kiss_fft_cpx> fIn, Span<kiss_fft_cpx> fOut, int fInStride)
        {
            Debug.Assert(fIn != fOut);

            kf_work(fOut, fIn, 1, fInStride, 0, in state);
        }

        internal static void kiss_fft(ref readonly kiss_fft_state state, ReadOnlySpan<kiss_fft_cpx> fIn, Span<kiss_fft_cpx> fOut)
        {
            kiss_fft_stride(in state, fIn, fOut, 1);
        }

#pragma warning restore IDE1006 // Naming Styles
    }
}
