// Copyright (C) 2004-2018 Max-Planck-Society
// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mreineck/pocketfft

using System;
using System.Runtime.CompilerServices;
using static PocketFft.PocketFft;

namespace PocketFft;

internal static partial class PocketFft<T>
{
#pragma warning disable IDE1006 // Naming Styles

    internal struct fftblue_plan
    {
        internal int n;
        internal int n2;
        internal cfftp_plan plan;
        internal cmplx[] bk;
        internal cmplx[] bkf;

        internal static fftblue_plan make_fftblue_plan(int length)
        {
            var plan = new fftblue_plan();
            plan.n = length;
            plan.n2 = good_size(plan.n * 2 - 1);
            plan.bk = new cmplx[plan.n];
            plan.bkf = new cmplx[plan.n2];

            var bk = plan.bk.AsSpan();
            var bkf = plan.bkf.AsSpan();

            // Initialize b_k
            var tmp = new T[4 * plan.n];
            sincos_2pibyn(2 * plan.n, tmp);
            bk[0].r = T.One;
            bk[0].i = T.Zero;

            int coeff = 0;
            for (int m = 1; m < plan.n; ++m)
            {
                coeff += 2 * m - 1;
                if (coeff >= 2 * plan.n)
                    coeff -= 2 * plan.n;
                bk[m].r = tmp[2 * coeff];
                bk[m].i = tmp[2 * coeff + 1];
            }

            // Initialize the zero-padded, normalized, Fourier-transformed b_k
            T xn2 = T.One / ToT(plan.n2);
            bkf[0].r = bk[0].r * xn2;
            bkf[0].i = bk[0].i * xn2;
            for (int m = 1; m < plan.n; ++m)
            {
                bkf[m].r = bkf[plan.n2 - m].r = bk[m].r * xn2;
                bkf[m].i = bkf[plan.n2 - m].i = bk[m].i * xn2;
            }
            for (int m = plan.n; m <= plan.n2 - plan.n; ++m)
            {
                bkf[m].r = T.Zero;
                bkf[m].i = T.Zero;
            }
            plan.plan = cfftp_plan.make_cfftp_plan(plan.n2);
            var ch = new cmplx[plan.n2];
            cfftp_plan.cfftp_forward(in plan.plan, bkf, ch);

            return plan;
        }

        internal static void cfftblue_backward(ref readonly fftblue_plan plan, Span<cmplx> c, Span<cmplx> akf, Span<cmplx> akfh)
        {
            int n = plan.n;
            int n2 = plan.n2;
            var bk = plan.bk.AsSpan();
            var bkf = plan.bkf.AsSpan();

            // Initialize a_k
            for (int m = 0; m < n; ++m)
                MultiplyAB(out akf[m], c[m], bk[m]);
            for (int m = n; m < n2; ++m)
            {
                akf[m].r = T.Zero;
                akf[m].i = T.Zero;
            }

            // FFT
            cfftp_plan.cfftp_forward(in plan.plan, akf, akfh);

            // Convolution
            for (int m = 0; m < n2; ++m)
                MultiplyAConjB(out akf[m], akf[m], bkf[m]);

            // Inverse FFT
            cfftp_plan.cfftp_backward(in plan.plan, akf, akfh);

            // Multiply by b_k
            for (int m = 0; m < n; ++m)
                MultiplyAB(out c[m], bk[m], akf[m]);
        }

        internal static void cfftblue_forward(ref readonly fftblue_plan plan, Span<cmplx> c, Span<cmplx> akf, Span<cmplx> akfh)
        {
            int n = plan.n;
            int n2 = plan.n2;
            var bk = plan.bk.AsSpan();
            var bkf = plan.bkf.AsSpan();

            // Initialize a_k
            for (int m = 0; m < n; ++m)
                MultiplyAConjB(out akf[m], c[m], bk[m]);
            for (int m = n; m < n2; ++m)
            {
                akf[m].r = T.Zero;
                akf[m].i = T.Zero;
            }

            // FFT
            cfftp_plan.cfftp_forward(in plan.plan, akf, akfh);

            // Convolution
            for (int m = 0; m < n2; ++m)
                MultiplyAB(out akf[m], akf[m], bkf[m]);

            // Inverse FFT
            cfftp_plan.cfftp_backward(in plan.plan, akf, akfh);

            // Multiply by b_k
            for (int m = 0; m < n; ++m)
                MultiplyConjAB(out c[m], bk[m], akf[m]);
        }

        // d = a * b
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void MultiplyAB(out cmplx d, cmplx a, cmplx b)
        {
            d.r = T.FusedMultiplyAdd(a.r, b.r, -a.i * b.i);
            d.i = T.FusedMultiplyAdd(a.r, b.i, a.i * b.r);
        }

        // d = conj(a) * b
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void MultiplyConjAB(out cmplx d, cmplx a, cmplx b)
        {
            d.r = T.FusedMultiplyAdd(a.r, b.r, a.i * b.i);
            d.i = T.FusedMultiplyAdd(a.r, b.i, -a.i * b.r);
        }

        // d = a * conj(b)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void MultiplyAConjB(out cmplx d, cmplx a, cmplx b)
        {
            d.r = T.FusedMultiplyAdd(a.r, b.r, a.i * b.i);
            d.i = T.FusedMultiplyAdd(a.i, b.r, -a.r * b.i);
        }
    }

#pragma warning restore IDE1006 // Naming Styles
}
