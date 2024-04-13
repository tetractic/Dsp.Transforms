// Copyright (C) 2004-2018 Max-Planck-Society
// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mreineck/pocketfft

using System;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.CompilerServices;
using static PocketFft.PocketFft;

namespace PocketFft;

internal static partial class PocketFft<T>
    where T : notnull, IFloatingPointIeee754<T>
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static T ToT(int value) => T.CreateTruncating(value);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static T FMA(T left, T right, T addend) => T.FusedMultiplyAdd(left, right, addend);

#pragma warning disable IDE1006 // Naming Styles

    private static void calc_first_octant(int den, Span<T> res)
    {
        int n = (den + 4) >> 3;
        if (n == 0)
            return;
        res[0] = T.One;
        res[1] = T.Zero;
        if (n == 1)
            return;
        for (int i = 1; i < n; ++i)
        {
            var (sinPi, cosPi) = T.SinCosPi(ToT(2) * ToT(i) / ToT(den));
            res[2 * i] = cosPi;
            res[2 * i + 1] = sinPi;
        }
    }

    private static void calc_first_quadrant(int n, Span<T> res)
    {
        Span<T> p = res.Slice(n);
        calc_first_octant(n << 1, p);
        int ndone = (n + 2) >> 2;
        int i = 0, idx1 = 0, idx2 = 2 * ndone - 2;
        for (; i + 1 < ndone; i += 2, idx1 += 2, idx2 -= 2)
        {
            res[idx1] = p[2 * i];
            res[idx1 + 1] = p[2 * i + 1];
            res[idx2] = p[2 * i + 3];
            res[idx2 + 1] = p[2 * i + 2];
        }
        if (i != ndone)
        {
            res[idx1] = p[2 * i];
            res[idx1 + 1] = p[2 * i + 1];
        }
    }

    private static void calc_first_half(int n, Span<T> res)
    {
        int ndone = (n + 1) >> 1;
        var p = res.Slice(n - 1);
        calc_first_octant(n << 2, p);
        int i4 = 0, @in = n, i = 0;
        for (; i4 <= @in - i4; ++i, i4 += 4) // octant 0
        {
            res[2 * i] = p[2 * i4];
            res[2 * i + 1] = p[2 * i4 + 1];
        }
        for (; i4 - @in <= 0; ++i, i4 += 4) // octant 1
        {
            int xm = @in - i4;
            res[2 * i] = p[2 * xm + 1];
            res[2 * i + 1] = p[2 * xm];
        }
        for (; i4 <= 3 * @in - i4; ++i, i4 += 4) // octant 2
        {
            int xm = i4 - @in;
            res[2 * i] = -p[2 * xm + 1];
            res[2 * i + 1] = p[2 * xm];
        }
        for (; i < ndone; ++i, i4 += 4) // octant 3
        {
            int xm = 2 * @in - i4;
            res[2 * i] = -p[2 * xm];
            res[2 * i + 1] = p[2 * xm + 1];
        }
    }

    private static void fill_first_quadrant(int n, Span<T> res)
    {
        T two = T.One + T.One;
        T hsqt2 = T.Sqrt(two) / two;
        int quart = n >> 2;
        if ((n & 7) == 0)
            res[quart] = res[quart + 1] = hsqt2;
        for (int i = 2, j = 2 * quart - 2; i < quart; i += 2, j -= 2)
        {
            res[j] = res[i + 1];
            res[j + 1] = res[i];
        }
    }

    private static void fill_first_half(int n, Span<T> res)
    {
        int half = n >> 1;
        if ((n & 3) == 0)
        {
            for (int i = 0; i < half; i += 2)
            {
                res[i + half] = -res[i + 1];
                res[i + half + 1] = res[i];
            }
        }
        else
        {
            for (int i = 2, j = 2 * half - 2; i < half; i += 2, j -= 2)
            {
                res[j] = -res[i];
                res[j + 1] = res[i + 1];
            }
        }
    }

    private static void fill_second_half(int n, Span<T> res)
    {
        if ((n & 1) == 0)
        {
            for (int i = 0; i < n; ++i)
                res[i + n] = -res[i];
        }
        else
        {
            for (int i = 2, j = 2 * n - 2; i < n; i += 2, j -= 2)
            {
                res[j] = res[i];
                res[j + 1] = -res[i + 1];
            }
        }
    }

    private static void sincos_2pibyn_half(int n, Span<T> res)
    {
        if ((n & 3) == 0)
        {
            calc_first_octant(n, res);
            fill_first_quadrant(n, res);
            fill_first_half(n, res);
        }
        else if ((n & 1) == 0)
        {
            calc_first_quadrant(n, res);
            fill_first_half(n, res);
        }
        else
        {
            calc_first_half(n, res);
        }
    }

    private static void sincos_2pibyn(int n, Span<T> res)
    {
        sincos_2pibyn_half(n, res);
        fill_second_half(n, res);
    }

    /// <exception cref="ArgumentOutOfRangeException"/>
    internal static cfft_plan make_cfft_plan(int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(length);

        cfft_plan plan = default;

        if (length < 50 || Sqr(largest_prime_factor(length)) <= length)
        {
            plan.packplan = cfftp_plan.make_cfftp_plan(length);
            return plan;
        }

        double comp1 = cost_guess(length);
        double comp2 = 2 * cost_guess(good_size(2 * length - 1));
        comp2 *= 1.5; // fudge factor that appears to give good overall performance
        if (comp2 < comp1) // use Bluestein
            plan.blueplan = fftblue_plan.make_fftblue_plan(length);
        else
            plan.packplan = cfftp_plan.make_cfftp_plan(length);
        return plan;

        static double Sqr(double x) => x * x;
    }

    internal static void cfft_backward(cfft_plan plan, Span<cmplx> c, Span<cmplx> akf, Span<cmplx> h)
    {
        if (plan.blueplan.n == 0)
            cfftp_plan.cfftp_backward(in plan.packplan, c, h);
        else
            fftblue_plan.cfftblue_backward(in plan.blueplan, c, akf, h);
    }

    internal static void cfft_forward(cfft_plan plan, Span<cmplx> c, Span<cmplx> akf, Span<cmplx> h)
    {
        if (plan.blueplan.n == 0)
            cfftp_plan.cfftp_forward(in plan.packplan, c, h);
        else
            fftblue_plan.cfftblue_forward(in plan.blueplan, c, akf, h);
    }

    internal struct cfft_plan
    {
        [UnscopedRef]
        internal ref cfftp_plan packplan => ref blueplan.plan;
        internal fftblue_plan blueplan;
    }

    [DebuggerDisplay("({r}, {i})")]
#pragma warning disable CS8981 // The type name only contains lower-cased ascii characters. Such names may become reserved for the language.
    internal struct cmplx
#pragma warning restore CS8981 // The type name only contains lower-cased ascii characters. Such names may become reserved for the language.
    {
        public T r;
        public T i;
    }

#pragma warning restore IDE1006 // Naming Styles
}
