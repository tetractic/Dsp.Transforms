// Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
// Copyright (c) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mborgerding/kissfft

using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace KissFft;

internal static partial class KissFft<T>
    where T : notnull, IFloatingPointIeee754<T>
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static T ToT(int value) => T.CreateTruncating(value);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static T ToT(double value) => T.CreateTruncating(value);

#pragma warning disable IDE1006 // Naming Styles

    private const int MAXFACTORS = 32;

    [StructLayout(LayoutKind.Sequential)]
    internal struct kiss_fft_cpx
    {
        public T r;
        public T i;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static kiss_fft_cpx operator *(kiss_fft_cpx a, kiss_fft_cpx b)
        {
            return new kiss_fft_cpx
            {
                r = a.r * b.r - a.i * b.i,
                i = a.r * b.i + a.i * b.r,
            };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static kiss_fft_cpx operator *(kiss_fft_cpx c, T s)
        {
            return new kiss_fft_cpx
            {
                r = c.r * s,
                i = c.i * s,
            };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static kiss_fft_cpx operator +(kiss_fft_cpx a, kiss_fft_cpx b)
        {
            return new kiss_fft_cpx
            {
                r = a.r + b.r,
                i = a.i + b.i,
            };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static kiss_fft_cpx operator -(kiss_fft_cpx a, kiss_fft_cpx b)
        {
            return new kiss_fft_cpx
            {
                r = a.r - b.r,
                i = a.i - b.i,
            };
        }
    }

    internal readonly struct kiss_fft_cfg
    {
        internal readonly Factors factors;
        internal readonly kiss_fft_cpx[] twiddles;

        public kiss_fft_cfg(int nfft)
        {
            twiddles = new kiss_fft_cpx[nfft];

            for (int i = 0; i < nfft; ++i)
            {
                T phase = ToT(-2) * ToT(i) / ToT(nfft);

                var (sinPi, cosPi) = T.SinCosPi(phase);
                twiddles[i] = new kiss_fft_cpx
                {
                    r = cosPi,
                    i = sinPi,
                };
            }

            kf_factor(nfft, ref factors);
        }

        private static void kf_factor(int n, ref Factors factors)
        {
            int maxFactor = (int)double.Floor(double.Sqrt(n));

            /*factor out powers of 4, powers of 2, then any remaining primes */
            int p = 4;
            int factorsIndex = 0;
            do
            {
                while (n % p != 0)
                {
                    switch (p)
                    {
                        case 4:
                            p = 2;
                            break;
                        case 2:
                            p = 3;
                            break;
                        default:
                            p += 2;
                            break;
                    }
                    if (p > maxFactor)
                        p = n;          /* no more factors, skip to end */
                }
                n /= p;
                factors[factorsIndex++] = p;
                factors[factorsIndex++] = n;
            } while (n > 1);
        }

        [InlineArray(2 * MAXFACTORS)]
        public struct Factors
        {
#pragma warning disable IDE0044 // Add readonly modifier
#pragma warning disable IDE0051 // Remove unused private members
            private int _element0;
#pragma warning restore IDE0044 // Add readonly modifier
#pragma warning restore IDE0051 // Remove unused private members
        }
    }

    internal readonly struct kiss_fft_state
    {
        internal readonly kiss_fft_cfg.Factors factors;
        internal readonly kiss_fft_cpx[] twiddles;
        internal readonly kiss_fft_cpx[] scratch;

        public kiss_fft_state(ref readonly kiss_fft_cfg cfg)
        {
            factors = cfg.factors;
            twiddles = cfg.twiddles;

            int maxP = 0;
            for (int i = 0; i < MAXFACTORS; ++i)
                if (maxP < cfg.factors[i * 2])
                    maxP = cfg.factors[i * 2];
            scratch = maxP > 5
                ? new kiss_fft_cpx[maxP]
                : null!;
        }
    }

#pragma warning restore IDE1006 // Naming Styles
}
