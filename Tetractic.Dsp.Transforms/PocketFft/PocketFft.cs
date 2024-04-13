// Copyright (C) 2004-2018 Max-Planck-Society
// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mreineck/pocketfft

using System.Runtime.CompilerServices;

namespace PocketFft;

internal static class PocketFft
{
#pragma warning disable IDE1006 // Naming Styles

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static void SWAP<T>(ref T left, ref T right)
    {
        T temp = left;
        left = right;
        right = temp;
    }

    internal static int largest_prime_factor(int n)
    {
        int res = 1;
        int tmp;

        while (((tmp = (n >> 1)) << 1) == n)
        {
            res = 2;
            n = tmp;
        }

        int limit = (int)double.Sqrt(n);
        for (int x = 3; x <= limit; x += 2)
        {
            while (((tmp = (n / x)) * x) == n)
            {
                res = x;
                n = tmp;
                limit = (int)double.Sqrt(n);
            }
        }
        if (n > 1)
            res = n;

        return res;
    }

    internal static double cost_guess(int n)
    {
        const double lfp = 1.1; // penalty for non-hardcoded larger factors
        int ni = n;
        double result = 0;
        int tmp;

        while (((tmp = (n >> 1)) << 1) == n)
        {
            result += 2;
            n = tmp;
        }

        int limit = (int)double.Sqrt(n);
        for (int x = 3; x <= limit; x += 2)
        {
            while ((tmp = (n / x)) * x == n)
            {
                result += (x <= 5) ? x : lfp * x; // penalize larger prime factors
                n = tmp;
                limit = (int)double.Sqrt(n);
            }
        }
        if (n > 1)
            result += (n <= 5) ? n : lfp * n;

        return result * ni;
    }

    internal static int good_size(int n)
    {
        if (n <= 6)
            return n;

        int bestfac = 2 * n;
        for (int f2 = 1; f2 < bestfac; f2 *= 2)
            for (int f23 = f2; f23 < bestfac; f23 *= 3)
                for (int f235 = f23; f235 < bestfac; f235 *= 5)
                    for (int f2357 = f235; f2357 < bestfac; f2357 *= 7)
                        for (int f235711 = f2357; f235711 < bestfac; f235711 *= 11)
                            if (f235711 >= n)
                                bestfac = f235711;
        return bestfac;
    }

#pragma warning restore IDE1006 // Naming Styles
}
