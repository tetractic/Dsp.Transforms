// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using System;
using System.Numerics;

namespace Tetractic.Dsp.Transforms;

/// <summary>
/// Provides methods for creating a discrete Fourier transform plan.
/// </summary>
public static class DftPlan
{
    /// <summary>
    /// Creates a discrete Fourier transform plan.
    /// </summary>
    /// <typeparam name="T">The type of values to transform.</typeparam>
    /// <param name="length">The number of complex values to transform.</param>
    /// <returns>The plan.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="length"/> is less than 1 or is
    ///     greater than half of <see cref="int.MaxValue"/>.</exception>
    public static DftPlan<T> Create<T>(int length)
        where T : notnull, IFloatingPointIeee754<T>
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(length, int.MaxValue / 2);

        return int.IsPow2(length)
            ? new OouraDftPlan<T>(length)
            : new PocketFftDftPlan<T>(length);
    }
}
