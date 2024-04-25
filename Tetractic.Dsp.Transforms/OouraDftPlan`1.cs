// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using Ooura;
using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Tetractic.Dsp.Transforms;

internal sealed class OouraDftPlan<T> : DftPlan<T>
    where T : notnull, IFloatingPointIeee754<T>
{
    private readonly int[] _ip;
    private readonly T[] _w;

    internal OouraDftPlan(int length)
        : base(length, canTransformInPlace: true)
    {
        Debug.Assert(length > 0 && int.IsPow2(length));

        _ip = new int[2 + (1 << (int)double.Log2(length + 0.5) / 2)];
        _w = new T[length / 2];

        int nw = (2 * length) >> 2;
        Fftsg<T>.makewt(nw, _ip, _w);
    }

    /// <inheritdoc/>
    public override void Transform(ReadOnlySpan<T> input, Span<T> output)
    {
        if (input.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(input));
        if (output.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(output));

        if (!Unsafe.AreSame(in MemoryMarshal.GetReference(input), in MemoryMarshal.GetReference(output)))
        {
            if (input.Overlaps(output))
                throw new ArgumentException("The output memory overlaps the input memory.");

            input.CopyTo(output);
        }

        Fftsg<T>.cdft(2 * Length, output, _ip, _w);
    }

    /// <inheritdoc/>
    public override void InverseTransform(ReadOnlySpan<T> input, Span<T> output)
    {
        if (input.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(input));
        if (output.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(output));

        if (!Unsafe.AreSame(in MemoryMarshal.GetReference(input), in MemoryMarshal.GetReference(output)))
        {
            if (input.Overlaps(output))
                throw new ArgumentException("The output memory overlaps the input memory.");

            input.CopyTo(output);
        }

        Fftsg<T>.icdft(2 * Length, output, _ip, _w);

        T scale = T.One / T.CreateTruncating(Length);
        for (int i = 0; i < output.Length; ++i)
            output[i] *= scale;
    }

    /// <inheritdoc/>
    public override void InverseTransformUnnormalized(ReadOnlySpan<T> input, Span<T> output)
    {
        if (input.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(input));
        if (output.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(output));

        if (!Unsafe.AreSame(in MemoryMarshal.GetReference(input), in MemoryMarshal.GetReference(output)))
        {
            if (input.Overlaps(output))
                throw new ArgumentException("The output memory overlaps the input memory.");

            input.CopyTo(output);
        }

        Fftsg<T>.icdft(2 * Length, output, _ip, _w);
    }
}
