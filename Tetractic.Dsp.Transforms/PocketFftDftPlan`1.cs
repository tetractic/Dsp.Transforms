// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using PocketFft;
using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Tetractic.Dsp.Transforms;

internal sealed class PocketFftDftPlan<T> : DftPlan<T>
    where T : IFloatingPointIeee754<T>
{
    private readonly PocketFft<T>.cfft_plan _plan;
    private readonly PocketFft<T>.cmplx[]? _akf;
    private readonly PocketFft<T>.cmplx[] _h;

    internal PocketFftDftPlan(int length)
        : base(length, canTransformInPlace: true)
    {
        Debug.Assert(length > 0);

        _plan = PocketFft<T>.make_cfft_plan(length);
        if (_plan.blueplan.n2 > 0)
        {
            _akf = new PocketFft<T>.cmplx[_plan.blueplan.n2];
            _h = new PocketFft<T>.cmplx[_plan.blueplan.n2];
        }
        else
        {
            _h = new PocketFft<T>.cmplx[_plan.packplan.length];
        }
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

        var cmplxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, PocketFft<T>.cmplx>(ref MemoryMarshal.GetReference(output)), Length);

        PocketFft<T>.cfft_forward(_plan, cmplxOutput, _akf, _h);
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

        var cmplxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, PocketFft<T>.cmplx>(ref MemoryMarshal.GetReference(output)), Length);

        PocketFft<T>.cfft_backward(_plan, cmplxOutput, _akf, _h);

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

        var cmplxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, PocketFft<T>.cmplx>(ref MemoryMarshal.GetReference(output)), Length);

        PocketFft<T>.cfft_backward(_plan, cmplxOutput, _akf, _h);
    }
}
