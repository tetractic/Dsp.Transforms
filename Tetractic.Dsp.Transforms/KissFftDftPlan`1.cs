// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using KissFft;
using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Tetractic.Dsp.Transforms;

internal sealed class KissFftDftPlan<T> : DftPlan<T>
    where T : notnull, IFloatingPointIeee754<T>
{
    private readonly KissFft<T>.kiss_fft_state _state;

    internal KissFftDftPlan(int length)
        : base(length, canTransformInPlace: false)
    {
        Debug.Assert(length > 0);

        var cfg = new KissFft<T>.kiss_fft_cfg(length);
        _state = new KissFft<T>.kiss_fft_state(in cfg);
    }

    /// <inheritdoc/>
    public override void Transform(ReadOnlySpan<T> input, Span<T> output)
    {
        if (input.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(input));
        if (output.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(output));

        if (input.Overlaps(output))
            throw new ArgumentException("The output memory overlaps the input memory.");

        var cpxInput = MemoryMarshal.CreateReadOnlySpan(in Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(input)), Length);
        var cpxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(output)), Length);

        KissFft<T>.Fft.kiss_fft(in _state, cpxInput, cpxOutput);
    }

    /// <inheritdoc/>
    public override void InverseTransform(ReadOnlySpan<T> input, Span<T> output)
    {
        if (input.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(input));
        if (output.Length != 2 * Length)
            throw new ArgumentException("The length is invalid.", nameof(output));

        if (input.Overlaps(output))
            throw new ArgumentException("The output memory overlaps the input memory.");

        var cpxInput = MemoryMarshal.CreateReadOnlySpan(in Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(input)), Length);
        var cpxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(output)), Length);

        KissFft<T>.Ifft.kiss_fft(in _state, cpxInput, cpxOutput);

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

        if (input.Overlaps(output))
            throw new ArgumentException("The output memory overlaps the input memory.");

        var cpxInput = MemoryMarshal.CreateReadOnlySpan(in Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(input)), Length);
        var cpxOutput = MemoryMarshal.CreateSpan(ref Unsafe.As<T, KissFft<T>.kiss_fft_cpx>(ref MemoryMarshal.GetReference(output)), Length);

        KissFft<T>.Ifft.kiss_fft(in _state, cpxInput, cpxOutput);
    }
}
