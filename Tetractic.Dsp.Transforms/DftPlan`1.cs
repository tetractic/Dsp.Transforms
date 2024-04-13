// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using System;
using System.Numerics;

namespace Tetractic.Dsp.Transforms;

/// <summary>
/// Represents a plan to perform a discrete Fourier transform or its inverse.
/// </summary>
/// <typeparam name="T">The type of values to transform.</typeparam>
public abstract class DftPlan<T>
    where T : notnull, IFloatingPointIeee754<T>
{
    /// <summary>
    /// Initializes a new instance of <see cref="DftPlan{T}"/>.
    /// </summary>
    /// <param name="length">The number of complex values that are transformed.</param>
    /// <param name="canTransformInPlace"><see langword="true"/> if the input and output of the
    ///     transform can be the same memory; otherwise, <see langword="false"/>.</param>
    protected DftPlan(int length, bool canTransformInPlace)
    {
        Length = length;
        CanTransformInPlace = canTransformInPlace;
    }

    /// <summary>
    /// Gets the number of complex values that are transformed.
    /// </summary>
    public int Length { get; }

    /// <summary>
    /// Gets a value that indicates whether the input and output of the transform can be the same
    /// memory.
    /// </summary>
    public bool CanTransformInPlace { get; }

    /// <summary>
    /// Computes the DFT of complex values.
    /// </summary>
    /// <param name="input">The interleaved real and imaginary components of the complex values to
    ///     transform.</param>
    /// <param name="output">Returns the transformed values.  If <see cref="CanTransformInPlace"/>
    ///     is <see langword="true"/> then <paramref name="output"/> can be the same memory as
    ///     <paramref name="input"/>.</param>
    /// <exception cref="ArgumentException">The length of <paramref name="input"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <exception cref="ArgumentException">The length of <paramref name="output"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <remarks>
    /// <para>
    /// Computes <c>X[k] = sum from n=0 to n=N-1 of x[n]*exp(-2*π*i*n*k/N)</c>.
    /// </para>
    /// <para>
    /// Input:
    /// <code>
    /// input[2*k] = Re(x[k])
    /// input[2*k+1] = Im(x[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// <para>
    /// Output:
    /// <code>
    /// output[2*k] = Re(X[k])
    /// output[2*k+1] = Im(X[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// </remarks>
    public abstract void Transform(ReadOnlySpan<T> input, Span<T> output);

    /// <summary>
    /// Computes the inverse DFT of complex values.
    /// </summary>
    /// <param name="input">The interleaved real and imaginary components of the complex values to
    ///     transform.</param>
    /// <param name="output">Returns the normalized transformed values.  If
    ///     <see cref="CanTransformInPlace"/> is <see langword="true"/> then
    ///     <paramref name="output"/> can be the same memory as <paramref name="input"/>.</param>
    /// <exception cref="ArgumentException">The length of <paramref name="input"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <exception cref="ArgumentException">The length of <paramref name="output"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <remarks>
    /// <para>
    /// Computes <c>x[k] = (1/N) * sum from n=0 to n=N-1 of X[n]*exp(2*π*i*n*k/N)</c>.
    /// </para>
    /// <para>
    /// Input:
    /// <code>
    /// input[2*k] = Re(X[k])
    /// input[2*k+1] = Im(X[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// <para>
    /// Output:
    /// <code>
    /// output[2*k] = Re(x[k])
    /// output[2*k+1] = Im(x[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// </remarks>
    public abstract void InverseTransform(ReadOnlySpan<T> input, Span<T> output);

    /// <summary>
    /// Computes the inverse DFT of complex values.
    /// </summary>
    /// <param name="input">The interleaved real and imaginary components of the complex values to
    ///     transform.</param>
    /// <param name="output">Returns the unnormalized transformed values.  If
    ///     <see cref="CanTransformInPlace"/> is <see langword="true"/> then
    ///     <paramref name="output"/> can be the same memory as <paramref name="input"/>.</param>
    /// <exception cref="ArgumentException">The length of <paramref name="input"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <exception cref="ArgumentException">The length of <paramref name="output"/> is not twice
    ///     <see cref="Length"/>.</exception>
    /// <remarks>
    /// <para>
    /// Computes <c>x[k] = sum from n=0 to n=N-1 of X[n]*exp(2*π*i*n*k/N)</c>.
    /// </para>
    /// <para>
    /// Input:
    /// <code>
    /// input[2*k] = Re(X[k])
    /// input[2*k+1] = Im(X[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// <para>
    /// Output:
    /// <code>
    /// output[2*k] = Re(x[k])
    /// output[2*k+1] = Im(x[k])
    /// where 0 ≤ k &lt; Length
    /// </code>
    /// </para>
    /// </remarks>
    public abstract void InverseTransformUnnormalized(ReadOnlySpan<T> input, Span<T> output);
}
