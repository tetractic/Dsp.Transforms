// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using System;
using System.Diagnostics;
using System.IO;
using Xunit;
using Xunit.Abstractions;

namespace Tetractic.Dsp.Transforms.Tests;

public class DftPlanTests
{
    private const string _fftwInputFormatString = "fftw.{0}.in";

    private const string _fftwOutputFormatString = "fftw.{0}.out";

    private static readonly double _epsilon = double.Pow(10, -15.55);

    public static readonly TheoryData<int> Pow2Lengths = new()
    {
        1,       // 2^0
        2,       // 2^1
        4,       // 2^2
        8,       // 2^3
        16,      // 2^4
        32,      // 2^5
        64,      // 2^6
        128,     // 2^7
        256,     // 2^8
        512,     // 2^9
        1024,    // 2^10
        2048,    // 2^11
        4096,    // 2^12
        8192,    // 2^13
        16384,   // 2^14
        32768,   // 2^15
        65536,   // 2^16
        131072,  // 2^17
    };

    public static readonly TheoryData<int> Lengths = new()
    {
        3,       // 3
        5,       // 5
        6,       // 2 × 3
        7,       // 7
        9,       // 3^2
        10,      // 2 × 5
        11,      // 11
        12,      // 3 × 4
        13,      // 13
        14,      // 2 × 7
        15,      // 3 × 5
        20,      // 4 × 5
        21,      // 3 × 7
        22,      // 2 × 11
        24,      // 2 × 3 × 4
        26,      // 2 × 13
        28,      // 4 × 7
        33,      // 3 × 11
        35,      // 5 × 7
        39,      // 3 × 13
        44,      // 4 × 11
        49,      // 7^2
        52,      // 4 × 13
        55,      // 5 × 11
        65,      // 5 × 13
        77,      // 7 × 11
        91,      // 7 × 13
        97,      // 97
        108,     // 2^2 × 3^3
        120,     // 2 × 3 × 4 × 5
        121,     // 11^2
        143,     // 11 × 13
        169,     // 13^2
        194,     // 2 × 97
        840,     // 2 × 3 × 4 × 5 × 7
        997,     // 997
        9240,    // 2 × 3 × 4 × 5 × 7 × 11
        9973,    // 9973
        120120,  // 2 × 3 × 4 × 5 × 7 × 11 × 13
    };

    private readonly ITestOutputHelper _output;

    public DftPlanTests(ITestOutputHelper output)
    {
        _output = output;
    }

    [Theory]
    [MemberData(nameof(Pow2Lengths))]
    [MemberData(nameof(Lengths))]
    public void Transform(int n)
    {
        var fft = DftPlan.Create<double>(n);
        double[] actual = ReadValues(string.Format(_fftwInputFormatString, n), n * 2);
        double[] expected = ReadValues(string.Format(_fftwOutputFormatString, n), n * 2);
        if (fft.CanTransformInPlace)
        {
            fft.Transform(actual, actual);
        }
        else
        {
            throw new NotImplementedException();
        }
        double error = Compare(actual, expected);
        double tolerance = _epsilon * double.Sqrt(Math.Log(n));
        _output.WriteLine($"Error = {error}; Tolerance = {tolerance}");
        Assert.Equal(0, error, tolerance);
    }

    [Theory]
    [MemberData(nameof(Pow2Lengths))]
    [MemberData(nameof(Lengths))]
    public void InverseTransform(int n)
    {
        var fft = DftPlan.Create<double>(n);
        double[] actual = ReadValues(string.Format(_fftwOutputFormatString, n), n * 2);
        double[] expected = ReadValues(string.Format(_fftwInputFormatString, n), n * 2);
        if (fft.CanTransformInPlace)
        {
            fft.InverseTransform(actual, actual);
        }
        else
        {
            throw new NotImplementedException();
        }
        double error = Compare(actual, expected);
        double tolerance = _epsilon * double.Sqrt(Math.Log(n));
        _output.WriteLine($"Error = {error}; Tolerance = {tolerance}");
        Assert.Equal(0, error, tolerance);
    }

    private static double[] ReadValues(string filename, int length)
    {
        string directoryPath = Path.GetDirectoryName(typeof(DftPlanTests).Assembly.Location)!;
        string path = Path.Combine(directoryPath, "Resources", filename);

        using (var stream = new FileStream(path, FileMode.Open, FileAccess.Read))
        using (var reader = new BinaryReader(stream))
        {
            double[] result = new double[length];

            for (int i = 0; i < result.Length; ++i)
                result[i] = reader.ReadDouble();

            if (stream.Position != stream.Length)
                throw new UnreachableException();

            return result;
        }
    }

    // http://www.fftw.org/accuracy/method.html
    public static double Compare(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Lengths are mismatched.");

        double sum1 = 0;
        double sum2 = 0;
        for (int i = 0; i < a.Length; ++i)
        {
            double x = a[i] - b[i];
            sum1 += x * x;
            sum2 += b[i] * b[i];
        }
        double num = double.Sqrt(sum1);
        double den = double.Sqrt(sum2);
        return num / den;
    }
}
