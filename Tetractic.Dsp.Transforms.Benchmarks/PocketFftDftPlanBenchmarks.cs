// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using BenchmarkDotNet.Attributes;

namespace Tetractic.Dsp.Transforms.Benchmarks;

public class PocketFftDftPlanBenchmarks
{
    private PocketFftDftPlan<double> _plan = null!;

    private double[] _data = null!;

    [Params(
    [
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
    ])]
    public int Length;

    [GlobalSetup]
    public void Setup()
    {
        _plan = new(Length);
        _data = GetRandom(Length);

        static double[] GetRandom(int n)
        {
            double[] result = new double[n * 2];
            for (int i = 0; i < result.Length; ++i)
                result[i] = Random.Shared.NextDouble();
            return result;
        }
    }

    [Benchmark]
    public void Transform()
    {
        _plan.Transform(_data, _data);
    }
}
