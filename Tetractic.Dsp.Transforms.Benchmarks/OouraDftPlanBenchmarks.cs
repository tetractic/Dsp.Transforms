// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using BenchmarkDotNet.Attributes;

namespace Tetractic.Dsp.Transforms.Benchmarks;

public class OouraDftPlanBenchmarks
{
    private OouraDftPlan<double> _plan = null!;
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
