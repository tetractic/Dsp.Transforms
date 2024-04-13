// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause

using BenchmarkDotNet.Running;

namespace Tetractic.Dsp.Transforms.Benchmarks;

public static class Program
{
    public static void Main()
    {
        _ = BenchmarkRunner.Run<OouraDftPlanBenchmarks>();
        _ = BenchmarkRunner.Run<PocketFftDftPlanBenchmarks>();
    }
}
