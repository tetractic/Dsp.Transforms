// Copyright (C) 2004-2018 Max-Planck-Society
// Copyright (C) 2024 Carl Reinke
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/mreineck/pocketfft

using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using static PocketFft.PocketFft;

namespace PocketFft;

internal static partial class PocketFft<T>
{
#pragma warning disable IDE1006 // Naming Styles

    internal struct cfftp_fctdata
    {
        internal int fct;
        internal cmplx[] tw;
        internal cmplx[] tws;
    }

    internal struct cfftp_plan
    {
        private const int NFCT = 25;

        internal int length;
        internal int nfct;
        internal Fct fct;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void PMC(out cmplx a, out cmplx b, cmplx c, cmplx d)
        {
            a.r = c.r + d.r;
            a.i = c.i + d.i;
            b.r = c.r - d.r;
            b.i = c.i - d.i;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ADDC(out cmplx a, cmplx b, cmplx c)
        {
            a.r = b.r + c.r;
            a.i = b.i + c.i;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ROT90(ref cmplx a)
        {
            T ar = a.r;
            a.r = -a.i;
            a.i = ar;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ROTM90(ref cmplx a)
        {
            T nar = -a.r;
            a.r = a.i;
            a.i = nar;
        }

        // a = b * c
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void A_EQ_B_MUL_C(out cmplx a, cmplx b, cmplx c)
        {
            a.r = FMA(b.r, c.r, -b.i * c.i);
            a.i = FMA(b.r, c.i, b.i * c.r);
        }

        // a = conj(b) * c
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void A_EQ_CB_MUL_C(out cmplx a, cmplx b, cmplx c)
        {
            a.r = FMA(b.r, c.r, b.i * c.i);
            a.i = FMA(b.r, c.i, -b.i * c.r);
        }

        internal interface IDirection
        {
        }

        internal readonly struct F : IDirection
        {
        }

        internal readonly struct B : IDirection
        {
        }

        private readonly ref struct Pass2<D>
            where D : struct, IDirection
        {
            private const int cdim = 2;

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            internal Pass2(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            internal void pass2()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                        PMC(out CH(0, k, 0), out CH(0, k, 1), CC(0, 0, k), CC(0, 1, k));
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        PMC(out CH(0, k, 0), out CH(0, k, 1), CC(0, 0, k), CC(0, 1, k));
                        for (int i = 1; i < ido; ++i)
                        {
                            cmplx t;
                            PMC(out CH(i, k, 0), out t, CC(i, 0, k), CC(i, 1, k));
                            if (typeof(D) == typeof(B))
                                A_EQ_B_MUL_C(out CH(i, k, 1), WA(1 - 1, i), t);
                            else
                                A_EQ_CB_MUL_C(out CH(i, k, 1), WA(1 - 1, i), t);

                        }
                    }
                }
            }
        }

        private readonly ref struct Pass3<D>
            where D : struct, IDirection
        {
            private const int cdim = 3;

            private static readonly T tw1r;
            private static readonly T tw1i;

            static Pass3()
            {
                (tw1i, tw1r) = T.SinCosPi(ToT(2) / ToT(3));
                if (typeof(D) == typeof(F))
                    tw1i = -tw1i;
            }

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            public Pass3(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP3a(int k, int i)
            {
                cmplx t0 = CC(i, 0, k), t1, t2;
                PMC(out t1, out t2, CC(i, 1, k), CC(i, 2, k));
                CH(i, k, 0).r = t0.r + t1.r;
                CH(i, k, 0).i = t0.i + t1.i;

                const int u1 = 1;
                const int u2 = 2;
                cmplx ca, cb;
                ca.r = FMA(tw1r, t1.r, t0.r);
                ca.i = FMA(tw1r, t1.i, t0.i);
                cb.i = tw1i * t2.r;
                cb.r = -(tw1i * t2.i);
                PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);  // TODO: FMA?
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP3(int k, int i)
            {
                cmplx t0 = CC(i, 0, k), t1, t2;
                PMC(out t1, out t2, CC(i, 1, k), CC(i, 2, k));
                CH(i, k, 0).r = t0.r + t1.r;
                CH(i, k, 0).i = t0.i + t1.i;

                const int u1 = 1;
                const int u2 = 2;
                cmplx ca, cb, da, db;
                ca.r = FMA(tw1r, t1.r, t0.r);
                ca.i = FMA(tw1r, t1.i, t0.i);
                cb.i = tw1i * t2.r;
                cb.r = -(tw1i * t2.i);
                PMC(out da, out db, ca, cb);  // TODO: FMA?
                if (typeof(D) == typeof(B))
                {
                    A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                    A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                }
                else
                {
                    A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                    A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                }
            }

            internal void pass3()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                        STEP3a(k, 0);
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        STEP3a(k, 0);
                        for (int i = 1; i < ido; ++i)
                            STEP3(k, i);
                    }
                }
            }
        }

        private readonly ref struct Pass4<D>
            where D : struct, IDirection
        {
            private const int cdim = 4;

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            public Pass4(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            internal void pass4()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        cmplx t1, t2, t3, t4;
                        PMC(out t2, out t1, CC(0, 0, k), CC(0, 2, k));
                        PMC(out t3, out t4, CC(0, 1, k), CC(0, 3, k));
                        if (typeof(D) == typeof(B))
                            ROT90(ref t4);
                        else
                            ROTM90(ref t4);
                        PMC(out CH(0, k, 0), out CH(0, k, 2), t2, t3);
                        PMC(out CH(0, k, 1), out CH(0, k, 3), t1, t4);
                    }
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        {
                            cmplx t1, t2, t3, t4;
                            PMC(out t2, out t1, CC(0, 0, k), CC(0, 2, k));
                            PMC(out t3, out t4, CC(0, 1, k), CC(0, 3, k));
                            if (typeof(D) == typeof(B))
                                ROT90(ref t4);
                            else
                                ROTM90(ref t4);
                            PMC(out CH(0, k, 0), out CH(0, k, 2), t2, t3);
                            PMC(out CH(0, k, 1), out CH(0, k, 3), t1, t4);
                        }
                        for (int i = 1; i < ido; ++i)
                        {
                            cmplx c2, c3, c4, t1, t2, t3, t4;
                            cmplx cc0 = CC(i, 0, k), cc1 = CC(i, 1, k), cc2 = CC(i, 2, k), cc3 = CC(i, 3, k);
                            PMC(out t2, out t1, cc0, cc2);
                            PMC(out t3, out t4, cc1, cc3);
                            if (typeof(D) == typeof(B))
                                ROT90(ref t4);
                            else
                                ROTM90(ref t4);
                            cmplx wa0 = WA(0, i), wa1 = WA(1, i), wa2 = WA(2, i);
                            PMC(out CH(i, k, 0), out c3, t2, t3);
                            PMC(out c2, out c4, t1, t4);
                            if (typeof(D) == typeof(B))
                            {
                                A_EQ_B_MUL_C(out CH(i, k, 1), wa0, c2);
                                A_EQ_B_MUL_C(out CH(i, k, 2), wa1, c3);
                                A_EQ_B_MUL_C(out CH(i, k, 3), wa2, c4);
                            }
                            else
                            {
                                A_EQ_CB_MUL_C(out CH(i, k, 1), wa0, c2);
                                A_EQ_CB_MUL_C(out CH(i, k, 2), wa1, c3);
                                A_EQ_CB_MUL_C(out CH(i, k, 3), wa2, c4);
                            }
                        }
                    }
                }
            }
        }

        private readonly ref struct Pass5<D>
            where D : struct, IDirection
        {
            private const int cdim = 5;

            private static readonly T tw1r;
            private static readonly T tw1i;
            private static readonly T tw2r;
            private static readonly T tw2i;

            static Pass5()
            {
                (tw1i, tw1r) = T.SinCosPi(ToT(2) / ToT(5));
                (tw2i, tw2r) = T.SinCosPi(ToT(4) / ToT(5));
                if (typeof(D) == typeof(F))
                {
                    tw1i = -tw1i;
                    tw2i = -tw2i;
                }
            }

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            public Pass5(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP5a(int k, int i)
            {
                cmplx t0 = CC(i, 0, k), t1, t2, t3, t4;
                PMC(out t1, out t4, CC(i, 1, k), CC(i, 4, k));
                PMC(out t2, out t3, CC(i, 2, k), CC(i, 3, k));
                CH(i, k, 0).r = t0.r + t1.r + t2.r;
                CH(i, k, 0).i = t0.i + t1.i + t2.i;

                {
                    const int u1 = 1;
                    const int u2 = 4;
                    cmplx ca, cb;
                    ca.r = FMA(tw2r, t2.r, FMA(tw1r, t1.r, t0.r));
                    ca.i = FMA(tw2r, t2.i, FMA(tw1r, t1.i, t0.i));
                    cb.i = FMA(tw1i, t4.r, tw2i * t3.r);
                    cb.r = -FMA(tw1i, t4.i, tw2i * t3.i);
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
                {
                    const int u1 = 2;
                    const int u2 = 3;
                    cmplx ca, cb;
                    ca.r = FMA(tw1r, t2.r, FMA(tw2r, t1.r, t0.r));
                    ca.i = FMA(tw1r, t2.i, FMA(tw2r, t1.i, t0.i));
                    cb.i = FMA(tw2i, t4.r, -tw1i * t3.r);
                    cb.r = -FMA(tw2i, t4.i, -tw1i * t3.i);
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP5(int k, int i)
            {
                cmplx t0 = CC(i, 0, k), t1, t2, t3, t4;
                PMC(out t1, out t4, CC(i, 1, k), CC(i, 4, k));
                PMC(out t2, out t3, CC(i, 2, k), CC(i, 3, k));
                CH(i, k, 0).r = t0.r + t1.r + t2.r;
                CH(i, k, 0).i = t0.i + t1.i + t2.i;

                {
                    const int u1 = 1;
                    const int u2 = 4;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw2r, t2.r, FMA(tw1r, t1.r, t0.r));
                    ca.i = FMA(tw2r, t2.i, FMA(tw1r, t1.i, t0.i));
                    cb.i = FMA(tw1i, t4.r, tw2i * t3.r);
                    cb.r = -FMA(tw1i, t4.i, tw2i * t3.i);
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 2;
                    const int u2 = 3;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw1r, t2.r, FMA(tw2r, t1.r, t0.r));
                    ca.i = FMA(tw1r, t2.i, FMA(tw2r, t1.i, t0.i));
                    cb.i = FMA(tw2i, t4.r, -tw1i * t3.r);
                    cb.r = -FMA(tw2i, t4.i, -tw1i * t3.i);
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
            }

            internal void pass5()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                        STEP5a(k, 0);
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        STEP5a(k, 0);
                        for (int i = 1; i < ido; ++i)
                            STEP5(k, i);
                    }
                }
            }
        }

        private readonly ref struct Pass7<D>
            where D : struct, IDirection
        {
            private const int cdim = 7;

            private static readonly T tw1r;
            private static readonly T tw1i;
            private static readonly T tw2r;
            private static readonly T tw2i;
            private static readonly T tw3r;
            private static readonly T tw3i;

            static Pass7()
            {
                (tw1i, tw1r) = T.SinCosPi(ToT(2) / ToT(7));
                (tw2i, tw2r) = T.SinCosPi(ToT(4) / ToT(7));
                (tw3i, tw3r) = T.SinCosPi(ToT(6) / ToT(7));
                if (typeof(D) == typeof(F))
                {
                    tw1i = -tw1i;
                    tw2i = -tw2i;
                    tw3i = -tw3i;
                }
            }

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            public Pass7(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP7a(int k, int i)
            {
                cmplx t1 = CC(i, 0, k), t2, t3, t4, t5, t6, t7;
                PMC(out t2, out t7, CC(i, 1, k), CC(i, 6, k));
                PMC(out t3, out t6, CC(i, 2, k), CC(i, 5, k));
                PMC(out t4, out t5, CC(i, 3, k), CC(i, 4, k));
                CH(i, k, 0).r = t1.r + t2.r + t3.r + t4.r;
                CH(i, k, 0).i = t1.i + t2.i + t3.i + t4.i;

                {
                    const int u1 = 1;
                    const int u2 = 6;
                    cmplx ca, cb;
                    ca.r = FMA(tw3r, t4.r, FMA(tw2r, t3.r, FMA(tw1r, t2.r, t1.r)));
                    ca.i = FMA(tw3r, t4.i, FMA(tw2r, t3.i, FMA(tw1r, t2.i, t1.i)));
                    cb.i = FMA(tw1i, t7.r, FMA(tw2i, t6.r, tw3i * t5.r));
                    cb.r = -FMA(tw1i, t7.i, FMA(tw2i, t6.i, tw3i * t5.i));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);  // TODO: FMA?
                }
                {
                    const int u1 = 2;
                    const int u2 = 5;
                    cmplx ca, cb;
                    ca.r = FMA(tw1r, t4.r, FMA(tw3r, t3.r, FMA(tw2r, t2.r, t1.r)));
                    ca.i = FMA(tw1r, t4.i, FMA(tw3r, t3.i, FMA(tw2r, t2.i, t1.i)));
                    cb.i = FMA(tw2i, t7.r, FMA(-tw3i, t6.r, -tw1i * t5.r));
                    cb.r = -FMA(tw2i, t7.i, FMA(-tw3i, t6.i, -tw1i * t5.i));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);  // TODO: FMA?
                }
                {
                    const int u1 = 3;
                    const int u2 = 4;
                    cmplx ca, cb;
                    ca.r = FMA(tw2r, t4.r, FMA(tw1r, t3.r, FMA(tw3r, t2.r, t1.r)));
                    ca.i = FMA(tw2r, t4.i, FMA(tw1r, t3.i, FMA(tw3r, t2.i, t1.i)));
                    cb.i = FMA(tw3i, t7.r, FMA(-tw1i, t6.r, tw2i * t5.r));
                    cb.r = -FMA(tw3i, t7.i, FMA(-tw1i, t6.i, tw2i * t5.i));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);  // TODO: FMA?
                }
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP7(int k, int i)
            {
                cmplx t1 = CC(i, 0, k), t2, t3, t4, t5, t6, t7;
                PMC(out t2, out t7, CC(i, 1, k), CC(i, 6, k));
                PMC(out t3, out t6, CC(i, 2, k), CC(i, 5, k));
                PMC(out t4, out t5, CC(i, 3, k), CC(i, 4, k));
                CH(i, k, 0).r = t1.r + t2.r + t3.r + t4.r;
                CH(i, k, 0).i = t1.i + t2.i + t3.i + t4.i;

                {
                    const int u1 = 1;
                    const int u2 = 6;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw3r, t4.r, FMA(tw2r, t3.r, FMA(tw1r, t2.r, t1.r)));
                    ca.i = FMA(tw3r, t4.i, FMA(tw2r, t3.i, FMA(tw1r, t2.i, t1.i)));
                    cb.i = FMA(tw1i, t7.r, FMA(tw2i, t6.r, tw3i * t5.r));
                    cb.r = -FMA(tw1i, t7.i, FMA(tw2i, t6.i, tw3i * t5.i));
                    PMC(out da, out db, ca, cb);  // TODO: FMA?
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 2;
                    const int u2 = 5;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw1r, t4.r, FMA(tw3r, t3.r, FMA(tw2r, t2.r, t1.r)));
                    ca.i = FMA(tw1r, t4.i, FMA(tw3r, t3.i, FMA(tw2r, t2.i, t1.i)));
                    cb.i = FMA(tw2i, t7.r, FMA(-tw3i, t6.r, -tw1i * t5.r));
                    cb.r = -FMA(tw2i, t7.i, FMA(-tw3i, t6.i, -tw1i * t5.i));
                    PMC(out da, out db, ca, cb);  // TODO: FMA?
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 3;
                    const int u2 = 4;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw2r, t4.r, FMA(tw1r, t3.r, FMA(tw3r, t2.r, t1.r)));
                    ca.i = FMA(tw2r, t4.i, FMA(tw1r, t3.i, FMA(tw3r, t2.i, t1.i)));
                    cb.i = FMA(tw3i, t7.r, FMA(-tw1i, t6.r, tw2i * t5.r));
                    cb.r = -FMA(tw3i, t7.i, FMA(-tw1i, t6.i, tw2i * t5.i));
                    PMC(out da, out db, ca, cb);  // TODO: FMA?
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
            }

            internal void pass7()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                        STEP7a(k, 0);
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        STEP7a(k, 0);
                        for (int i = 1; i < ido; ++i)
                            STEP7(k, i);
                    }
                }
            }
        }

        private readonly ref struct Pass11<D>
            where D : struct, IDirection
        {
            private const int cdim = 11;

            private static readonly T tw1r;
            private static readonly T tw1i;
            private static readonly T tw2r;
            private static readonly T tw2i;
            private static readonly T tw3r;
            private static readonly T tw3i;
            private static readonly T tw4r;
            private static readonly T tw4i;
            private static readonly T tw5r;
            private static readonly T tw5i;

            static Pass11()
            {
                (tw1i, tw1r) = T.SinCosPi(ToT(2) / ToT(11));
                (tw2i, tw2r) = T.SinCosPi(ToT(4) / ToT(11));
                (tw3i, tw3r) = T.SinCosPi(ToT(6) / ToT(11));
                (tw4i, tw4r) = T.SinCosPi(ToT(8) / ToT(11));
                (tw5i, tw5r) = T.SinCosPi(ToT(10) / ToT(11));
                if (typeof(D) == typeof(F))
                {
                    tw1i = -tw1i;
                    tw2i = -tw2i;
                    tw3i = -tw3i;
                    tw4i = -tw4i;
                    tw5i = -tw5i;
                }
            }

            private readonly int ido;
            private readonly int l1;
            private readonly ReadOnlySpan<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;

            public Pass11(int ido, int l1, ReadOnlySpan<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa)
            {
                this.ido = ido;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref readonly cmplx WA(int x, int i) => ref wa[i - 1 + x * (ido - 1)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP11a(int k, int i)
            {
                cmplx t1 = CC(i, 0, k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                PMC(out t2, out t11, CC(i, 1, k), CC(i, 10, k));
                PMC(out t3, out t10, CC(i, 2, k), CC(i, 9, k));
                PMC(out t4, out t9, CC(i, 3, k), CC(i, 8, k));
                PMC(out t5, out t8, CC(i, 4, k), CC(i, 7, k));
                PMC(out t6, out t7, CC(i, 5, k), CC(i, 6, k));
                CH(i, k, 0).r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r;
                CH(i, k, 0).i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;

                {
                    const int u1 = 1;
                    const int u2 = 10;
                    cmplx ca, cb;
                    ca.r = FMA(tw5r, t6.r, FMA(tw4r, t5.r, FMA(tw3r, t4.r, FMA(tw2r, t3.r, FMA(tw1r, t2.r, t1.r)))));
                    ca.i = FMA(tw5r, t6.i, FMA(tw4r, t5.i, FMA(tw3r, t4.i, FMA(tw2r, t3.i, FMA(tw1r, t2.i, t1.i)))));
                    cb.i = FMA(tw1i, t11.r, FMA(tw2i, t10.r, FMA(tw3i, t9.r, FMA(tw4i, t8.r, tw5i * t7.r))));
                    cb.r = -FMA(tw1i, t11.i, FMA(tw2i, t10.i, FMA(tw3i, t9.i, FMA(tw4i, t8.i, tw5i * t7.i))));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
                {
                    const int u1 = 2;
                    const int u2 = 9;
                    cmplx ca, cb;
                    ca.r = FMA(tw1r, t6.r, FMA(tw3r, t5.r, FMA(tw5r, t4.r, FMA(tw4r, t3.r, FMA(tw2r, t2.r, t1.r)))));
                    ca.i = FMA(tw1r, t6.i, FMA(tw3r, t5.i, FMA(tw5r, t4.i, FMA(tw4r, t3.i, FMA(tw2r, t2.i, t1.i)))));
                    cb.i = FMA(tw2i, t11.r, FMA(tw4i, t10.r, FMA(-tw5i, t9.r, FMA(-tw3i, t8.r, -tw1i * t7.r))));
                    cb.r = -FMA(tw2i, t11.i, FMA(tw4i, t10.i, FMA(-tw5i, t9.i, FMA(-tw3i, t8.i, -tw1i * t7.i))));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
                {
                    const int u1 = 3;
                    const int u2 = 8;
                    cmplx ca, cb;
                    ca.r = FMA(tw4r, t6.r, FMA(tw1r, t5.r, FMA(tw2r, t4.r, FMA(tw5r, t3.r, FMA(tw3r, t2.r, t1.r)))));
                    ca.i = FMA(tw4r, t6.i, FMA(tw1r, t5.i, FMA(tw2r, t4.i, FMA(tw5r, t3.i, FMA(tw3r, t2.i, t1.i)))));
                    cb.i = FMA(tw3i, t11.r, FMA(-tw5i, t10.r, FMA(-tw2i, t9.r, FMA(tw1i, t8.r, tw4i * t7.r))));
                    cb.r = -FMA(tw3i, t11.i, FMA(-tw5i, t10.i, FMA(-tw2i, t9.i, FMA(tw1i, t8.i, tw4i * t7.i))));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
                {
                    const int u1 = 4;
                    const int u2 = 7;
                    cmplx ca, cb;
                    ca.r = FMA(tw2r, t6.r, FMA(tw5r, t5.r, FMA(tw1r, t4.r, FMA(tw3r, t3.r, FMA(tw4r, t2.r, t1.r)))));
                    ca.i = FMA(tw2r, t6.i, FMA(tw5r, t5.i, FMA(tw1r, t4.i, FMA(tw3r, t3.i, FMA(tw4r, t2.i, t1.i)))));
                    cb.i = FMA(tw4i, t11.r, FMA(-tw3i, t10.r, FMA(tw1i, t9.r, FMA(tw5i, t8.r, -tw2i * t7.r))));
                    cb.r = -FMA(tw4i, t11.i, FMA(-tw3i, t10.i, FMA(tw1i, t9.i, FMA(tw5i, t8.i, -tw2i * t7.i))));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
                {
                    const int u1 = 5;
                    const int u2 = 6;
                    cmplx ca, cb;
                    ca.r = FMA(tw3r, t6.r, FMA(tw2r, t5.r, FMA(tw4r, t4.r, FMA(tw1r, t3.r, FMA(tw5r, t2.r, t1.r)))));
                    ca.i = FMA(tw3r, t6.i, FMA(tw2r, t5.i, FMA(tw4r, t4.i, FMA(tw1r, t3.i, FMA(tw5r, t2.i, t1.i)))));
                    cb.i = FMA(tw5i, t11.r, FMA(-tw1i, t10.r, FMA(tw4i, t9.r, FMA(-tw2i, t8.r, tw3i * t7.r))));
                    cb.r = -FMA(tw5i, t11.i, FMA(-tw1i, t10.i, FMA(tw4i, t9.i, FMA(-tw2i, t8.i, tw3i * t7.i))));
                    PMC(out CH(0, k, u1), out CH(0, k, u2), ca, cb);
                }
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private void STEP11(int k, int i)
            {
                cmplx t1 = CC(i, 0, k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                PMC(out t2, out t11, CC(i, 1, k), CC(i, 10, k));
                PMC(out t3, out t10, CC(i, 2, k), CC(i, 9, k));
                PMC(out t4, out t9, CC(i, 3, k), CC(i, 8, k));
                PMC(out t5, out t8, CC(i, 4, k), CC(i, 7, k));
                PMC(out t6, out t7, CC(i, 5, k), CC(i, 6, k));
                CH(i, k, 0).r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r;
                CH(i, k, 0).i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;

                {
                    const int u1 = 1;
                    const int u2 = 10;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw1r, t2.r, FMA(tw2r, t3.r, FMA(tw3r, t4.r, FMA(tw4r, t5.r, FMA(tw5r, t6.r, t1.r)))));
                    ca.i = FMA(tw1r, t2.i, FMA(tw2r, t3.i, FMA(tw3r, t4.i, FMA(tw4r, t5.i, FMA(tw5r, t6.i, t1.i)))));
                    cb.i = FMA(tw1i, t11.r, FMA(tw2i, t10.r, FMA(tw3i, t9.r, FMA(tw4i, t8.r, tw5i * t7.r))));
                    cb.r = -FMA(tw1i, t11.i, FMA(tw2i, t10.i, FMA(tw3i, t9.i, FMA(tw4i, t8.i, tw5i * t7.i))));
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 2;
                    const int u2 = 9;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw1r, t6.r, FMA(tw3r, t5.r, FMA(tw5r, t4.r, FMA(tw4r, t3.r, FMA(tw2r, t2.r, t1.r)))));
                    ca.i = FMA(tw1r, t6.i, FMA(tw3r, t5.i, FMA(tw5r, t4.i, FMA(tw4r, t3.i, FMA(tw2r, t2.i, t1.i)))));
                    cb.i = FMA(tw2i, t11.r, FMA(tw4i, t10.r, FMA(-tw5i, t9.r, FMA(-tw3i, t8.r, -tw1i * t7.r))));
                    cb.r = -FMA(tw2i, t11.i, FMA(tw4i, t10.i, FMA(-tw5i, t9.i, FMA(-tw3i, t8.i, -tw1i * t7.i))));
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 3;
                    const int u2 = 8;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw4r, t6.r, FMA(tw1r, t5.r, FMA(tw2r, t4.r, FMA(tw5r, t3.r, FMA(tw3r, t2.r, t1.r)))));
                    ca.i = FMA(tw4r, t6.i, FMA(tw1r, t5.i, FMA(tw2r, t4.i, FMA(tw5r, t3.i, FMA(tw3r, t2.i, t1.i)))));
                    cb.i = FMA(tw3i, t11.r, FMA(-tw5i, t10.r, FMA(-tw2i, t9.r, FMA(tw1i, t8.r, tw4i * t7.r))));
                    cb.r = -FMA(tw3i, t11.i, FMA(-tw5i, t10.i, FMA(-tw2i, t9.i, FMA(tw1i, t8.i, tw4i * t7.i))));
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 4;
                    const int u2 = 7;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw2r, t6.r, FMA(tw5r, t5.r, FMA(tw1r, t4.r, FMA(tw3r, t3.r, FMA(tw4r, t2.r, t1.r)))));
                    ca.i = FMA(tw2r, t6.i, FMA(tw5r, t5.i, FMA(tw1r, t4.i, FMA(tw3r, t3.i, FMA(tw4r, t2.i, t1.i)))));
                    cb.i = FMA(tw4i, t11.r, FMA(-tw3i, t10.r, FMA(tw1i, t9.r, FMA(tw5i, t8.r, -tw2i * t7.r))));
                    cb.r = -FMA(tw4i, t11.i, FMA(-tw3i, t10.i, FMA(tw1i, t9.i, FMA(tw5i, t8.i, -tw2i * t7.i))));
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
                {
                    const int u1 = 5;
                    const int u2 = 6;
                    cmplx ca, cb, da, db;
                    ca.r = FMA(tw3r, t6.r, FMA(tw2r, t5.r, FMA(tw4r, t4.r, FMA(tw1r, t3.r, FMA(tw5r, t2.r, t1.r)))));
                    ca.i = FMA(tw3r, t6.i, FMA(tw2r, t5.i, FMA(tw4r, t4.i, FMA(tw1r, t3.i, FMA(tw5r, t2.i, t1.i)))));
                    cb.i = FMA(tw5i, t11.r, FMA(-tw1i, t10.r, FMA(tw4i, t9.r, FMA(-tw2i, t8.r, tw3i * t7.r))));
                    cb.r = -FMA(tw5i, t11.i, FMA(-tw1i, t10.i, FMA(tw4i, t9.i, FMA(-tw2i, t8.i, tw3i * t7.i))));
                    PMC(out da, out db, ca, cb);
                    if (typeof(D) == typeof(B))
                    {
                        A_EQ_B_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_B_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                    else
                    {
                        A_EQ_CB_MUL_C(out CH(i, k, u1), WA(u1 - 1, i), da);
                        A_EQ_CB_MUL_C(out CH(i, k, u2), WA(u2 - 1, i), db);
                    }
                }
            }

            internal void pass11()
            {
                if (ido == 1)
                {
                    for (int k = 0; k < l1; ++k)
                        STEP11a(k, 0);
                }
                else
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        STEP11a(k, 0);
                        for (int i = 1; i < ido; ++i)
                            STEP11(k, i);
                    }
                }
            }
        }

        private readonly ref struct PassG<D>
            where D : struct, IDirection
        {
            private int cdim => ip;

            private readonly int ido;
            private readonly int ip;
            private readonly int l1;
            private readonly Span<cmplx> cc;
            private readonly Span<cmplx> ch;
            private readonly ReadOnlySpan<cmplx> wa;
            private readonly ReadOnlySpan<cmplx> csarr;
            private readonly int idl1;

            public PassG(int ido, int ip, int l1, Span<cmplx> cc, Span<cmplx> ch, ReadOnlySpan<cmplx> wa, ReadOnlySpan<cmplx> csarr)
            {
                this.ido = ido;
                this.ip = ip;
                this.l1 = l1;
                this.cc = cc;
                this.ch = ch;
                this.wa = wa;
                this.csarr = csarr;

                idl1 = ido * l1;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH(int a, int b, int c) => ref ch[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CC(int a, int b, int c) => ref cc[a + ido * (b + cdim * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CX(int a, int b, int c) => ref cc[a + ido * (b + l1 * c)];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CX2(int a, int b) => ref cc[a + idl1 * b];

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            private ref cmplx CH2(int a, int b) => ref ch[a + idl1 * b];

            internal void passg()
            {
                int ipph = (ip + 1) / 2;

                for (int k = 0; k < l1; ++k)
                    for (int i = 0; i < ido; ++i)
                        CH(i, k, 0) = CC(i, 0, k);
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                    for (int k = 0; k < l1; ++k)
                        for (int i = 0; i < ido; ++i)
                            PMC(out CH(i, k, j), out CH(i, k, jc), CC(i, j, k), CC(i, jc, k));
                for (int k = 0; k < l1; ++k)
                {
                    for (int i = 0; i < ido; ++i)
                    {
                        cmplx tmp = CH(i, k, 0);
                        for (int j = 1; j < ipph; ++j)
                            ADDC(out tmp, tmp, CH(i, k, j));
                        CX(i, k, 0) = tmp;
                    }
                }
                for (int l = 1, lc = ip - 1; l < ipph; ++l, --lc)
                {
                    // j=0
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cmplx wal = csarr[l];
                        cmplx wa2l = csarr[2 * l];
                        Op1(out CX2(ik, l), wa2l, CH2(ik, 2), wal, CH2(ik, 1), CH2(ik, 0));
                        if (typeof(D) == typeof(B))
                            Op2B(out CX2(ik, lc), wa2l, CH2(ik, ip - 2), wal, CH2(ik, ip - 1));
                        else
                            Op2F(out CX2(ik, lc), wa2l, CH2(ik, ip - 2), wal, CH2(ik, ip - 1));

                        [MethodImpl(MethodImplOptions.AggressiveInlining)]
                        static void Op1(out cmplx d, cmplx wa2l, cmplx ch2_2, cmplx wal, cmplx ch2_1, cmplx ch2_0)
                        {
                            d.r = FMA(wa2l.r, ch2_2.r, FMA(wal.r, ch2_1.r, ch2_0.r));
                            d.i = FMA(wa2l.r, ch2_2.i, FMA(wal.r, ch2_1.i, ch2_0.i));
                        }

                        [MethodImpl(MethodImplOptions.AggressiveInlining)]
                        static void Op2B(out cmplx d, cmplx wa2l, cmplx ch2_ipm2, cmplx wal, cmplx ch2_ipm1)
                        {
                            d.r = -FMA(wa2l.i, ch2_ipm2.i, wal.i * ch2_ipm1.i);
                            d.i = FMA(wa2l.i, ch2_ipm2.r, wal.i * ch2_ipm1.r);
                        }

                        [MethodImpl(MethodImplOptions.AggressiveInlining)]
                        static void Op2F(out cmplx d, cmplx wa2l, cmplx ch2_ipm2, cmplx wal, cmplx ch2_ipm1)
                        {
                            d.r = FMA(wa2l.i, ch2_ipm2.i, wal.i * ch2_ipm1.i);
                            d.i = -FMA(wa2l.i, ch2_ipm2.r, wal.i * ch2_ipm1.r);
                        }
                    }

                    int iwal = 2 * l;
                    int j = 3, jc = ip - 3;
                    for (; j < ipph - 1; j += 2, jc -= 2)
                    {
                        iwal += l;
                        if (iwal > ip)
                            iwal -= ip;
                        cmplx xwal = csarr[iwal];
                        iwal += l;
                        if (iwal > ip)
                            iwal -= ip;
                        cmplx xwal2 = csarr[iwal];
                        for (int ik = 0; ik < idl1; ++ik)
                        {
                            Op1(ref CX2(ik, l), xwal2, CH2(ik, j + 1), xwal, CH2(ik, j));
                            if (typeof(D) == typeof(B))
                                Op2B(ref CX2(ik, lc), xwal2, CH2(ik, jc - 1), xwal, CH2(ik, jc));
                            else
                                Op2F(ref CX2(ik, lc), xwal2, CH2(ik, jc - 1), xwal, CH2(ik, jc));

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op1(ref cmplx d, cmplx xwal2, cmplx ch2_jp1, cmplx xwal, cmplx ch2_j)
                            {
                                d.r = FMA(xwal2.r, ch2_jp1.r, FMA(xwal.r, ch2_j.r, d.r));
                                d.i = FMA(xwal2.r, ch2_jp1.i, FMA(xwal.r, ch2_j.i, d.i));
                            }

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op2B(ref cmplx d, cmplx xwal2, cmplx ch2_jcm1, cmplx xwal, cmplx ch2_jc)
                            {
                                d.r = FMA(-xwal2.i, ch2_jcm1.i, FMA(-xwal.i, ch2_jc.i, d.r));
                                d.i = FMA(xwal2.i, ch2_jcm1.r, FMA(xwal.i, ch2_jc.r, d.i));
                            }

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op2F(ref cmplx d, cmplx xwal2, cmplx ch2_jcm1, cmplx xwal, cmplx ch2_jc)
                            {
                                d.r = FMA(xwal2.i, ch2_jcm1.i, FMA(xwal.i, ch2_jc.i, d.r));
                                d.i = FMA(-xwal2.i, ch2_jcm1.r, FMA(-xwal.i, ch2_jc.r, d.i));
                            }
                        }
                    }
                    for (; j < ipph; ++j, --jc)
                    {
                        iwal += l;
                        if (iwal > ip)
                            iwal -= ip;
                        cmplx xwal = csarr[iwal];
                        for (int ik = 0; ik < idl1; ++ik)
                        {
                            Op1(ref CX2(ik, l), xwal, CH2(ik, j));
                            if (typeof(D) == typeof(B))
                                Op2B(ref CX2(ik, lc), xwal, CH2(ik, jc));
                            else
                                Op2F(ref CX2(ik, lc), xwal, CH2(ik, jc));

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op1(ref cmplx d, cmplx xwal, cmplx ch2_j)
                            {
                                d.r = FMA(xwal.r, ch2_j.r, d.r);
                                d.i = FMA(xwal.r, ch2_j.i, d.i);
                            }

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op2B(ref cmplx d, cmplx xwal, cmplx ch2_jc)
                            {
                                d.r = FMA(-xwal.i, ch2_jc.i, d.r);
                                d.i = FMA(xwal.i, ch2_jc.r, d.i);
                            }

                            [MethodImpl(MethodImplOptions.AggressiveInlining)]
                            static void Op2F(ref cmplx d, cmplx xwal, cmplx ch2_jc)
                            {
                                d.r = FMA(xwal.i, ch2_jc.i, d.r);
                                d.i = FMA(-xwal.i, ch2_jc.r, d.i);
                            }
                        }
                    }
                }

                // Shuffling and twiddling
                if (ido == 1)
                {
                    for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                    {
                        for (int ik = 0; ik < idl1; ++ik)
                        {
                            cmplx t1 = CX2(ik, j), t2 = CX2(ik, jc);
                            PMC(out CX2(ik, j), out CX2(ik, jc), t1, t2);
                        }
                    }
                }
                else
                {
                    for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                    {
                        for (int k = 0; k < l1; ++k)
                        {
                            cmplx t1 = CX(0, k, j), t2 = CX(0, k, jc);
                            PMC(out CX(0, k, j), out CX(0, k, jc), t1, t2);
                            for (int i = 1; i < ido; ++i)
                            {
                                cmplx x1, x2;
                                PMC(out x1, out x2, CX(i, k, j), CX(i, k, jc));
                                int idij = (j - 1) * (ido - 1) + i - 1;  // TODO: Strength reduction?
                                if (typeof(D) == typeof(B))
                                    A_EQ_B_MUL_C(out CX(i, k, j), wa[idij], x1);
                                else
                                    A_EQ_CB_MUL_C(out CX(i, k, j), wa[idij], x1);
                                idij = (jc - 1) * (ido - 1) + i - 1;  // TODO: Strength reduction?
                                if (typeof(D) == typeof(B))
                                    A_EQ_B_MUL_C(out CX(i, k, jc), wa[idij], x2);
                                else
                                    A_EQ_CB_MUL_C(out CX(i, k, jc), wa[idij], x2);
                            }
                        }
                    }
                }
            }
        }

        internal static void pass_all<D>(ref readonly cfftp_plan plan, Span<cmplx> c, Span<cmplx> ch)
            where D : struct, IDirection
        {
            if (plan.length == 1)
                return;

            int len = plan.length;
            int l1 = 1, nf = plan.nfct;
            Span<cmplx> p1 = c, p2 = ch;

            for (int k1 = 0; k1 < nf; k1++)
            {
                int ip = plan.fct[k1].fct;
                int l2 = ip * l1;
                int ido = len / l2;
                switch (ip)
                {
                    case 4:
                        new Pass4<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass4();
                        break;
                    case 2:
                        new Pass2<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass2();
                        break;
                    case 3:
                        new Pass3<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass3();
                        break;
                    case 5:
                        new Pass5<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass5();
                        break;
                    case 7:
                        new Pass7<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass7();
                        break;
                    case 11:
                        new Pass11<D>(ido, l1, p1, p2, plan.fct[k1].tw).pass11();
                        break;
                    default:
                    {
                        new PassG<D>(ido, ip, l1, p1, p2, plan.fct[k1].tw, plan.fct[k1].tws).passg();

                        var temp2 = p1;
                        p1 = p2;
                        p2 = temp2;
                        break;
                    }
                }

                var temp = p1;
                p1 = p2;
                p2 = temp;

                l1 = l2;
            }
            if (p1 != c)
                p1.Slice(0, len).CopyTo(c);
        }

        internal static void cfftp_forward(ref readonly cfftp_plan plan, Span<cmplx> c, Span<cmplx> ch)
        {
            pass_all<F>(in plan, c, ch);
        }

        internal static void cfftp_backward(ref readonly cfftp_plan plan, Span<cmplx> c, Span<cmplx> ch)
        {
            pass_all<B>(in plan, c, ch);
        }

        private static void cfftp_factorize(ref cfftp_plan plan)
        {
            int length = plan.length;
            int nfct = 0;
            while ((length % 4) == 0)
            {
                if (nfct >= NFCT)
                    throw new UnreachableException();
                plan.fct[nfct++].fct = 4;
                length >>= 2;
            }
            if ((length % 2) == 0)
            {
                length >>= 1;
                // factor 2 should be at the front of the factor list
                if (nfct >= NFCT)
                    throw new UnreachableException();
                plan.fct[nfct++].fct = 2;
                SWAP(ref plan.fct[0].fct, ref plan.fct[nfct - 1].fct);
            }
            int maxl = (int)double.Sqrt(length) + 1;
            for (int divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
            {
                if ((length % divisor) == 0)
                {
                    while ((length % divisor) == 0)
                    {
                        if (nfct >= NFCT)
                            throw new UnreachableException();
                        plan.fct[nfct++].fct = divisor;
                        length /= divisor;
                    }
                    maxl = (int)double.Sqrt(length) + 1;
                }
            }
            if (length > 1)
                plan.fct[nfct++].fct = length;
            plan.nfct = nfct;
        }

        private static void cfftp_comp_twiddle(ref cfftp_plan plan)
        {
            int length = plan.length;
            var twid = new T[2 * length];
            sincos_2pibyn(length, twid);
            int l1 = 1;
            for (int k = 0; k < plan.nfct; ++k)
            {
                int ip = plan.fct[k].fct, ido = length / (l1 * ip);
                plan.fct[k].tw = new cmplx[(ip - 1) * (ido - 1)];
                for (int j = 1; j < ip; ++j)
                {
                    for (int i = 1; i < ido; ++i)
                    {
                        ref var tw = ref plan.fct[k].tw[(j - 1) * (ido - 1) + i - 1];
                        tw.r = twid[2 * j * l1 * i];
                        tw.i = twid[2 * j * l1 * i + 1];
                    }
                }
                if (ip > 11)
                {
                    plan.fct[k].tws = new cmplx[ip];
                    for (int j = 0; j < ip; ++j)
                    {
                        ref var tws = ref plan.fct[k].tws[j];
                        tws.r = twid[2 * j * l1 * ido];
                        tws.i = twid[2 * j * l1 * ido + 1];
                    }
                }
                l1 *= ip;
            }
        }

        internal static cfftp_plan make_cfftp_plan(int length)
        {
            var plan = new cfftp_plan();
            plan.length = length;
            if (length > 1)
            {
                cfftp_factorize(ref plan);
                cfftp_comp_twiddle(ref plan);
            }
            return plan;
        }

        [InlineArray(NFCT)]
        internal struct Fct
        {
#pragma warning disable IDE0044 // Add readonly modifier
#pragma warning disable IDE0051 // Remove unused private members
            private cfftp_fctdata _element0;
#pragma warning restore IDE0044 // Add readonly modifier
#pragma warning restore IDE0051 // Remove unused private members
        }
    }

#pragma warning restore IDE1006 // Naming Styles
}
