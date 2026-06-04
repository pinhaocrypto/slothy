/*
 * Symbolic Slothy source for NTRU+768 Phase 4 8-way parallel 32-point NTT.
 *
 * This is not a production assembly file yet.  It is the source of truth for
 * Slothy register allocation and scheduling experiments.
 *
 * Forward kernel contract:
 *
 *   natural input -> complete bit-reversed output in the row buffer
 *   input comes only from my_ntt.s Phase123 row buffers
 *   radix-2 Cooley-Tukey butterfly:
 *     t  = fqmul(high, twiddle)
 *     lo = low + t
 *     hi = low - t
 *
 * This source intentionally does not reduce the 32 loaded input vectors.
 * Phase123 feeds raw 3-point DFT outputs bounded by 3*(q-1).  The five lazy
 * CT stages stay below signed int16 range.  The final stage345 stores reduce
 * to canonical range and writes the 32 Q vectors back to row_base.  Do not use
 * this kernel as a standalone arbitrary-int16 NTT32 without restoring input
 * normalization or tightening the caller contract.
 *
 * Layout consumed by this file:
 *
 *   Q0  = [row_base + 16*0]
 *   Q1  = [row_base + 16*1]
 *   ...
 *   Q31 = [row_base + 16*31]
 *
 * Each Q vector contains 8 independent lanes.  The work index is carried by
 * the Q register number, so a CT butterfly is a vector operation between
 * work[i] and work[j].
 *
 * Region cut for practical Slothy RA:
 *
 *   1. stage12_stripe0..7 each load four work vectors, compute stage 1 and
 *      stage 2 for one stripe, then store the stage-2 values back to row_base.
 *   2. stage345_block0..3 each load one 8-vector block of stage-2 values,
 *      finish stage 3/4/5 inside that block, then reduce and store the
 *      bit-reversed output back to row_base.
 *
 * Every Slothy label boundary is also a memory boundary.  This avoids carrying
 * 32 symbolic work vectors across regions and lets each region be allocated
 * independently.
 *
 * Register contract for these symbolic regions:
 *   row_base      = x4 points at the 32-Q staged input/output block.
 *   tw_ptr        = x12 is reset before each region that loads twiddle vectors.
 *   v0.h[0]       = q = 3457 and must be reserved in Slothy config.
 *   v0.h[1]       = Barrett reduce constant used by final output reduction.
 */

row_base      .req x4
tw_ptr        .req x12

/*
 * Stage 1/2 stripes.  Stripe s touches work[s], work[s+8], work[s+16],
 * and work[s+24].  The store at the end is the handoff to stage345 blocks.
 */

    .global ntt32_8way
    .global _ntt32_8way
ntt32_8way:
_ntt32_8way:
    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe0_slothy_start:
    // Stage 1/2 stripe 0: work[0], work[8], work[16], work[24].
    ldr Q<w00_s0_r0>, [row_base, #16*0]
    ldr Q<w08_s0_r0>, [row_base, #16*8]
    ldr Q<w16_s0_r0>, [row_base, #16*16]
    ldr Q<w24_s0_r0>, [row_base, #16*24]

    ldr Q<tw_s12_r0_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r0_p00_07>, [tw_ptr], #16

    // work[0] with work[16], CT twiddle power 0.
    sqrdmulh V<q16_s1_r0>.8h, V<w16_s0_r0>.8h, V<pre_s12_r0_p00_07>.h[0]
    mul V<t16_s1_r0>.8h, V<w16_s0_r0>.8h, V<tw_s12_r0_p00_07>.h[0]
    mls V<t16_s1_r0>.8h, V<q16_s1_r0>.8h, v0.h[0]
    add V<w00_s1_r0>.8h, V<w00_s0_r0>.8h, V<t16_s1_r0>.8h
    sub V<w16_s1_r0>.8h, V<w00_s0_r0>.8h, V<t16_s1_r0>.8h
    // work[8] with work[24], CT twiddle power 0.
    sqrdmulh V<q24_s1_r0>.8h, V<w24_s0_r0>.8h, V<pre_s12_r0_p00_07>.h[0]
    mul V<t24_s1_r0>.8h, V<w24_s0_r0>.8h, V<tw_s12_r0_p00_07>.h[0]
    mls V<t24_s1_r0>.8h, V<q24_s1_r0>.8h, v0.h[0]
    add V<w08_s1_r0>.8h, V<w08_s0_r0>.8h, V<t24_s1_r0>.8h
    sub V<w24_s1_r0>.8h, V<w08_s0_r0>.8h, V<t24_s1_r0>.8h
    // work[0] with work[8], CT twiddle power 0.
    sqrdmulh V<q08_s2_r0>.8h, V<w08_s1_r0>.8h, V<pre_s12_r0_p00_07>.h[0]
    mul V<t08_s2_r0>.8h, V<w08_s1_r0>.8h, V<tw_s12_r0_p00_07>.h[0]
    mls V<t08_s2_r0>.8h, V<q08_s2_r0>.8h, v0.h[0]
    add V<w00_s2_r0>.8h, V<w00_s1_r0>.8h, V<t08_s2_r0>.8h
    sub V<w08_s2_r0>.8h, V<w00_s1_r0>.8h, V<t08_s2_r0>.8h
    ldr Q<tw_s12_r0_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r0_p08_15>, [tw_ptr], #16

    // work[16] with work[24], CT twiddle power 8.
    sqrdmulh V<q24_s2_r0>.8h, V<w24_s1_r0>.8h, V<pre_s12_r0_p08_15>.h[0]
    mul V<t24_s2_r0>.8h, V<w24_s1_r0>.8h, V<tw_s12_r0_p08_15>.h[0]
    mls V<t24_s2_r0>.8h, V<q24_s2_r0>.8h, v0.h[0]
    add V<w16_s2_r0>.8h, V<w16_s1_r0>.8h, V<t24_s2_r0>.8h
    sub V<w24_s2_r0>.8h, V<w16_s1_r0>.8h, V<t24_s2_r0>.8h
    str Q<w00_s2_r0>, [row_base, #16*0]
    str Q<w08_s2_r0>, [row_base, #16*8]
    str Q<w16_s2_r0>, [row_base, #16*16]
    str Q<w24_s2_r0>, [row_base, #16*24]
_ntt32_stage12_stripe0_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe1_slothy_start:
    // Stage 1/2 stripe 1: work[1], work[9], work[17], work[25].
    ldr Q<w01_s0_r1>, [row_base, #16*1]
    ldr Q<w09_s0_r1>, [row_base, #16*9]
    ldr Q<w17_s0_r1>, [row_base, #16*17]
    ldr Q<w25_s0_r1>, [row_base, #16*25]

    ldr Q<tw_s12_r1_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r1_p00_07>, [tw_ptr], #16

    // work[1] with work[17], CT twiddle power 0.
    sqrdmulh V<q17_s1_r1>.8h, V<w17_s0_r1>.8h, V<pre_s12_r1_p00_07>.h[0]
    mul V<t17_s1_r1>.8h, V<w17_s0_r1>.8h, V<tw_s12_r1_p00_07>.h[0]
    mls V<t17_s1_r1>.8h, V<q17_s1_r1>.8h, v0.h[0]
    add V<w01_s1_r1>.8h, V<w01_s0_r1>.8h, V<t17_s1_r1>.8h
    sub V<w17_s1_r1>.8h, V<w01_s0_r1>.8h, V<t17_s1_r1>.8h
    // work[9] with work[25], CT twiddle power 0.
    sqrdmulh V<q25_s1_r1>.8h, V<w25_s0_r1>.8h, V<pre_s12_r1_p00_07>.h[0]
    mul V<t25_s1_r1>.8h, V<w25_s0_r1>.8h, V<tw_s12_r1_p00_07>.h[0]
    mls V<t25_s1_r1>.8h, V<q25_s1_r1>.8h, v0.h[0]
    add V<w09_s1_r1>.8h, V<w09_s0_r1>.8h, V<t25_s1_r1>.8h
    sub V<w25_s1_r1>.8h, V<w09_s0_r1>.8h, V<t25_s1_r1>.8h
    // work[1] with work[9], CT twiddle power 0.
    sqrdmulh V<q09_s2_r1>.8h, V<w09_s1_r1>.8h, V<pre_s12_r1_p00_07>.h[0]
    mul V<t09_s2_r1>.8h, V<w09_s1_r1>.8h, V<tw_s12_r1_p00_07>.h[0]
    mls V<t09_s2_r1>.8h, V<q09_s2_r1>.8h, v0.h[0]
    add V<w01_s2_r1>.8h, V<w01_s1_r1>.8h, V<t09_s2_r1>.8h
    sub V<w09_s2_r1>.8h, V<w01_s1_r1>.8h, V<t09_s2_r1>.8h
    ldr Q<tw_s12_r1_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r1_p08_15>, [tw_ptr], #16

    // work[17] with work[25], CT twiddle power 8.
    sqrdmulh V<q25_s2_r1>.8h, V<w25_s1_r1>.8h, V<pre_s12_r1_p08_15>.h[0]
    mul V<t25_s2_r1>.8h, V<w25_s1_r1>.8h, V<tw_s12_r1_p08_15>.h[0]
    mls V<t25_s2_r1>.8h, V<q25_s2_r1>.8h, v0.h[0]
    add V<w17_s2_r1>.8h, V<w17_s1_r1>.8h, V<t25_s2_r1>.8h
    sub V<w25_s2_r1>.8h, V<w17_s1_r1>.8h, V<t25_s2_r1>.8h
    str Q<w01_s2_r1>, [row_base, #16*1]
    str Q<w09_s2_r1>, [row_base, #16*9]
    str Q<w17_s2_r1>, [row_base, #16*17]
    str Q<w25_s2_r1>, [row_base, #16*25]
_ntt32_stage12_stripe1_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe2_slothy_start:
    // Stage 1/2 stripe 2: work[2], work[10], work[18], work[26].
    ldr Q<w02_s0_r2>, [row_base, #16*2]
    ldr Q<w10_s0_r2>, [row_base, #16*10]
    ldr Q<w18_s0_r2>, [row_base, #16*18]
    ldr Q<w26_s0_r2>, [row_base, #16*26]

    ldr Q<tw_s12_r2_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r2_p00_07>, [tw_ptr], #16

    // work[2] with work[18], CT twiddle power 0.
    sqrdmulh V<q18_s1_r2>.8h, V<w18_s0_r2>.8h, V<pre_s12_r2_p00_07>.h[0]
    mul V<t18_s1_r2>.8h, V<w18_s0_r2>.8h, V<tw_s12_r2_p00_07>.h[0]
    mls V<t18_s1_r2>.8h, V<q18_s1_r2>.8h, v0.h[0]
    add V<w02_s1_r2>.8h, V<w02_s0_r2>.8h, V<t18_s1_r2>.8h
    sub V<w18_s1_r2>.8h, V<w02_s0_r2>.8h, V<t18_s1_r2>.8h
    // work[10] with work[26], CT twiddle power 0.
    sqrdmulh V<q26_s1_r2>.8h, V<w26_s0_r2>.8h, V<pre_s12_r2_p00_07>.h[0]
    mul V<t26_s1_r2>.8h, V<w26_s0_r2>.8h, V<tw_s12_r2_p00_07>.h[0]
    mls V<t26_s1_r2>.8h, V<q26_s1_r2>.8h, v0.h[0]
    add V<w10_s1_r2>.8h, V<w10_s0_r2>.8h, V<t26_s1_r2>.8h
    sub V<w26_s1_r2>.8h, V<w10_s0_r2>.8h, V<t26_s1_r2>.8h
    // work[2] with work[10], CT twiddle power 0.
    sqrdmulh V<q10_s2_r2>.8h, V<w10_s1_r2>.8h, V<pre_s12_r2_p00_07>.h[0]
    mul V<t10_s2_r2>.8h, V<w10_s1_r2>.8h, V<tw_s12_r2_p00_07>.h[0]
    mls V<t10_s2_r2>.8h, V<q10_s2_r2>.8h, v0.h[0]
    add V<w02_s2_r2>.8h, V<w02_s1_r2>.8h, V<t10_s2_r2>.8h
    sub V<w10_s2_r2>.8h, V<w02_s1_r2>.8h, V<t10_s2_r2>.8h
    ldr Q<tw_s12_r2_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r2_p08_15>, [tw_ptr], #16

    // work[18] with work[26], CT twiddle power 8.
    sqrdmulh V<q26_s2_r2>.8h, V<w26_s1_r2>.8h, V<pre_s12_r2_p08_15>.h[0]
    mul V<t26_s2_r2>.8h, V<w26_s1_r2>.8h, V<tw_s12_r2_p08_15>.h[0]
    mls V<t26_s2_r2>.8h, V<q26_s2_r2>.8h, v0.h[0]
    add V<w18_s2_r2>.8h, V<w18_s1_r2>.8h, V<t26_s2_r2>.8h
    sub V<w26_s2_r2>.8h, V<w18_s1_r2>.8h, V<t26_s2_r2>.8h
    str Q<w02_s2_r2>, [row_base, #16*2]
    str Q<w10_s2_r2>, [row_base, #16*10]
    str Q<w18_s2_r2>, [row_base, #16*18]
    str Q<w26_s2_r2>, [row_base, #16*26]
_ntt32_stage12_stripe2_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe3_slothy_start:
    // Stage 1/2 stripe 3: work[3], work[11], work[19], work[27].
    ldr Q<w03_s0_r3>, [row_base, #16*3]
    ldr Q<w11_s0_r3>, [row_base, #16*11]
    ldr Q<w19_s0_r3>, [row_base, #16*19]
    ldr Q<w27_s0_r3>, [row_base, #16*27]

    ldr Q<tw_s12_r3_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r3_p00_07>, [tw_ptr], #16

    // work[3] with work[19], CT twiddle power 0.
    sqrdmulh V<q19_s1_r3>.8h, V<w19_s0_r3>.8h, V<pre_s12_r3_p00_07>.h[0]
    mul V<t19_s1_r3>.8h, V<w19_s0_r3>.8h, V<tw_s12_r3_p00_07>.h[0]
    mls V<t19_s1_r3>.8h, V<q19_s1_r3>.8h, v0.h[0]
    add V<w03_s1_r3>.8h, V<w03_s0_r3>.8h, V<t19_s1_r3>.8h
    sub V<w19_s1_r3>.8h, V<w03_s0_r3>.8h, V<t19_s1_r3>.8h
    // work[11] with work[27], CT twiddle power 0.
    sqrdmulh V<q27_s1_r3>.8h, V<w27_s0_r3>.8h, V<pre_s12_r3_p00_07>.h[0]
    mul V<t27_s1_r3>.8h, V<w27_s0_r3>.8h, V<tw_s12_r3_p00_07>.h[0]
    mls V<t27_s1_r3>.8h, V<q27_s1_r3>.8h, v0.h[0]
    add V<w11_s1_r3>.8h, V<w11_s0_r3>.8h, V<t27_s1_r3>.8h
    sub V<w27_s1_r3>.8h, V<w11_s0_r3>.8h, V<t27_s1_r3>.8h
    // work[3] with work[11], CT twiddle power 0.
    sqrdmulh V<q11_s2_r3>.8h, V<w11_s1_r3>.8h, V<pre_s12_r3_p00_07>.h[0]
    mul V<t11_s2_r3>.8h, V<w11_s1_r3>.8h, V<tw_s12_r3_p00_07>.h[0]
    mls V<t11_s2_r3>.8h, V<q11_s2_r3>.8h, v0.h[0]
    add V<w03_s2_r3>.8h, V<w03_s1_r3>.8h, V<t11_s2_r3>.8h
    sub V<w11_s2_r3>.8h, V<w03_s1_r3>.8h, V<t11_s2_r3>.8h
    ldr Q<tw_s12_r3_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r3_p08_15>, [tw_ptr], #16

    // work[19] with work[27], CT twiddle power 8.
    sqrdmulh V<q27_s2_r3>.8h, V<w27_s1_r3>.8h, V<pre_s12_r3_p08_15>.h[0]
    mul V<t27_s2_r3>.8h, V<w27_s1_r3>.8h, V<tw_s12_r3_p08_15>.h[0]
    mls V<t27_s2_r3>.8h, V<q27_s2_r3>.8h, v0.h[0]
    add V<w19_s2_r3>.8h, V<w19_s1_r3>.8h, V<t27_s2_r3>.8h
    sub V<w27_s2_r3>.8h, V<w19_s1_r3>.8h, V<t27_s2_r3>.8h
    str Q<w03_s2_r3>, [row_base, #16*3]
    str Q<w11_s2_r3>, [row_base, #16*11]
    str Q<w19_s2_r3>, [row_base, #16*19]
    str Q<w27_s2_r3>, [row_base, #16*27]
_ntt32_stage12_stripe3_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe4_slothy_start:
    // Stage 1/2 stripe 4: work[4], work[12], work[20], work[28].
    ldr Q<w04_s0_r4>, [row_base, #16*4]
    ldr Q<w12_s0_r4>, [row_base, #16*12]
    ldr Q<w20_s0_r4>, [row_base, #16*20]
    ldr Q<w28_s0_r4>, [row_base, #16*28]

    ldr Q<tw_s12_r4_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r4_p00_07>, [tw_ptr], #16

    // work[4] with work[20], CT twiddle power 0.
    sqrdmulh V<q20_s1_r4>.8h, V<w20_s0_r4>.8h, V<pre_s12_r4_p00_07>.h[0]
    mul V<t20_s1_r4>.8h, V<w20_s0_r4>.8h, V<tw_s12_r4_p00_07>.h[0]
    mls V<t20_s1_r4>.8h, V<q20_s1_r4>.8h, v0.h[0]
    add V<w04_s1_r4>.8h, V<w04_s0_r4>.8h, V<t20_s1_r4>.8h
    sub V<w20_s1_r4>.8h, V<w04_s0_r4>.8h, V<t20_s1_r4>.8h
    // work[12] with work[28], CT twiddle power 0.
    sqrdmulh V<q28_s1_r4>.8h, V<w28_s0_r4>.8h, V<pre_s12_r4_p00_07>.h[0]
    mul V<t28_s1_r4>.8h, V<w28_s0_r4>.8h, V<tw_s12_r4_p00_07>.h[0]
    mls V<t28_s1_r4>.8h, V<q28_s1_r4>.8h, v0.h[0]
    add V<w12_s1_r4>.8h, V<w12_s0_r4>.8h, V<t28_s1_r4>.8h
    sub V<w28_s1_r4>.8h, V<w12_s0_r4>.8h, V<t28_s1_r4>.8h
    // work[4] with work[12], CT twiddle power 0.
    sqrdmulh V<q12_s2_r4>.8h, V<w12_s1_r4>.8h, V<pre_s12_r4_p00_07>.h[0]
    mul V<t12_s2_r4>.8h, V<w12_s1_r4>.8h, V<tw_s12_r4_p00_07>.h[0]
    mls V<t12_s2_r4>.8h, V<q12_s2_r4>.8h, v0.h[0]
    add V<w04_s2_r4>.8h, V<w04_s1_r4>.8h, V<t12_s2_r4>.8h
    sub V<w12_s2_r4>.8h, V<w04_s1_r4>.8h, V<t12_s2_r4>.8h
    ldr Q<tw_s12_r4_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r4_p08_15>, [tw_ptr], #16

    // work[20] with work[28], CT twiddle power 8.
    sqrdmulh V<q28_s2_r4>.8h, V<w28_s1_r4>.8h, V<pre_s12_r4_p08_15>.h[0]
    mul V<t28_s2_r4>.8h, V<w28_s1_r4>.8h, V<tw_s12_r4_p08_15>.h[0]
    mls V<t28_s2_r4>.8h, V<q28_s2_r4>.8h, v0.h[0]
    add V<w20_s2_r4>.8h, V<w20_s1_r4>.8h, V<t28_s2_r4>.8h
    sub V<w28_s2_r4>.8h, V<w20_s1_r4>.8h, V<t28_s2_r4>.8h
    str Q<w04_s2_r4>, [row_base, #16*4]
    str Q<w12_s2_r4>, [row_base, #16*12]
    str Q<w20_s2_r4>, [row_base, #16*20]
    str Q<w28_s2_r4>, [row_base, #16*28]
_ntt32_stage12_stripe4_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe5_slothy_start:
    // Stage 1/2 stripe 5: work[5], work[13], work[21], work[29].
    ldr Q<w05_s0_r5>, [row_base, #16*5]
    ldr Q<w13_s0_r5>, [row_base, #16*13]
    ldr Q<w21_s0_r5>, [row_base, #16*21]
    ldr Q<w29_s0_r5>, [row_base, #16*29]

    ldr Q<tw_s12_r5_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r5_p00_07>, [tw_ptr], #16

    // work[5] with work[21], CT twiddle power 0.
    sqrdmulh V<q21_s1_r5>.8h, V<w21_s0_r5>.8h, V<pre_s12_r5_p00_07>.h[0]
    mul V<t21_s1_r5>.8h, V<w21_s0_r5>.8h, V<tw_s12_r5_p00_07>.h[0]
    mls V<t21_s1_r5>.8h, V<q21_s1_r5>.8h, v0.h[0]
    add V<w05_s1_r5>.8h, V<w05_s0_r5>.8h, V<t21_s1_r5>.8h
    sub V<w21_s1_r5>.8h, V<w05_s0_r5>.8h, V<t21_s1_r5>.8h
    // work[13] with work[29], CT twiddle power 0.
    sqrdmulh V<q29_s1_r5>.8h, V<w29_s0_r5>.8h, V<pre_s12_r5_p00_07>.h[0]
    mul V<t29_s1_r5>.8h, V<w29_s0_r5>.8h, V<tw_s12_r5_p00_07>.h[0]
    mls V<t29_s1_r5>.8h, V<q29_s1_r5>.8h, v0.h[0]
    add V<w13_s1_r5>.8h, V<w13_s0_r5>.8h, V<t29_s1_r5>.8h
    sub V<w29_s1_r5>.8h, V<w13_s0_r5>.8h, V<t29_s1_r5>.8h
    // work[5] with work[13], CT twiddle power 0.
    sqrdmulh V<q13_s2_r5>.8h, V<w13_s1_r5>.8h, V<pre_s12_r5_p00_07>.h[0]
    mul V<t13_s2_r5>.8h, V<w13_s1_r5>.8h, V<tw_s12_r5_p00_07>.h[0]
    mls V<t13_s2_r5>.8h, V<q13_s2_r5>.8h, v0.h[0]
    add V<w05_s2_r5>.8h, V<w05_s1_r5>.8h, V<t13_s2_r5>.8h
    sub V<w13_s2_r5>.8h, V<w05_s1_r5>.8h, V<t13_s2_r5>.8h
    ldr Q<tw_s12_r5_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r5_p08_15>, [tw_ptr], #16

    // work[21] with work[29], CT twiddle power 8.
    sqrdmulh V<q29_s2_r5>.8h, V<w29_s1_r5>.8h, V<pre_s12_r5_p08_15>.h[0]
    mul V<t29_s2_r5>.8h, V<w29_s1_r5>.8h, V<tw_s12_r5_p08_15>.h[0]
    mls V<t29_s2_r5>.8h, V<q29_s2_r5>.8h, v0.h[0]
    add V<w21_s2_r5>.8h, V<w21_s1_r5>.8h, V<t29_s2_r5>.8h
    sub V<w29_s2_r5>.8h, V<w21_s1_r5>.8h, V<t29_s2_r5>.8h
    str Q<w05_s2_r5>, [row_base, #16*5]
    str Q<w13_s2_r5>, [row_base, #16*13]
    str Q<w21_s2_r5>, [row_base, #16*21]
    str Q<w29_s2_r5>, [row_base, #16*29]
_ntt32_stage12_stripe5_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe6_slothy_start:
    // Stage 1/2 stripe 6: work[6], work[14], work[22], work[30].
    ldr Q<w06_s0_r6>, [row_base, #16*6]
    ldr Q<w14_s0_r6>, [row_base, #16*14]
    ldr Q<w22_s0_r6>, [row_base, #16*22]
    ldr Q<w30_s0_r6>, [row_base, #16*30]

    ldr Q<tw_s12_r6_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r6_p00_07>, [tw_ptr], #16

    // work[6] with work[22], CT twiddle power 0.
    sqrdmulh V<q22_s1_r6>.8h, V<w22_s0_r6>.8h, V<pre_s12_r6_p00_07>.h[0]
    mul V<t22_s1_r6>.8h, V<w22_s0_r6>.8h, V<tw_s12_r6_p00_07>.h[0]
    mls V<t22_s1_r6>.8h, V<q22_s1_r6>.8h, v0.h[0]
    add V<w06_s1_r6>.8h, V<w06_s0_r6>.8h, V<t22_s1_r6>.8h
    sub V<w22_s1_r6>.8h, V<w06_s0_r6>.8h, V<t22_s1_r6>.8h
    // work[14] with work[30], CT twiddle power 0.
    sqrdmulh V<q30_s1_r6>.8h, V<w30_s0_r6>.8h, V<pre_s12_r6_p00_07>.h[0]
    mul V<t30_s1_r6>.8h, V<w30_s0_r6>.8h, V<tw_s12_r6_p00_07>.h[0]
    mls V<t30_s1_r6>.8h, V<q30_s1_r6>.8h, v0.h[0]
    add V<w14_s1_r6>.8h, V<w14_s0_r6>.8h, V<t30_s1_r6>.8h
    sub V<w30_s1_r6>.8h, V<w14_s0_r6>.8h, V<t30_s1_r6>.8h
    // work[6] with work[14], CT twiddle power 0.
    sqrdmulh V<q14_s2_r6>.8h, V<w14_s1_r6>.8h, V<pre_s12_r6_p00_07>.h[0]
    mul V<t14_s2_r6>.8h, V<w14_s1_r6>.8h, V<tw_s12_r6_p00_07>.h[0]
    mls V<t14_s2_r6>.8h, V<q14_s2_r6>.8h, v0.h[0]
    add V<w06_s2_r6>.8h, V<w06_s1_r6>.8h, V<t14_s2_r6>.8h
    sub V<w14_s2_r6>.8h, V<w06_s1_r6>.8h, V<t14_s2_r6>.8h
    ldr Q<tw_s12_r6_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r6_p08_15>, [tw_ptr], #16

    // work[22] with work[30], CT twiddle power 8.
    sqrdmulh V<q30_s2_r6>.8h, V<w30_s1_r6>.8h, V<pre_s12_r6_p08_15>.h[0]
    mul V<t30_s2_r6>.8h, V<w30_s1_r6>.8h, V<tw_s12_r6_p08_15>.h[0]
    mls V<t30_s2_r6>.8h, V<q30_s2_r6>.8h, v0.h[0]
    add V<w22_s2_r6>.8h, V<w22_s1_r6>.8h, V<t30_s2_r6>.8h
    sub V<w30_s2_r6>.8h, V<w22_s1_r6>.8h, V<t30_s2_r6>.8h
    str Q<w06_s2_r6>, [row_base, #16*6]
    str Q<w14_s2_r6>, [row_base, #16*14]
    str Q<w22_s2_r6>, [row_base, #16*22]
    str Q<w30_s2_r6>, [row_base, #16*30]
_ntt32_stage12_stripe6_slothy_end:

    adr tw_ptr, ntt32_twiddle_vecs
_ntt32_stage12_stripe7_slothy_start:
    // Stage 1/2 stripe 7: work[7], work[15], work[23], work[31].
    ldr Q<w07_s0_r7>, [row_base, #16*7]
    ldr Q<w15_s0_r7>, [row_base, #16*15]
    ldr Q<w23_s0_r7>, [row_base, #16*23]
    ldr Q<w31_s0_r7>, [row_base, #16*31]

    ldr Q<tw_s12_r7_p00_07>,  [tw_ptr], #16
    ldr Q<pre_s12_r7_p00_07>, [tw_ptr], #16

    // work[7] with work[23], CT twiddle power 0.
    sqrdmulh V<q23_s1_r7>.8h, V<w23_s0_r7>.8h, V<pre_s12_r7_p00_07>.h[0]
    mul V<t23_s1_r7>.8h, V<w23_s0_r7>.8h, V<tw_s12_r7_p00_07>.h[0]
    mls V<t23_s1_r7>.8h, V<q23_s1_r7>.8h, v0.h[0]
    add V<w07_s1_r7>.8h, V<w07_s0_r7>.8h, V<t23_s1_r7>.8h
    sub V<w23_s1_r7>.8h, V<w07_s0_r7>.8h, V<t23_s1_r7>.8h
    // work[15] with work[31], CT twiddle power 0.
    sqrdmulh V<q31_s1_r7>.8h, V<w31_s0_r7>.8h, V<pre_s12_r7_p00_07>.h[0]
    mul V<t31_s1_r7>.8h, V<w31_s0_r7>.8h, V<tw_s12_r7_p00_07>.h[0]
    mls V<t31_s1_r7>.8h, V<q31_s1_r7>.8h, v0.h[0]
    add V<w15_s1_r7>.8h, V<w15_s0_r7>.8h, V<t31_s1_r7>.8h
    sub V<w31_s1_r7>.8h, V<w15_s0_r7>.8h, V<t31_s1_r7>.8h
    // work[7] with work[15], CT twiddle power 0.
    sqrdmulh V<q15_s2_r7>.8h, V<w15_s1_r7>.8h, V<pre_s12_r7_p00_07>.h[0]
    mul V<t15_s2_r7>.8h, V<w15_s1_r7>.8h, V<tw_s12_r7_p00_07>.h[0]
    mls V<t15_s2_r7>.8h, V<q15_s2_r7>.8h, v0.h[0]
    add V<w07_s2_r7>.8h, V<w07_s1_r7>.8h, V<t15_s2_r7>.8h
    sub V<w15_s2_r7>.8h, V<w07_s1_r7>.8h, V<t15_s2_r7>.8h
    ldr Q<tw_s12_r7_p08_15>,  [tw_ptr], #16
    ldr Q<pre_s12_r7_p08_15>, [tw_ptr], #16

    // work[23] with work[31], CT twiddle power 8.
    sqrdmulh V<q31_s2_r7>.8h, V<w31_s1_r7>.8h, V<pre_s12_r7_p08_15>.h[0]
    mul V<t31_s2_r7>.8h, V<w31_s1_r7>.8h, V<tw_s12_r7_p08_15>.h[0]
    mls V<t31_s2_r7>.8h, V<q31_s2_r7>.8h, v0.h[0]
    add V<w23_s2_r7>.8h, V<w23_s1_r7>.8h, V<t31_s2_r7>.8h
    sub V<w31_s2_r7>.8h, V<w23_s1_r7>.8h, V<t31_s2_r7>.8h
    str Q<w07_s2_r7>, [row_base, #16*7]
    str Q<w15_s2_r7>, [row_base, #16*15]
    str Q<w23_s2_r7>, [row_base, #16*23]
    str Q<w31_s2_r7>, [row_base, #16*31]
_ntt32_stage12_stripe7_slothy_end:

/*
 * Stage 3/4/5 blocks.  After stage 2, each 8-work block is independent:
 * block0 = work[0..7], block1 = work[8..15], block2 = work[16..23],
 * block3 = work[24..31].
 */

    adr tw_ptr, ntt32_twiddle_stage3
_ntt32_stage345_block0_slothy_start:
    // Stage 3/4/5 block 0: work[0..7].
    ldr Q<w00_s2_b0>, [row_base, #16*0]
    ldr Q<w01_s2_b0>, [row_base, #16*1]
    ldr Q<w02_s2_b0>, [row_base, #16*2]
    ldr Q<w03_s2_b0>, [row_base, #16*3]
    ldr Q<w04_s2_b0>, [row_base, #16*4]
    ldr Q<w05_s2_b0>, [row_base, #16*5]
    ldr Q<w06_s2_b0>, [row_base, #16*6]
    ldr Q<w07_s2_b0>, [row_base, #16*7]

    ldr Q<tw_s345_b0_stage3>,  [tw_ptr], #16
    ldr Q<pre_s345_b0_stage3>, [tw_ptr], #16

    // work[0] with work[4], CT twiddle power 0.
    sqrdmulh V<q04_s3_b0>.8h, V<w04_s2_b0>.8h, V<pre_s345_b0_stage3>.h[0]
    mul V<t04_s3_b0>.8h, V<w04_s2_b0>.8h, V<tw_s345_b0_stage3>.h[0]
    mls V<t04_s3_b0>.8h, V<q04_s3_b0>.8h, v0.h[0]
    add V<w00_s3_b0>.8h, V<w00_s2_b0>.8h, V<t04_s3_b0>.8h
    sub V<w04_s3_b0>.8h, V<w00_s2_b0>.8h, V<t04_s3_b0>.8h
    // work[1] with work[5], CT twiddle power 0.
    sqrdmulh V<q05_s3_b0>.8h, V<w05_s2_b0>.8h, V<pre_s345_b0_stage3>.h[0]
    mul V<t05_s3_b0>.8h, V<w05_s2_b0>.8h, V<tw_s345_b0_stage3>.h[0]
    mls V<t05_s3_b0>.8h, V<q05_s3_b0>.8h, v0.h[0]
    add V<w01_s3_b0>.8h, V<w01_s2_b0>.8h, V<t05_s3_b0>.8h
    sub V<w05_s3_b0>.8h, V<w01_s2_b0>.8h, V<t05_s3_b0>.8h
    // work[2] with work[6], CT twiddle power 0.
    sqrdmulh V<q06_s3_b0>.8h, V<w06_s2_b0>.8h, V<pre_s345_b0_stage3>.h[0]
    mul V<t06_s3_b0>.8h, V<w06_s2_b0>.8h, V<tw_s345_b0_stage3>.h[0]
    mls V<t06_s3_b0>.8h, V<q06_s3_b0>.8h, v0.h[0]
    add V<w02_s3_b0>.8h, V<w02_s2_b0>.8h, V<t06_s3_b0>.8h
    sub V<w06_s3_b0>.8h, V<w02_s2_b0>.8h, V<t06_s3_b0>.8h
    // work[3] with work[7], CT twiddle power 0.
    sqrdmulh V<q07_s3_b0>.8h, V<w07_s2_b0>.8h, V<pre_s345_b0_stage3>.h[0]
    mul V<t07_s3_b0>.8h, V<w07_s2_b0>.8h, V<tw_s345_b0_stage3>.h[0]
    mls V<t07_s3_b0>.8h, V<q07_s3_b0>.8h, v0.h[0]
    add V<w03_s3_b0>.8h, V<w03_s2_b0>.8h, V<t07_s3_b0>.8h
    sub V<w07_s3_b0>.8h, V<w03_s2_b0>.8h, V<t07_s3_b0>.8h
    ldr Q<tw_s345_b0_stage45_low>,  [tw_ptr], #16
    ldr Q<pre_s345_b0_stage45_low>, [tw_ptr], #16

    // work[0] with work[2], CT twiddle power 0.
    sqrdmulh V<q02_s4_b0>.8h, V<w02_s3_b0>.8h, V<pre_s345_b0_stage45_low>.h[0]
    mul V<t02_s4_b0>.8h, V<w02_s3_b0>.8h, V<tw_s345_b0_stage45_low>.h[0]
    mls V<t02_s4_b0>.8h, V<q02_s4_b0>.8h, v0.h[0]
    add V<w00_s4_b0>.8h, V<w00_s3_b0>.8h, V<t02_s4_b0>.8h
    sub V<w02_s4_b0>.8h, V<w00_s3_b0>.8h, V<t02_s4_b0>.8h
    // work[1] with work[3], CT twiddle power 0.
    sqrdmulh V<q03_s4_b0>.8h, V<w03_s3_b0>.8h, V<pre_s345_b0_stage45_low>.h[0]
    mul V<t03_s4_b0>.8h, V<w03_s3_b0>.8h, V<tw_s345_b0_stage45_low>.h[0]
    mls V<t03_s4_b0>.8h, V<q03_s4_b0>.8h, v0.h[0]
    add V<w01_s4_b0>.8h, V<w01_s3_b0>.8h, V<t03_s4_b0>.8h
    sub V<w03_s4_b0>.8h, V<w01_s3_b0>.8h, V<t03_s4_b0>.8h
    // work[4] with work[6], CT twiddle power 8.
    sqrdmulh V<q06_s4_b0>.8h, V<w06_s3_b0>.8h, V<pre_s345_b0_stage45_low>.h[1]
    mul V<t06_s4_b0>.8h, V<w06_s3_b0>.8h, V<tw_s345_b0_stage45_low>.h[1]
    mls V<t06_s4_b0>.8h, V<q06_s4_b0>.8h, v0.h[0]
    add V<w04_s4_b0>.8h, V<w04_s3_b0>.8h, V<t06_s4_b0>.8h
    sub V<w06_s4_b0>.8h, V<w04_s3_b0>.8h, V<t06_s4_b0>.8h
    // work[5] with work[7], CT twiddle power 8.
    sqrdmulh V<q07_s4_b0>.8h, V<w07_s3_b0>.8h, V<pre_s345_b0_stage45_low>.h[1]
    mul V<t07_s4_b0>.8h, V<w07_s3_b0>.8h, V<tw_s345_b0_stage45_low>.h[1]
    mls V<t07_s4_b0>.8h, V<q07_s4_b0>.8h, v0.h[0]
    add V<w05_s4_b0>.8h, V<w05_s3_b0>.8h, V<t07_s4_b0>.8h
    sub V<w07_s4_b0>.8h, V<w05_s3_b0>.8h, V<t07_s4_b0>.8h
    // work[0] with work[1], CT twiddle power 0.
    sqrdmulh V<q01_s5_b0>.8h, V<w01_s4_b0>.8h, V<pre_s345_b0_stage45_low>.h[0]
    mul V<t01_s5_b0>.8h, V<w01_s4_b0>.8h, V<tw_s345_b0_stage45_low>.h[0]
    mls V<t01_s5_b0>.8h, V<q01_s5_b0>.8h, v0.h[0]
    add V<w00_s5_b0>.8h, V<w00_s4_b0>.8h, V<t01_s5_b0>.8h
    sub V<w01_s5_b0>.8h, V<w00_s4_b0>.8h, V<t01_s5_b0>.8h
    // work[2] with work[3], CT twiddle power 8.
    sqrdmulh V<q03_s5_b0>.8h, V<w03_s4_b0>.8h, V<pre_s345_b0_stage45_low>.h[1]
    mul V<t03_s5_b0>.8h, V<w03_s4_b0>.8h, V<tw_s345_b0_stage45_low>.h[1]
    mls V<t03_s5_b0>.8h, V<q03_s5_b0>.8h, v0.h[0]
    add V<w02_s5_b0>.8h, V<w02_s4_b0>.8h, V<t03_s5_b0>.8h
    sub V<w03_s5_b0>.8h, V<w02_s4_b0>.8h, V<t03_s5_b0>.8h
    // work[4] with work[5], CT twiddle power 4.
    sqrdmulh V<q05_s5_b0>.8h, V<w05_s4_b0>.8h, V<pre_s345_b0_stage45_low>.h[2]
    mul V<t05_s5_b0>.8h, V<w05_s4_b0>.8h, V<tw_s345_b0_stage45_low>.h[2]
    mls V<t05_s5_b0>.8h, V<q05_s5_b0>.8h, v0.h[0]
    add V<w04_s5_b0>.8h, V<w04_s4_b0>.8h, V<t05_s5_b0>.8h
    sub V<w05_s5_b0>.8h, V<w04_s4_b0>.8h, V<t05_s5_b0>.8h
    // work[6] with work[7], CT twiddle power 12.
    sqrdmulh V<q07_s5_b0>.8h, V<w07_s4_b0>.8h, V<pre_s345_b0_stage45_low>.h[3]
    mul V<t07_s5_b0>.8h, V<w07_s4_b0>.8h, V<tw_s345_b0_stage45_low>.h[3]
    mls V<t07_s5_b0>.8h, V<q07_s5_b0>.8h, v0.h[0]
    add V<w06_s5_b0>.8h, V<w06_s4_b0>.8h, V<t07_s5_b0>.8h
    sub V<w07_s5_b0>.8h, V<w06_s4_b0>.8h, V<t07_s5_b0>.8h
    sqdmulh V<red00_s5_b0>.8h, V<w00_s5_b0>.8h, v0.h[1]
    srshr V<red00_s5_b0>.8h, V<red00_s5_b0>.8h, #11
    mls V<w00_s5_b0>.8h, V<red00_s5_b0>.8h, v0.h[0]
    str Q<w00_s5_b0>, [row_base, #16*0]
    sqdmulh V<red01_s5_b0>.8h, V<w01_s5_b0>.8h, v0.h[1]
    srshr V<red01_s5_b0>.8h, V<red01_s5_b0>.8h, #11
    mls V<w01_s5_b0>.8h, V<red01_s5_b0>.8h, v0.h[0]
    str Q<w01_s5_b0>, [row_base, #16*1]
    sqdmulh V<red02_s5_b0>.8h, V<w02_s5_b0>.8h, v0.h[1]
    srshr V<red02_s5_b0>.8h, V<red02_s5_b0>.8h, #11
    mls V<w02_s5_b0>.8h, V<red02_s5_b0>.8h, v0.h[0]
    str Q<w02_s5_b0>, [row_base, #16*2]
    sqdmulh V<red03_s5_b0>.8h, V<w03_s5_b0>.8h, v0.h[1]
    srshr V<red03_s5_b0>.8h, V<red03_s5_b0>.8h, #11
    mls V<w03_s5_b0>.8h, V<red03_s5_b0>.8h, v0.h[0]
    str Q<w03_s5_b0>, [row_base, #16*3]
    sqdmulh V<red04_s5_b0>.8h, V<w04_s5_b0>.8h, v0.h[1]
    srshr V<red04_s5_b0>.8h, V<red04_s5_b0>.8h, #11
    mls V<w04_s5_b0>.8h, V<red04_s5_b0>.8h, v0.h[0]
    str Q<w04_s5_b0>, [row_base, #16*4]
    sqdmulh V<red05_s5_b0>.8h, V<w05_s5_b0>.8h, v0.h[1]
    srshr V<red05_s5_b0>.8h, V<red05_s5_b0>.8h, #11
    mls V<w05_s5_b0>.8h, V<red05_s5_b0>.8h, v0.h[0]
    str Q<w05_s5_b0>, [row_base, #16*5]
    sqdmulh V<red06_s5_b0>.8h, V<w06_s5_b0>.8h, v0.h[1]
    srshr V<red06_s5_b0>.8h, V<red06_s5_b0>.8h, #11
    mls V<w06_s5_b0>.8h, V<red06_s5_b0>.8h, v0.h[0]
    str Q<w06_s5_b0>, [row_base, #16*6]
    sqdmulh V<red07_s5_b0>.8h, V<w07_s5_b0>.8h, v0.h[1]
    srshr V<red07_s5_b0>.8h, V<red07_s5_b0>.8h, #11
    mls V<w07_s5_b0>.8h, V<red07_s5_b0>.8h, v0.h[0]
    str Q<w07_s5_b0>, [row_base, #16*7]
_ntt32_stage345_block0_slothy_end:

    adr tw_ptr, ntt32_twiddle_stage3
_ntt32_stage345_block1_slothy_start:
    // Stage 3/4/5 block 1: work[8..15].
    ldr Q<w08_s2_b1>, [row_base, #16*8]
    ldr Q<w09_s2_b1>, [row_base, #16*9]
    ldr Q<w10_s2_b1>, [row_base, #16*10]
    ldr Q<w11_s2_b1>, [row_base, #16*11]
    ldr Q<w12_s2_b1>, [row_base, #16*12]
    ldr Q<w13_s2_b1>, [row_base, #16*13]
    ldr Q<w14_s2_b1>, [row_base, #16*14]
    ldr Q<w15_s2_b1>, [row_base, #16*15]

    ldr Q<tw_s345_b1_stage3>,  [tw_ptr], #16
    ldr Q<pre_s345_b1_stage3>, [tw_ptr], #16

    // work[8] with work[12], CT twiddle power 8.
    sqrdmulh V<q12_s3_b1>.8h, V<w12_s2_b1>.8h, V<pre_s345_b1_stage3>.h[1]
    mul V<t12_s3_b1>.8h, V<w12_s2_b1>.8h, V<tw_s345_b1_stage3>.h[1]
    mls V<t12_s3_b1>.8h, V<q12_s3_b1>.8h, v0.h[0]
    add V<w08_s3_b1>.8h, V<w08_s2_b1>.8h, V<t12_s3_b1>.8h
    sub V<w12_s3_b1>.8h, V<w08_s2_b1>.8h, V<t12_s3_b1>.8h
    // work[9] with work[13], CT twiddle power 8.
    sqrdmulh V<q13_s3_b1>.8h, V<w13_s2_b1>.8h, V<pre_s345_b1_stage3>.h[1]
    mul V<t13_s3_b1>.8h, V<w13_s2_b1>.8h, V<tw_s345_b1_stage3>.h[1]
    mls V<t13_s3_b1>.8h, V<q13_s3_b1>.8h, v0.h[0]
    add V<w09_s3_b1>.8h, V<w09_s2_b1>.8h, V<t13_s3_b1>.8h
    sub V<w13_s3_b1>.8h, V<w09_s2_b1>.8h, V<t13_s3_b1>.8h
    // work[10] with work[14], CT twiddle power 8.
    sqrdmulh V<q14_s3_b1>.8h, V<w14_s2_b1>.8h, V<pre_s345_b1_stage3>.h[1]
    mul V<t14_s3_b1>.8h, V<w14_s2_b1>.8h, V<tw_s345_b1_stage3>.h[1]
    mls V<t14_s3_b1>.8h, V<q14_s3_b1>.8h, v0.h[0]
    add V<w10_s3_b1>.8h, V<w10_s2_b1>.8h, V<t14_s3_b1>.8h
    sub V<w14_s3_b1>.8h, V<w10_s2_b1>.8h, V<t14_s3_b1>.8h
    // work[11] with work[15], CT twiddle power 8.
    sqrdmulh V<q15_s3_b1>.8h, V<w15_s2_b1>.8h, V<pre_s345_b1_stage3>.h[1]
    mul V<t15_s3_b1>.8h, V<w15_s2_b1>.8h, V<tw_s345_b1_stage3>.h[1]
    mls V<t15_s3_b1>.8h, V<q15_s3_b1>.8h, v0.h[0]
    add V<w11_s3_b1>.8h, V<w11_s2_b1>.8h, V<t15_s3_b1>.8h
    sub V<w15_s3_b1>.8h, V<w11_s2_b1>.8h, V<t15_s3_b1>.8h
    ldr Q<tw_s345_b1_stage45_low>,  [tw_ptr], #16
    ldr Q<pre_s345_b1_stage45_low>, [tw_ptr], #16

    // work[8] with work[10], CT twiddle power 4.
    sqrdmulh V<q10_s4_b1>.8h, V<w10_s3_b1>.8h, V<pre_s345_b1_stage45_low>.h[2]
    mul V<t10_s4_b1>.8h, V<w10_s3_b1>.8h, V<tw_s345_b1_stage45_low>.h[2]
    mls V<t10_s4_b1>.8h, V<q10_s4_b1>.8h, v0.h[0]
    add V<w08_s4_b1>.8h, V<w08_s3_b1>.8h, V<t10_s4_b1>.8h
    sub V<w10_s4_b1>.8h, V<w08_s3_b1>.8h, V<t10_s4_b1>.8h
    // work[9] with work[11], CT twiddle power 4.
    sqrdmulh V<q11_s4_b1>.8h, V<w11_s3_b1>.8h, V<pre_s345_b1_stage45_low>.h[2]
    mul V<t11_s4_b1>.8h, V<w11_s3_b1>.8h, V<tw_s345_b1_stage45_low>.h[2]
    mls V<t11_s4_b1>.8h, V<q11_s4_b1>.8h, v0.h[0]
    add V<w09_s4_b1>.8h, V<w09_s3_b1>.8h, V<t11_s4_b1>.8h
    sub V<w11_s4_b1>.8h, V<w09_s3_b1>.8h, V<t11_s4_b1>.8h
    // work[12] with work[14], CT twiddle power 12.
    sqrdmulh V<q14_s4_b1>.8h, V<w14_s3_b1>.8h, V<pre_s345_b1_stage45_low>.h[3]
    mul V<t14_s4_b1>.8h, V<w14_s3_b1>.8h, V<tw_s345_b1_stage45_low>.h[3]
    mls V<t14_s4_b1>.8h, V<q14_s4_b1>.8h, v0.h[0]
    add V<w12_s4_b1>.8h, V<w12_s3_b1>.8h, V<t14_s4_b1>.8h
    sub V<w14_s4_b1>.8h, V<w12_s3_b1>.8h, V<t14_s4_b1>.8h
    // work[13] with work[15], CT twiddle power 12.
    sqrdmulh V<q15_s4_b1>.8h, V<w15_s3_b1>.8h, V<pre_s345_b1_stage45_low>.h[3]
    mul V<t15_s4_b1>.8h, V<w15_s3_b1>.8h, V<tw_s345_b1_stage45_low>.h[3]
    mls V<t15_s4_b1>.8h, V<q15_s4_b1>.8h, v0.h[0]
    add V<w13_s4_b1>.8h, V<w13_s3_b1>.8h, V<t15_s4_b1>.8h
    sub V<w15_s4_b1>.8h, V<w13_s3_b1>.8h, V<t15_s4_b1>.8h
    // work[8] with work[9], CT twiddle power 2.
    sqrdmulh V<q09_s5_b1>.8h, V<w09_s4_b1>.8h, V<pre_s345_b1_stage45_low>.h[4]
    mul V<t09_s5_b1>.8h, V<w09_s4_b1>.8h, V<tw_s345_b1_stage45_low>.h[4]
    mls V<t09_s5_b1>.8h, V<q09_s5_b1>.8h, v0.h[0]
    add V<w08_s5_b1>.8h, V<w08_s4_b1>.8h, V<t09_s5_b1>.8h
    sub V<w09_s5_b1>.8h, V<w08_s4_b1>.8h, V<t09_s5_b1>.8h
    // work[10] with work[11], CT twiddle power 10.
    sqrdmulh V<q11_s5_b1>.8h, V<w11_s4_b1>.8h, V<pre_s345_b1_stage45_low>.h[5]
    mul V<t11_s5_b1>.8h, V<w11_s4_b1>.8h, V<tw_s345_b1_stage45_low>.h[5]
    mls V<t11_s5_b1>.8h, V<q11_s5_b1>.8h, v0.h[0]
    add V<w10_s5_b1>.8h, V<w10_s4_b1>.8h, V<t11_s5_b1>.8h
    sub V<w11_s5_b1>.8h, V<w10_s4_b1>.8h, V<t11_s5_b1>.8h
    // work[12] with work[13], CT twiddle power 6.
    sqrdmulh V<q13_s5_b1>.8h, V<w13_s4_b1>.8h, V<pre_s345_b1_stage45_low>.h[6]
    mul V<t13_s5_b1>.8h, V<w13_s4_b1>.8h, V<tw_s345_b1_stage45_low>.h[6]
    mls V<t13_s5_b1>.8h, V<q13_s5_b1>.8h, v0.h[0]
    add V<w12_s5_b1>.8h, V<w12_s4_b1>.8h, V<t13_s5_b1>.8h
    sub V<w13_s5_b1>.8h, V<w12_s4_b1>.8h, V<t13_s5_b1>.8h
    // work[14] with work[15], CT twiddle power 14.
    sqrdmulh V<q15_s5_b1>.8h, V<w15_s4_b1>.8h, V<pre_s345_b1_stage45_low>.h[7]
    mul V<t15_s5_b1>.8h, V<w15_s4_b1>.8h, V<tw_s345_b1_stage45_low>.h[7]
    mls V<t15_s5_b1>.8h, V<q15_s5_b1>.8h, v0.h[0]
    add V<w14_s5_b1>.8h, V<w14_s4_b1>.8h, V<t15_s5_b1>.8h
    sub V<w15_s5_b1>.8h, V<w14_s4_b1>.8h, V<t15_s5_b1>.8h
    sqdmulh V<red08_s5_b1>.8h, V<w08_s5_b1>.8h, v0.h[1]
    srshr V<red08_s5_b1>.8h, V<red08_s5_b1>.8h, #11
    mls V<w08_s5_b1>.8h, V<red08_s5_b1>.8h, v0.h[0]
    str Q<w08_s5_b1>, [row_base, #16*8]
    sqdmulh V<red09_s5_b1>.8h, V<w09_s5_b1>.8h, v0.h[1]
    srshr V<red09_s5_b1>.8h, V<red09_s5_b1>.8h, #11
    mls V<w09_s5_b1>.8h, V<red09_s5_b1>.8h, v0.h[0]
    str Q<w09_s5_b1>, [row_base, #16*9]
    sqdmulh V<red10_s5_b1>.8h, V<w10_s5_b1>.8h, v0.h[1]
    srshr V<red10_s5_b1>.8h, V<red10_s5_b1>.8h, #11
    mls V<w10_s5_b1>.8h, V<red10_s5_b1>.8h, v0.h[0]
    str Q<w10_s5_b1>, [row_base, #16*10]
    sqdmulh V<red11_s5_b1>.8h, V<w11_s5_b1>.8h, v0.h[1]
    srshr V<red11_s5_b1>.8h, V<red11_s5_b1>.8h, #11
    mls V<w11_s5_b1>.8h, V<red11_s5_b1>.8h, v0.h[0]
    str Q<w11_s5_b1>, [row_base, #16*11]
    sqdmulh V<red12_s5_b1>.8h, V<w12_s5_b1>.8h, v0.h[1]
    srshr V<red12_s5_b1>.8h, V<red12_s5_b1>.8h, #11
    mls V<w12_s5_b1>.8h, V<red12_s5_b1>.8h, v0.h[0]
    str Q<w12_s5_b1>, [row_base, #16*12]
    sqdmulh V<red13_s5_b1>.8h, V<w13_s5_b1>.8h, v0.h[1]
    srshr V<red13_s5_b1>.8h, V<red13_s5_b1>.8h, #11
    mls V<w13_s5_b1>.8h, V<red13_s5_b1>.8h, v0.h[0]
    str Q<w13_s5_b1>, [row_base, #16*13]
    sqdmulh V<red14_s5_b1>.8h, V<w14_s5_b1>.8h, v0.h[1]
    srshr V<red14_s5_b1>.8h, V<red14_s5_b1>.8h, #11
    mls V<w14_s5_b1>.8h, V<red14_s5_b1>.8h, v0.h[0]
    str Q<w14_s5_b1>, [row_base, #16*14]
    sqdmulh V<red15_s5_b1>.8h, V<w15_s5_b1>.8h, v0.h[1]
    srshr V<red15_s5_b1>.8h, V<red15_s5_b1>.8h, #11
    mls V<w15_s5_b1>.8h, V<red15_s5_b1>.8h, v0.h[0]
    str Q<w15_s5_b1>, [row_base, #16*15]
_ntt32_stage345_block1_slothy_end:

    adr tw_ptr, ntt32_twiddle_stage3
_ntt32_stage345_block2_slothy_start:
    // Stage 3/4/5 block 2: work[16..23].
    ldr Q<w16_s2_b2>, [row_base, #16*16]
    ldr Q<w17_s2_b2>, [row_base, #16*17]
    ldr Q<w18_s2_b2>, [row_base, #16*18]
    ldr Q<w19_s2_b2>, [row_base, #16*19]
    ldr Q<w20_s2_b2>, [row_base, #16*20]
    ldr Q<w21_s2_b2>, [row_base, #16*21]
    ldr Q<w22_s2_b2>, [row_base, #16*22]
    ldr Q<w23_s2_b2>, [row_base, #16*23]

    ldr Q<tw_s345_b2_stage3>,  [tw_ptr], #16
    ldr Q<pre_s345_b2_stage3>, [tw_ptr], #16

    // work[16] with work[20], CT twiddle power 4.
    sqrdmulh V<q20_s3_b2>.8h, V<w20_s2_b2>.8h, V<pre_s345_b2_stage3>.h[2]
    mul V<t20_s3_b2>.8h, V<w20_s2_b2>.8h, V<tw_s345_b2_stage3>.h[2]
    mls V<t20_s3_b2>.8h, V<q20_s3_b2>.8h, v0.h[0]
    add V<w16_s3_b2>.8h, V<w16_s2_b2>.8h, V<t20_s3_b2>.8h
    sub V<w20_s3_b2>.8h, V<w16_s2_b2>.8h, V<t20_s3_b2>.8h
    // work[17] with work[21], CT twiddle power 4.
    sqrdmulh V<q21_s3_b2>.8h, V<w21_s2_b2>.8h, V<pre_s345_b2_stage3>.h[2]
    mul V<t21_s3_b2>.8h, V<w21_s2_b2>.8h, V<tw_s345_b2_stage3>.h[2]
    mls V<t21_s3_b2>.8h, V<q21_s3_b2>.8h, v0.h[0]
    add V<w17_s3_b2>.8h, V<w17_s2_b2>.8h, V<t21_s3_b2>.8h
    sub V<w21_s3_b2>.8h, V<w17_s2_b2>.8h, V<t21_s3_b2>.8h
    // work[18] with work[22], CT twiddle power 4.
    sqrdmulh V<q22_s3_b2>.8h, V<w22_s2_b2>.8h, V<pre_s345_b2_stage3>.h[2]
    mul V<t22_s3_b2>.8h, V<w22_s2_b2>.8h, V<tw_s345_b2_stage3>.h[2]
    mls V<t22_s3_b2>.8h, V<q22_s3_b2>.8h, v0.h[0]
    add V<w18_s3_b2>.8h, V<w18_s2_b2>.8h, V<t22_s3_b2>.8h
    sub V<w22_s3_b2>.8h, V<w18_s2_b2>.8h, V<t22_s3_b2>.8h
    // work[19] with work[23], CT twiddle power 4.
    sqrdmulh V<q23_s3_b2>.8h, V<w23_s2_b2>.8h, V<pre_s345_b2_stage3>.h[2]
    mul V<t23_s3_b2>.8h, V<w23_s2_b2>.8h, V<tw_s345_b2_stage3>.h[2]
    mls V<t23_s3_b2>.8h, V<q23_s3_b2>.8h, v0.h[0]
    add V<w19_s3_b2>.8h, V<w19_s2_b2>.8h, V<t23_s3_b2>.8h
    sub V<w23_s3_b2>.8h, V<w19_s2_b2>.8h, V<t23_s3_b2>.8h
    ldr Q<tw_s345_b2_stage45_low>,  [tw_ptr], #16
    ldr Q<pre_s345_b2_stage45_low>, [tw_ptr], #16

    // work[16] with work[18], CT twiddle power 2.
    sqrdmulh V<q18_s4_b2>.8h, V<w18_s3_b2>.8h, V<pre_s345_b2_stage45_low>.h[4]
    mul V<t18_s4_b2>.8h, V<w18_s3_b2>.8h, V<tw_s345_b2_stage45_low>.h[4]
    mls V<t18_s4_b2>.8h, V<q18_s4_b2>.8h, v0.h[0]
    add V<w16_s4_b2>.8h, V<w16_s3_b2>.8h, V<t18_s4_b2>.8h
    sub V<w18_s4_b2>.8h, V<w16_s3_b2>.8h, V<t18_s4_b2>.8h
    // work[17] with work[19], CT twiddle power 2.
    sqrdmulh V<q19_s4_b2>.8h, V<w19_s3_b2>.8h, V<pre_s345_b2_stage45_low>.h[4]
    mul V<t19_s4_b2>.8h, V<w19_s3_b2>.8h, V<tw_s345_b2_stage45_low>.h[4]
    mls V<t19_s4_b2>.8h, V<q19_s4_b2>.8h, v0.h[0]
    add V<w17_s4_b2>.8h, V<w17_s3_b2>.8h, V<t19_s4_b2>.8h
    sub V<w19_s4_b2>.8h, V<w17_s3_b2>.8h, V<t19_s4_b2>.8h
    // work[20] with work[22], CT twiddle power 10.
    sqrdmulh V<q22_s4_b2>.8h, V<w22_s3_b2>.8h, V<pre_s345_b2_stage45_low>.h[5]
    mul V<t22_s4_b2>.8h, V<w22_s3_b2>.8h, V<tw_s345_b2_stage45_low>.h[5]
    mls V<t22_s4_b2>.8h, V<q22_s4_b2>.8h, v0.h[0]
    add V<w20_s4_b2>.8h, V<w20_s3_b2>.8h, V<t22_s4_b2>.8h
    sub V<w22_s4_b2>.8h, V<w20_s3_b2>.8h, V<t22_s4_b2>.8h
    // work[21] with work[23], CT twiddle power 10.
    sqrdmulh V<q23_s4_b2>.8h, V<w23_s3_b2>.8h, V<pre_s345_b2_stage45_low>.h[5]
    mul V<t23_s4_b2>.8h, V<w23_s3_b2>.8h, V<tw_s345_b2_stage45_low>.h[5]
    mls V<t23_s4_b2>.8h, V<q23_s4_b2>.8h, v0.h[0]
    add V<w21_s4_b2>.8h, V<w21_s3_b2>.8h, V<t23_s4_b2>.8h
    sub V<w23_s4_b2>.8h, V<w21_s3_b2>.8h, V<t23_s4_b2>.8h
    ldr Q<tw_s345_b2_stage5_high>,  [tw_ptr], #16
    ldr Q<pre_s345_b2_stage5_high>, [tw_ptr], #16

    // work[16] with work[17], CT twiddle power 1.
    sqrdmulh V<q17_s5_b2>.8h, V<w17_s4_b2>.8h, V<pre_s345_b2_stage5_high>.h[0]
    mul V<t17_s5_b2>.8h, V<w17_s4_b2>.8h, V<tw_s345_b2_stage5_high>.h[0]
    mls V<t17_s5_b2>.8h, V<q17_s5_b2>.8h, v0.h[0]
    add V<w16_s5_b2>.8h, V<w16_s4_b2>.8h, V<t17_s5_b2>.8h
    sub V<w17_s5_b2>.8h, V<w16_s4_b2>.8h, V<t17_s5_b2>.8h
    // work[18] with work[19], CT twiddle power 9.
    sqrdmulh V<q19_s5_b2>.8h, V<w19_s4_b2>.8h, V<pre_s345_b2_stage5_high>.h[1]
    mul V<t19_s5_b2>.8h, V<w19_s4_b2>.8h, V<tw_s345_b2_stage5_high>.h[1]
    mls V<t19_s5_b2>.8h, V<q19_s5_b2>.8h, v0.h[0]
    add V<w18_s5_b2>.8h, V<w18_s4_b2>.8h, V<t19_s5_b2>.8h
    sub V<w19_s5_b2>.8h, V<w18_s4_b2>.8h, V<t19_s5_b2>.8h
    // work[20] with work[21], CT twiddle power 5.
    sqrdmulh V<q21_s5_b2>.8h, V<w21_s4_b2>.8h, V<pre_s345_b2_stage5_high>.h[2]
    mul V<t21_s5_b2>.8h, V<w21_s4_b2>.8h, V<tw_s345_b2_stage5_high>.h[2]
    mls V<t21_s5_b2>.8h, V<q21_s5_b2>.8h, v0.h[0]
    add V<w20_s5_b2>.8h, V<w20_s4_b2>.8h, V<t21_s5_b2>.8h
    sub V<w21_s5_b2>.8h, V<w20_s4_b2>.8h, V<t21_s5_b2>.8h
    // work[22] with work[23], CT twiddle power 13.
    sqrdmulh V<q23_s5_b2>.8h, V<w23_s4_b2>.8h, V<pre_s345_b2_stage5_high>.h[3]
    mul V<t23_s5_b2>.8h, V<w23_s4_b2>.8h, V<tw_s345_b2_stage5_high>.h[3]
    mls V<t23_s5_b2>.8h, V<q23_s5_b2>.8h, v0.h[0]
    add V<w22_s5_b2>.8h, V<w22_s4_b2>.8h, V<t23_s5_b2>.8h
    sub V<w23_s5_b2>.8h, V<w22_s4_b2>.8h, V<t23_s5_b2>.8h
    sqdmulh V<red16_s5_b2>.8h, V<w16_s5_b2>.8h, v0.h[1]
    srshr V<red16_s5_b2>.8h, V<red16_s5_b2>.8h, #11
    mls V<w16_s5_b2>.8h, V<red16_s5_b2>.8h, v0.h[0]
    str Q<w16_s5_b2>, [row_base, #16*16]
    sqdmulh V<red17_s5_b2>.8h, V<w17_s5_b2>.8h, v0.h[1]
    srshr V<red17_s5_b2>.8h, V<red17_s5_b2>.8h, #11
    mls V<w17_s5_b2>.8h, V<red17_s5_b2>.8h, v0.h[0]
    str Q<w17_s5_b2>, [row_base, #16*17]
    sqdmulh V<red18_s5_b2>.8h, V<w18_s5_b2>.8h, v0.h[1]
    srshr V<red18_s5_b2>.8h, V<red18_s5_b2>.8h, #11
    mls V<w18_s5_b2>.8h, V<red18_s5_b2>.8h, v0.h[0]
    str Q<w18_s5_b2>, [row_base, #16*18]
    sqdmulh V<red19_s5_b2>.8h, V<w19_s5_b2>.8h, v0.h[1]
    srshr V<red19_s5_b2>.8h, V<red19_s5_b2>.8h, #11
    mls V<w19_s5_b2>.8h, V<red19_s5_b2>.8h, v0.h[0]
    str Q<w19_s5_b2>, [row_base, #16*19]
    sqdmulh V<red20_s5_b2>.8h, V<w20_s5_b2>.8h, v0.h[1]
    srshr V<red20_s5_b2>.8h, V<red20_s5_b2>.8h, #11
    mls V<w20_s5_b2>.8h, V<red20_s5_b2>.8h, v0.h[0]
    str Q<w20_s5_b2>, [row_base, #16*20]
    sqdmulh V<red21_s5_b2>.8h, V<w21_s5_b2>.8h, v0.h[1]
    srshr V<red21_s5_b2>.8h, V<red21_s5_b2>.8h, #11
    mls V<w21_s5_b2>.8h, V<red21_s5_b2>.8h, v0.h[0]
    str Q<w21_s5_b2>, [row_base, #16*21]
    sqdmulh V<red22_s5_b2>.8h, V<w22_s5_b2>.8h, v0.h[1]
    srshr V<red22_s5_b2>.8h, V<red22_s5_b2>.8h, #11
    mls V<w22_s5_b2>.8h, V<red22_s5_b2>.8h, v0.h[0]
    str Q<w22_s5_b2>, [row_base, #16*22]
    sqdmulh V<red23_s5_b2>.8h, V<w23_s5_b2>.8h, v0.h[1]
    srshr V<red23_s5_b2>.8h, V<red23_s5_b2>.8h, #11
    mls V<w23_s5_b2>.8h, V<red23_s5_b2>.8h, v0.h[0]
    str Q<w23_s5_b2>, [row_base, #16*23]
_ntt32_stage345_block2_slothy_end:

    adr tw_ptr, ntt32_twiddle_stage3
_ntt32_stage345_block3_slothy_start:
    // Stage 3/4/5 block 3: work[24..31].
    ldr Q<w24_s2_b3>, [row_base, #16*24]
    ldr Q<w25_s2_b3>, [row_base, #16*25]
    ldr Q<w26_s2_b3>, [row_base, #16*26]
    ldr Q<w27_s2_b3>, [row_base, #16*27]
    ldr Q<w28_s2_b3>, [row_base, #16*28]
    ldr Q<w29_s2_b3>, [row_base, #16*29]
    ldr Q<w30_s2_b3>, [row_base, #16*30]
    ldr Q<w31_s2_b3>, [row_base, #16*31]

    ldr Q<tw_s345_b3_stage3>,  [tw_ptr], #16
    ldr Q<pre_s345_b3_stage3>, [tw_ptr], #16

    // work[24] with work[28], CT twiddle power 12.
    sqrdmulh V<q28_s3_b3>.8h, V<w28_s2_b3>.8h, V<pre_s345_b3_stage3>.h[3]
    mul V<t28_s3_b3>.8h, V<w28_s2_b3>.8h, V<tw_s345_b3_stage3>.h[3]
    mls V<t28_s3_b3>.8h, V<q28_s3_b3>.8h, v0.h[0]
    add V<w24_s3_b3>.8h, V<w24_s2_b3>.8h, V<t28_s3_b3>.8h
    sub V<w28_s3_b3>.8h, V<w24_s2_b3>.8h, V<t28_s3_b3>.8h
    // work[25] with work[29], CT twiddle power 12.
    sqrdmulh V<q29_s3_b3>.8h, V<w29_s2_b3>.8h, V<pre_s345_b3_stage3>.h[3]
    mul V<t29_s3_b3>.8h, V<w29_s2_b3>.8h, V<tw_s345_b3_stage3>.h[3]
    mls V<t29_s3_b3>.8h, V<q29_s3_b3>.8h, v0.h[0]
    add V<w25_s3_b3>.8h, V<w25_s2_b3>.8h, V<t29_s3_b3>.8h
    sub V<w29_s3_b3>.8h, V<w25_s2_b3>.8h, V<t29_s3_b3>.8h
    // work[26] with work[30], CT twiddle power 12.
    sqrdmulh V<q30_s3_b3>.8h, V<w30_s2_b3>.8h, V<pre_s345_b3_stage3>.h[3]
    mul V<t30_s3_b3>.8h, V<w30_s2_b3>.8h, V<tw_s345_b3_stage3>.h[3]
    mls V<t30_s3_b3>.8h, V<q30_s3_b3>.8h, v0.h[0]
    add V<w26_s3_b3>.8h, V<w26_s2_b3>.8h, V<t30_s3_b3>.8h
    sub V<w30_s3_b3>.8h, V<w26_s2_b3>.8h, V<t30_s3_b3>.8h
    // work[27] with work[31], CT twiddle power 12.
    sqrdmulh V<q31_s3_b3>.8h, V<w31_s2_b3>.8h, V<pre_s345_b3_stage3>.h[3]
    mul V<t31_s3_b3>.8h, V<w31_s2_b3>.8h, V<tw_s345_b3_stage3>.h[3]
    mls V<t31_s3_b3>.8h, V<q31_s3_b3>.8h, v0.h[0]
    add V<w27_s3_b3>.8h, V<w27_s2_b3>.8h, V<t31_s3_b3>.8h
    sub V<w31_s3_b3>.8h, V<w27_s2_b3>.8h, V<t31_s3_b3>.8h
    ldr Q<tw_s345_b3_stage45_low>,  [tw_ptr], #16
    ldr Q<pre_s345_b3_stage45_low>, [tw_ptr], #16

    // work[24] with work[26], CT twiddle power 6.
    sqrdmulh V<q26_s4_b3>.8h, V<w26_s3_b3>.8h, V<pre_s345_b3_stage45_low>.h[6]
    mul V<t26_s4_b3>.8h, V<w26_s3_b3>.8h, V<tw_s345_b3_stage45_low>.h[6]
    mls V<t26_s4_b3>.8h, V<q26_s4_b3>.8h, v0.h[0]
    add V<w24_s4_b3>.8h, V<w24_s3_b3>.8h, V<t26_s4_b3>.8h
    sub V<w26_s4_b3>.8h, V<w24_s3_b3>.8h, V<t26_s4_b3>.8h
    // work[25] with work[27], CT twiddle power 6.
    sqrdmulh V<q27_s4_b3>.8h, V<w27_s3_b3>.8h, V<pre_s345_b3_stage45_low>.h[6]
    mul V<t27_s4_b3>.8h, V<w27_s3_b3>.8h, V<tw_s345_b3_stage45_low>.h[6]
    mls V<t27_s4_b3>.8h, V<q27_s4_b3>.8h, v0.h[0]
    add V<w25_s4_b3>.8h, V<w25_s3_b3>.8h, V<t27_s4_b3>.8h
    sub V<w27_s4_b3>.8h, V<w25_s3_b3>.8h, V<t27_s4_b3>.8h
    // work[28] with work[30], CT twiddle power 14.
    sqrdmulh V<q30_s4_b3>.8h, V<w30_s3_b3>.8h, V<pre_s345_b3_stage45_low>.h[7]
    mul V<t30_s4_b3>.8h, V<w30_s3_b3>.8h, V<tw_s345_b3_stage45_low>.h[7]
    mls V<t30_s4_b3>.8h, V<q30_s4_b3>.8h, v0.h[0]
    add V<w28_s4_b3>.8h, V<w28_s3_b3>.8h, V<t30_s4_b3>.8h
    sub V<w30_s4_b3>.8h, V<w28_s3_b3>.8h, V<t30_s4_b3>.8h
    // work[29] with work[31], CT twiddle power 14.
    sqrdmulh V<q31_s4_b3>.8h, V<w31_s3_b3>.8h, V<pre_s345_b3_stage45_low>.h[7]
    mul V<t31_s4_b3>.8h, V<w31_s3_b3>.8h, V<tw_s345_b3_stage45_low>.h[7]
    mls V<t31_s4_b3>.8h, V<q31_s4_b3>.8h, v0.h[0]
    add V<w29_s4_b3>.8h, V<w29_s3_b3>.8h, V<t31_s4_b3>.8h
    sub V<w31_s4_b3>.8h, V<w29_s3_b3>.8h, V<t31_s4_b3>.8h
    ldr Q<tw_s345_b3_stage5_high>,  [tw_ptr], #16
    ldr Q<pre_s345_b3_stage5_high>, [tw_ptr], #16

    // work[24] with work[25], CT twiddle power 3.
    sqrdmulh V<q25_s5_b3>.8h, V<w25_s4_b3>.8h, V<pre_s345_b3_stage5_high>.h[4]
    mul V<t25_s5_b3>.8h, V<w25_s4_b3>.8h, V<tw_s345_b3_stage5_high>.h[4]
    mls V<t25_s5_b3>.8h, V<q25_s5_b3>.8h, v0.h[0]
    add V<w24_s5_b3>.8h, V<w24_s4_b3>.8h, V<t25_s5_b3>.8h
    sub V<w25_s5_b3>.8h, V<w24_s4_b3>.8h, V<t25_s5_b3>.8h
    // work[26] with work[27], CT twiddle power 11.
    sqrdmulh V<q27_s5_b3>.8h, V<w27_s4_b3>.8h, V<pre_s345_b3_stage5_high>.h[5]
    mul V<t27_s5_b3>.8h, V<w27_s4_b3>.8h, V<tw_s345_b3_stage5_high>.h[5]
    mls V<t27_s5_b3>.8h, V<q27_s5_b3>.8h, v0.h[0]
    add V<w26_s5_b3>.8h, V<w26_s4_b3>.8h, V<t27_s5_b3>.8h
    sub V<w27_s5_b3>.8h, V<w26_s4_b3>.8h, V<t27_s5_b3>.8h
    // work[28] with work[29], CT twiddle power 7.
    sqrdmulh V<q29_s5_b3>.8h, V<w29_s4_b3>.8h, V<pre_s345_b3_stage5_high>.h[6]
    mul V<t29_s5_b3>.8h, V<w29_s4_b3>.8h, V<tw_s345_b3_stage5_high>.h[6]
    mls V<t29_s5_b3>.8h, V<q29_s5_b3>.8h, v0.h[0]
    add V<w28_s5_b3>.8h, V<w28_s4_b3>.8h, V<t29_s5_b3>.8h
    sub V<w29_s5_b3>.8h, V<w28_s4_b3>.8h, V<t29_s5_b3>.8h
    // work[30] with work[31], CT twiddle power 15.
    sqrdmulh V<q31_s5_b3>.8h, V<w31_s4_b3>.8h, V<pre_s345_b3_stage5_high>.h[7]
    mul V<t31_s5_b3>.8h, V<w31_s4_b3>.8h, V<tw_s345_b3_stage5_high>.h[7]
    mls V<t31_s5_b3>.8h, V<q31_s5_b3>.8h, v0.h[0]
    add V<w30_s5_b3>.8h, V<w30_s4_b3>.8h, V<t31_s5_b3>.8h
    sub V<w31_s5_b3>.8h, V<w30_s4_b3>.8h, V<t31_s5_b3>.8h
    sqdmulh V<red24_s5_b3>.8h, V<w24_s5_b3>.8h, v0.h[1]
    srshr V<red24_s5_b3>.8h, V<red24_s5_b3>.8h, #11
    mls V<w24_s5_b3>.8h, V<red24_s5_b3>.8h, v0.h[0]
    str Q<w24_s5_b3>, [row_base, #16*24]
    sqdmulh V<red25_s5_b3>.8h, V<w25_s5_b3>.8h, v0.h[1]
    srshr V<red25_s5_b3>.8h, V<red25_s5_b3>.8h, #11
    mls V<w25_s5_b3>.8h, V<red25_s5_b3>.8h, v0.h[0]
    str Q<w25_s5_b3>, [row_base, #16*25]
    sqdmulh V<red26_s5_b3>.8h, V<w26_s5_b3>.8h, v0.h[1]
    srshr V<red26_s5_b3>.8h, V<red26_s5_b3>.8h, #11
    mls V<w26_s5_b3>.8h, V<red26_s5_b3>.8h, v0.h[0]
    str Q<w26_s5_b3>, [row_base, #16*26]
    sqdmulh V<red27_s5_b3>.8h, V<w27_s5_b3>.8h, v0.h[1]
    srshr V<red27_s5_b3>.8h, V<red27_s5_b3>.8h, #11
    mls V<w27_s5_b3>.8h, V<red27_s5_b3>.8h, v0.h[0]
    str Q<w27_s5_b3>, [row_base, #16*27]
    sqdmulh V<red28_s5_b3>.8h, V<w28_s5_b3>.8h, v0.h[1]
    srshr V<red28_s5_b3>.8h, V<red28_s5_b3>.8h, #11
    mls V<w28_s5_b3>.8h, V<red28_s5_b3>.8h, v0.h[0]
    str Q<w28_s5_b3>, [row_base, #16*28]
    sqdmulh V<red29_s5_b3>.8h, V<w29_s5_b3>.8h, v0.h[1]
    srshr V<red29_s5_b3>.8h, V<red29_s5_b3>.8h, #11
    mls V<w29_s5_b3>.8h, V<red29_s5_b3>.8h, v0.h[0]
    str Q<w29_s5_b3>, [row_base, #16*29]
    sqdmulh V<red30_s5_b3>.8h, V<w30_s5_b3>.8h, v0.h[1]
    srshr V<red30_s5_b3>.8h, V<red30_s5_b3>.8h, #11
    mls V<w30_s5_b3>.8h, V<red30_s5_b3>.8h, v0.h[0]
    str Q<w30_s5_b3>, [row_base, #16*30]
    sqdmulh V<red31_s5_b3>.8h, V<w31_s5_b3>.8h, v0.h[1]
    srshr V<red31_s5_b3>.8h, V<red31_s5_b3>.8h, #11
    mls V<w31_s5_b3>.8h, V<red31_s5_b3>.8h, v0.h[0]
    str Q<w31_s5_b3>, [row_base, #16*31]
_ntt32_stage345_block3_slothy_end:

    ret

.align 4
ntt32_twiddle_vecs:
    // Powers 0..7, centered normal multipliers.
ntt32_twiddle_p00_07:
    .hword      1,  -1673,  -1241,  -1464,   1716,  -1558,    -44,   1015
    .hword      9, -15858, -11763, -13877,  16266, -14768,   -417,   9621

    // Powers 8..15, centered normal multipliers.
ntt32_twiddle_p08_15:
    .hword   -708,  -1267,    550,   -588,  -1521,    281,     39,    436
    .hword  -6711, -12010,   5213,  -5573, -14417,   2664,    370,   4133

    // CT stage 3 powers: 0,8,4,12.
ntt32_twiddle_stage3:
    .hword      1,   -708,   1716,  -1521,      0,      0,      0,      0
    .hword      9,  -6711,  16266, -14417,      0,      0,      0,      0

    // CT stage 4 powers and CT stage 5 low half: 0,8,4,12,2,10,6,14.
ntt32_twiddle_stage45_low:
    .hword      1,   -708,   1716,  -1521,  -1241,    550,    -44,     39
    .hword      9,  -6711,  16266, -14417, -11763,   5213,   -417,    370

    // CT stage 5 high half powers: 1,9,5,13,3,11,7,15.
ntt32_twiddle_stage5_high:
    .hword  -1673,  -1267,  -1558,    281,  -1464,   -588,   1015,    436
    .hword -15858, -12010, -14768,   2664, -13877,  -5573,   9621,   4133

    .unreq row_base
    .unreq tw_ptr

