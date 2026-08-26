#
# Copyright (c) 2026 Slothy contributors
# SPDX-License-Identifier: MIT
#

"""Cortex-A76 microarchitecture model.

The latency, throughput, pipeline, dispatch, and special-case data in this
module is derived from version 10.0 of the Arm Cortex-A76 Software Optimization
Guide (PJDOC-466751330-7215).

Unlike the older Neoverse-N1 approximation, this model represents heterogeneous
multi-uop resource use. In particular, ASIMD structure loads and ASIMD stores
use load/store and vector resources at the same time. It also models the
Cortex-A76's 4-Mop/8-uop dispatch limits, vector load forwarding latency,
partial-register dispatch repair, and recommended 32-byte loop alignment.

Multi-cycle resource durations are a scheduling abstraction calibrated to the
guide's documented steady-state throughput. They describe pipeline pressure;
they do not imply that one physical uop remains resident for the full duration.

Signed AdvSIMD SMULL/SMLAL timing is additionally validated on a Raspberry
Pi 5 Cortex-A76 r4p1 by two independent benchmark paths:
https://github.com/cheng-wei-huang0612/cortex-a76-smull-forwarding

Those measurements establish V0 execution, reciprocal throughput 1, ordinary
latency 4, SMULL-to-SMLAL accumulator latency 3, and SMLAL accumulator
self-forwarding latency 1 for both 16x16-to-32 and 32x32-to-64 forms.

The AArch64 architecture model itself is intentionally partial, so this target
model covers the instruction classes currently understood by
``aarch64_neon.py`` rather than the complete Armv8.2-A ISA.
"""

from enum import Enum
from itertools import product

from slothy.helper import lookup_multidict
from slothy.targets.aarch64.aarch64_neon import (
    RegisterType,
    find_class,
    all_subclass_leaves,
    Ldp_X,
    Ldp_W,
    Ldr_X,
    Str_X,
    Stp_X,
    Stp_W,
    Ldr_D,
    Ldr_Q,
    Str_Q,
    Stp_Q,
    Ldp_Q,
    Vrev,
    uaddlp,
    vmov,
    vmovi,
    vadd,
    vxtn,
    vusra,
    vdup,
    vdup_w,
    AESInstruction,
    Transpose,
    AArch64NeonLogical,
    VShiftImmediateBasic,
    VShiftRegBasic,
    AArch64BasicArithmetic,
    AArch64ConditionalSelect,
    AArch64ConditionalCompare,
    AArch64Logical,
    AArch64LogicalShifted,
    AArch64Move,
    AArch64Shift,
    Tst,
    AArch64ShiftedArithmetic,
    AArch64HighMultiply,
    AArch64Multiply,
    AArch64CRC32,
    VecToGprMov,
    q_st1_4_with_postinc,
    St3,
    St4,
    Vzip,
    vsub,
    vsrshr,
    Vmul,
    Vmla,
    Vqdmulh,
    Vmull,
    Vmlal,
    vsmull,
    vsmull2,
    vsmull_lane,
    vsmull2_lane,
    vsmlal,
    vsmlal2,
    vsmlal_lane,
    vsmlal2_lane,
    umull_wform,
    vmul,
    vmul_lane,
    vmla,
    vmla_lane,
    vmls,
    vmls_lane,
    AArch64NeonCount,
    ASimdCompare,
    vext,
    AArch64NeonShiftInsert,
    vtbl,
    vtbl_2,
    vuaddlv_sform,
    fmov_s_form,
    fmov_d_form,
    fmov_0,
    fmov_0_force_output,
    fmov_1,
    fmov_1_force_output,
    q_ldr1_stack,
    Q_Ld2_Lane_Post_Inc,
    q_ld2_lane_s,
    Ld4,
    Ld3,
    Ld2,
    q_ld1_2,
    St2,
    mov_wtov_s,
    mov_vtov_d,
    lsr,
    movk_imm_lsl,
)

# Four architectural instructions (Mops) can be decoded/renamed per cycle, and
# those Mops may expand to at most eight uops at dispatch.
issue_rate = 4
dispatch_rate = 8
llvm_mca_target = "cortex-a76"

# Section 4.9 recommends aligning loop branch targets to a 32-byte instruction
# memory region. Slothy.optimize_loop() consumes this target hook.
loop_alignment_bytes = 32


class ExecutionUnit(Enum):
    """Execution pipelines in the Cortex-A76 model."""

    SCALAR_S0 = 0
    SCALAR_S1 = 1
    SCALAR_M = 2
    LSU0 = 3
    LSU1 = 4
    VEC0 = 5
    VEC1 = 6

    def __repr__(self):
        return self.name

    @classmethod
    def I(cls):  # noqa: E743
        return [cls.SCALAR_S0, cls.SCALAR_S1, cls.SCALAR_M]

    @classmethod
    def M(cls):
        return [cls.SCALAR_M]

    @classmethod
    def L(cls):
        return [cls.LSU0, cls.LSU1]

    @classmethod
    def V(cls):
        return [cls.VEC0, cls.VEC1]

    @classmethod
    def V0(cls):
        return [cls.VEC0]

    @classmethod
    def V1(cls):
        return [cls.VEC1]


def _one_of(units, duration=1):
    return [{unit: duration} for unit in units]


def _one_from_each(*groups):
    """Build alternatives selecting one ``(unit, duration)`` from each group."""
    alternatives = []
    for choices in product(*groups):
        alternatives.append({unit: duration for unit, duration in choices})
    return alternatives


def _fixed(**usages):
    return [
        {
            getattr(ExecutionUnit, unit_name): duration
            for unit_name, duration in usages.items()
        }
    ]


def _with_integer_writeback(inst, alternatives):
    has_writeback = inst.increment is not None or isinstance(inst, Q_Ld2_Lane_Post_Inc)
    if not has_writeback:
        return alternatives

    expanded = []
    for usage in alternatives:
        for unit in ExecutionUnit.I():
            choice = usage.copy()
            choice[unit] = max(choice.get(unit, 0), 1)
            expanded.append(choice)
    return expanded


def _datatype(inst):
    datatype = inst.datatype
    if isinstance(datatype, list):
        datatype = datatype[0]
    return datatype.lower() if isinstance(datatype, str) else None


def _is_q_form(inst):
    try:
        return inst.is_q_form_vector_instruction()
    except Exception:
        return False


def _is_q_d_form(inst):
    return _is_q_form(inst) and _datatype(inst) == "2d"


# The regular single-uop instruction map. Memory instructions with heterogeneous
# uops and form-dependent costs are handled in get_resource_usages().
execution_units = {
    (Ldp_X, Ldp_W, Ldr_X, Str_X, Stp_X, Stp_W, Ldr_D, Ldr_Q): ExecutionUnit.L(),
    ASimdCompare: ExecutionUnit.V(),
    (vtbl, vtbl_2): ExecutionUnit.V(),
    (Vzip, Vrev, uaddlp): ExecutionUnit.V(),
    AArch64NeonCount: ExecutionUnit.V(),
    vmov: ExecutionUnit.V(),
    VecToGprMov: ExecutionUnit.V1(),
    Transpose: ExecutionUnit.V(),
    vmovi: ExecutionUnit.V(),
    (vadd, vsub): ExecutionUnit.V(),
    vxtn: ExecutionUnit.V(),
    VShiftImmediateBasic: ExecutionUnit.V1(),
    VShiftRegBasic: ExecutionUnit.V1(),
    (AArch64NeonShiftInsert, vsrshr, vusra): ExecutionUnit.V1(),
    AESInstruction: ExecutionUnit.V0(),
    (Vmul, Vmla, Vqdmulh, Vmull, Vmlal): ExecutionUnit.V0(),
    AArch64NeonLogical: ExecutionUnit.V(),
    vext: ExecutionUnit.V(),
    (
        AArch64BasicArithmetic,
        AArch64ConditionalSelect,
        AArch64ConditionalCompare,
        AArch64Logical,
        AArch64LogicalShifted,
        AArch64Move,
        AArch64Shift,
        Tst,
    ): ExecutionUnit.I(),
    AArch64ShiftedArithmetic: ExecutionUnit.M(),
    (fmov_0, fmov_0_force_output, fmov_1, fmov_1_force_output): ExecutionUnit.M(),
    (fmov_s_form, fmov_d_form): ExecutionUnit.V1(),
    umull_wform: ExecutionUnit.M(),
    (AArch64HighMultiply, AArch64Multiply, AArch64CRC32): ExecutionUnit.M(),
    (vdup, vdup_w): ExecutionUnit.M(),
    # The 8B/8H reduction uses both vector pipelines.
    vuaddlv_sform: [[ExecutionUnit.VEC0, ExecutionUnit.VEC1]],
    mov_wtov_s: ExecutionUnit.V(),
    mov_vtov_d: ExecutionUnit.V(),
    lsr: ExecutionUnit.I(),
    movk_imm_lsl: ExecutionUnit.I(),
}


inverse_throughput = {
    (Ldr_X, Str_X, Ldr_D, Ldr_Q): 1,
    (Ldp_X, Stp_X): 2,
    (Ldp_W, Stp_W): 1,
    AArch64NeonCount: 1,
    (Vzip, uaddlp, Vrev): 1,
    VecToGprMov: 1,
    (vadd, vsub, vmov, ASimdCompare, Transpose): 1,
    AESInstruction: 1,
    AArch64NeonLogical: 1,
    vext: 1,
    (vmovi, vxtn): 1,
    (VShiftImmediateBasic, VShiftRegBasic): 1,
    (AArch64NeonShiftInsert, vsrshr): 1,
    (Vmul, Vmla, Vqdmulh): 2,
    vusra: 1,
    (vtbl, vtbl_2, Vmull, Vmlal): 1,
    (
        AArch64BasicArithmetic,
        AArch64ConditionalSelect,
        AArch64ConditionalCompare,
        AArch64Logical,
        AArch64LogicalShifted,
        AArch64Move,
        AArch64Shift,
        Tst,
        AArch64ShiftedArithmetic,
    ): 1,
    (fmov_0, fmov_0_force_output, fmov_1, fmov_1_force_output): 1,
    (fmov_s_form, fmov_d_form): 1,
    AArch64HighMultiply: 4,
    AArch64Multiply: 3,
    AArch64CRC32: 1,
    (vdup, vdup_w, umull_wform, vuaddlv_sform): 1,
    (mov_wtov_s, mov_vtov_d, lsr, movk_imm_lsl): 1,
}


default_latencies = {
    (Ldp_X, Ldp_W, Ldr_X): 4,
    (Str_X, Stp_X): 1,
    Stp_W: 1,
    (Vzip, Vrev, uaddlp): 2,
    VecToGprMov: 2,
    ASimdCompare: 2,
    vxtn: 2,
    AArch64NeonCount: 2,
    AESInstruction: 2,
    AArch64NeonLogical: 2,
    vext: 2,
    Transpose: 2,
    (vadd, vsub, vmov, vmovi): 2,
    (Vmul, Vmla, Vqdmulh): 5,
    vusra: 4,
    (Vmull, Vmlal): 4,
    (VShiftImmediateBasic, VShiftRegBasic, AArch64NeonShiftInsert): 2,
    vsrshr: 4,
    (
        AArch64BasicArithmetic,
        AArch64ConditionalSelect,
        AArch64ConditionalCompare,
        AArch64Logical,
        AArch64LogicalShifted,
        AArch64Move,
        AArch64Shift,
        Tst,
    ): 1,
    AArch64ShiftedArithmetic: 2,
    (fmov_0, fmov_0_force_output, fmov_1, fmov_1_force_output): 3,
    (fmov_s_form, fmov_d_form): 2,
    AArch64HighMultiply: 5,
    AArch64Multiply: 4,
    AArch64CRC32: 2,
    (vdup, vdup_w): 3,
    umull_wform: 2,
    (vtbl, vtbl_2): 2,
    vuaddlv_sform: 5,
    mov_wtov_s: 5,
    mov_vtov_d: 2,
    (lsr, movk_imm_lsl): 1,
}


def _structured_load_usage(inst):
    """Return L/V pipeline-cycle pressure for LD2/LD3/LD4."""
    q_form = _is_q_form(inst)
    if isinstance(inst, Ld2):
        if q_form:
            return _fixed(LSU0=1, LSU1=1, VEC0=1, VEC1=1)
        return _one_from_each(
            [(unit, 2) for unit in ExecutionUnit.L()],
            [(unit, 2) for unit in ExecutionUnit.V()],
        )
    if isinstance(inst, Ld3):
        return _fixed(LSU0=2, LSU1=2, VEC0=2, VEC1=2)
    if q_form:
        # Four 128-bit load/deinterleave operations: the V side is the
        # documented 1/5-throughput bottleneck.
        return _fixed(LSU0=2, LSU1=2, VEC0=5, VEC1=5)
    # D-form LD4 has throughput 2/7. Selecting one of the two V pipes for
    # seven cycles represents that fractional throughput exactly.
    alternatives = []
    for vector_unit in ExecutionUnit.V():
        alternatives.append(
            {
                ExecutionUnit.LSU0: 1,
                ExecutionUnit.LSU1: 1,
                vector_unit: 7,
            }
        )
    return alternatives


def _structured_store_usage(inst):
    """Return L/V pipeline-cycle pressure for ST2/ST3/ST4."""
    q_form = _is_q_form(inst)
    if isinstance(inst, St2):
        duration = 2 if q_form else 1
    elif isinstance(inst, St3):
        duration = 3 if q_form else 2
    elif q_form:
        # Q-form ST4 is 1/6 for B/H/S elements and 1/4 for D elements.
        duration = 4 if _is_q_d_form(inst) else 6
    else:
        duration = 3
    return _fixed(
        LSU0=max(1, (duration + 1) // 3),
        LSU1=max(1, (duration + 1) // 3),
        VEC0=duration,
        VEC1=duration,
    )


def get_resource_usages(inst):
    """Return alternative per-pipeline occupancy mappings for ``inst``."""
    if isinstance(inst, (Ld2, Ld3, Ld4)):
        usages = _structured_load_usage(inst)
    elif isinstance(inst, (St2, St3, St4)):
        usages = _structured_store_usage(inst)
    elif isinstance(inst, q_st1_4_with_postinc):
        usages = _fixed(LSU0=2, LSU1=2, VEC0=4, VEC1=4)
    elif isinstance(inst, (Q_Ld2_Lane_Post_Inc, q_ld2_lane_s)):
        usages = _one_from_each(
            [(unit, 2) for unit in ExecutionUnit.L()],
            [(unit, 2) for unit in ExecutionUnit.V()],
        )
    elif isinstance(inst, q_ldr1_stack):
        usages = _one_from_each(
            [(unit, 1) for unit in ExecutionUnit.L()],
            [(unit, 1) for unit in ExecutionUnit.V()],
        )
    elif isinstance(inst, q_ld1_2):
        usages = _fixed(LSU0=1, LSU1=1)
    elif isinstance(inst, Ldp_Q):
        usages = _fixed(LSU0=1, LSU1=1)
    elif isinstance(inst, Stp_Q):
        vector_duration = 2 if inst.increment is not None else 4
        usages = _one_from_each(
            [(unit, 2) for unit in ExecutionUnit.L()],
            [(unit, vector_duration) for unit in ExecutionUnit.V()],
        )
    elif isinstance(inst, Str_Q):
        # One address uop and one store-data uop. Two cycles of aggregate V
        # occupancy reproduce the documented one-Q-store-per-cycle limit.
        usages = _one_from_each(
            [(unit, 1) for unit in ExecutionUnit.L()],
            [(unit, 2) for unit in ExecutionUnit.V()],
        )
    else:
        instclass = find_class(inst)
        units = lookup_multidict(execution_units, inst, instclass)
        duration = get_inverse_throughput(inst)
        if len(units) == 1 and isinstance(units[0], list):
            usages = [{unit: duration for unit in units[0]}]
        else:
            usages = _one_of(units, duration)

    if isinstance(
        inst,
        (
            Ldp_X,
            Ldp_W,
            Ldr_X,
            Str_X,
            Stp_X,
            Stp_W,
            Ldr_D,
            Ldr_Q,
            Ldp_Q,
            Str_Q,
            Stp_Q,
            Ld2,
            Ld3,
            Ld4,
            St2,
            St3,
            St4,
            q_st1_4_with_postinc,
            Q_Ld2_Lane_Post_Inc,
        ),
    ) and not isinstance(inst, q_ld1_2):
        usages = _with_integer_writeback(inst, usages)
    return usages


def get_dispatch_uops(inst):
    """Return dispatch uops, including the optional address-update uop."""
    split_classes = (
        Str_X,
        Stp_X,
        Stp_W,
        Str_Q,
        Stp_Q,
        Ldp_X,
        Ldp_W,
        Ldp_Q,
        Ld2,
        Ld3,
        Ld4,
        St2,
        St3,
        St4,
        q_st1_4_with_postinc,
        Q_Ld2_Lane_Post_Inc,
        q_ld2_lane_s,
        q_ldr1_stack,
    )
    uops = 2 if isinstance(inst, split_classes) else 1
    if inst.increment is not None or isinstance(inst, Q_Ld2_Lane_Post_Inc):
        uops += 1
    return uops


def get_latency(src, out_idx, dst):
    """Return producer-to-consumer latency, including A76 forwarding paths."""
    _ = out_idx
    instclass_src = find_class(src)
    instclass_dst = find_class(dst)

    # FP/ASIMD loads need one more cycle than a standard integer load to
    # forward into the V pipelines (SWOG sections 3.13 and 3.18).
    if isinstance(src, (Ldr_D, Ldr_Q)):
        return 5
    if isinstance(src, Ldp_Q):
        return 7
    if isinstance(src, q_ld1_2):
        return 5
    if isinstance(src, (Q_Ld2_Lane_Post_Inc, q_ld2_lane_s, q_ldr1_stack)):
        return 7
    if isinstance(src, Ld2):
        return 7
    if isinstance(src, Ld3):
        return 8
    if isinstance(src, Ld4):
        return 10 if _is_q_form(src) else 8
    if isinstance(src, Str_Q):
        return 2
    if isinstance(src, Stp_Q):
        return 3
    if isinstance(src, St2):
        return 5 if _is_q_form(src) else 4
    if isinstance(src, St3):
        return 6 if _is_q_form(src) else 5
    if isinstance(src, St4):
        if not _is_q_form(src):
            return 7
        return 6 if _is_q_d_form(src) else 9
    if isinstance(src, q_st1_4_with_postinc):
        return 5

    latency = lookup_multidict(default_latencies, src, instclass_src)

    # D-form integer multiplies/MLAs are one cycle shorter and twice as
    # throughput-friendly as Q-form operations on Cortex-A76.
    if isinstance(src, (Vmul, Vmla, Vqdmulh)) and not _is_q_form(src):
        latency = 4

    # Late forwarding into an accumulate operand.
    if (
        instclass_src in [vmul, vmul_lane]
        and instclass_dst in [vmla, vmla_lane, vmls, vmls_lane]
        and src.args_out[0] == dst.args_in_out[0]
    ):
        return 2 if _is_q_form(src) else 1
    if (
        instclass_src in [vmla, vmla_lane, vmls, vmls_lane]
        and instclass_dst in [vmla, vmla_lane, vmls, vmls_lane]
        and src.args_in_out[0] == dst.args_in_out[0]
    ):
        return 2 if _is_q_form(src) else 1
    # Pi 5 Cortex-A76 r4p1 measurements establish a 3-cycle signed
    # SMULL->SMLAL edge when the SMULL destination is the SMLAL accumulator.
    if (
        instclass_src in [vsmull, vsmull2, vsmull_lane, vsmull2_lane]
        and instclass_dst in [vsmlal, vsmlal2, vsmlal_lane, vsmlal2_lane]
        and src.args_out[0] == dst.args_in_out[0]
    ):
        return 3
    # Retain the related N1 approximation for unmeasured unsigned forms.
    if (
        instclass_src in all_subclass_leaves(Vmull)
        and instclass_dst in all_subclass_leaves(Vmlal)
        and src.args_out[0] == dst.args_in_out[0]
    ):
        return 1
    if (
        instclass_src in all_subclass_leaves(Vmlal)
        and instclass_dst in all_subclass_leaves(Vmlal)
        and src.args_in_out[0] == dst.args_in_out[0]
    ):
        return 1
    if (
        isinstance(src, vusra)
        and isinstance(dst, vusra)
        and src.args_in_out[0] == dst.args_in_out[0]
    ):
        return 1
    if (
        instclass_src in all_subclass_leaves(AArch64CRC32)
        and instclass_dst in all_subclass_leaves(AArch64CRC32)
        and src.args_out[0] == dst.args_in[0]
    ):
        return 1
    return latency


def get_inverse_throughput(inst):
    """Return legacy integer inverse throughput.

    The rich resource interface is authoritative. This compatibility function
    returns the reciprocal throughput for integral cases and the conservative
    ceiling (4) for D-form LD4's documented 2/7 throughput.
    """
    if isinstance(inst, Ld4):
        return 5 if _is_q_form(inst) else 4
    if isinstance(inst, Ld3):
        return 2
    if isinstance(inst, Ld2):
        return 1
    if isinstance(inst, St4):
        if not _is_q_form(inst):
            return 3
        return 4 if _is_q_d_form(inst) else 6
    if isinstance(inst, St3):
        return 3 if _is_q_form(inst) else 2
    if isinstance(inst, St2):
        return 2 if _is_q_form(inst) else 1
    if isinstance(inst, q_st1_4_with_postinc):
        return 4
    if isinstance(inst, (Q_Ld2_Lane_Post_Inc, q_ld2_lane_s)):
        return 1
    if isinstance(inst, q_ldr1_stack):
        return 1
    if isinstance(inst, q_ld1_2):
        return 1
    if isinstance(inst, Ldp_Q):
        return 1
    if isinstance(inst, Str_Q):
        return 1
    if isinstance(inst, Stp_Q):
        return 1 if inst.increment is not None else 2
    if isinstance(inst, (Vmul, Vmla, Vqdmulh)):
        return 2 if _is_q_form(inst) else 1
    instclass = find_class(inst)
    return lookup_multidict(inverse_throughput, inst, instclass)


def get_units(inst):
    """Return the legacy unit alternatives derived from resource usages."""
    usages = get_resource_usages(inst)
    alternatives = [list(usage) for usage in usages]
    if len(alternatives) == 1:
        return [alternatives[0]]
    return alternatives


def _partial_single_word_outputs(inst):
    """Return output/in-out indices written as 32-bit pieces of a V register."""
    output_indices = []
    inout_indices = []
    pattern = getattr(inst, "pattern", "").lower()

    for idx, (name, ty) in enumerate(getattr(inst, "pattern_outputs", [])):
        if ty != RegisterType.NEON:
            continue
        if name.lower().startswith("s") or f"<{name.lower()}>.s" in pattern:
            output_indices.append(idx)
    for idx, (name, ty) in enumerate(getattr(inst, "pattern_in_outs", [])):
        if ty != RegisterType.NEON:
            continue
        if name.lower().startswith("s") or f"<{name.lower()}>.s" in pattern:
            inout_indices.append(idx)

    # Datatype placeholders in lane forms are resolved only after parsing.
    if getattr(inst, "index", None) is not None and _datatype(inst) == "s":
        output_indices.extend(
            idx for idx, ty in enumerate(inst.arg_types_out) if ty == RegisterType.NEON
        )
        inout_indices.extend(
            idx
            for idx, ty in enumerate(inst.arg_types_in_out)
            if ty == RegisterType.NEON
        )
    return sorted(set(output_indices)), sorted(set(inout_indices))


def is_partial_register_write(inst):
    outputs, inouts = _partial_single_word_outputs(inst)
    return bool(outputs or inouts)


def is_partial_register_hazard_consumer(inst):
    """Whether ``inst`` has a V uop with more than one Q-register source."""
    # The guide explicitly excludes store-data and MOV uops from this
    # forwarding hazard.
    if isinstance(
        inst,
        (
            Str_Q,
            Stp_Q,
            St2,
            St3,
            St4,
            q_st1_4_with_postinc,
            vmov,
            mov_vtov_d,
        ),
    ):
        return False
    q_form = _is_q_form(inst) or isinstance(inst, (Stp_Q, q_st1_4_with_postinc))
    if not q_form:
        return False
    vector_sources = sum(
        ty == RegisterType.NEON for ty in inst.arg_types_in + inst.arg_types_in_out
    )
    if vector_sources <= 1:
        return False
    return any(
        unit in ExecutionUnit.V()
        for usage in get_resource_usages(inst)
        for unit in usage
    )


def _enforce(constraint, condition):
    if condition is not None:
        constraint.OnlyEnforceIf(condition)


def _add_dispatch_repair_for_consumers(slothy, consumers):
    consumers = list(
        dict.fromkeys(
            consumer
            for consumer in consumers
            if is_partial_register_hazard_consumer(consumer.inst)
        )
    )
    if not consumers:
        return

    first_vars = [None]
    if len(consumers) > 1:
        first_vars = [
            slothy._NewBoolVar(f"{consumer.varname()}_first_partial_q_consumer")
            for consumer in consumers
        ]
        slothy._AddExactlyOne(first_vars)
        for consumer, first_var in zip(consumers, first_vars):
            for other in consumers:
                if other is consumer:
                    continue
                slothy._Add(
                    consumer.program_start_var < other.program_start_var
                ).OnlyEnforceIf(first_var)

    for consumer, first_var in zip(consumers, first_vars):
        # The repair occupies three complete dispatch cycles before the first
        # qualifying consumer. Instructions may still share the successful
        # dispatch cycle with that consumer.
        _enforce(slothy._Add(consumer.cycle_start_var >= 3), first_var)
        for other in slothy._get_nodes():
            if other is consumer:
                continue
            before_window = slothy._NewBoolVar(
                f"{other.varname()}_outside_partial_stall_{consumer.varname()}"
            )
            before = slothy._Add(other.cycle_start_var <= consumer.cycle_start_var - 4)
            after = slothy._Add(other.cycle_start_var >= consumer.cycle_start_var)
            if first_var is None:
                before.OnlyEnforceIf(before_window)
                after.OnlyEnforceIf(before_window.Not())
            else:
                before.OnlyEnforceIf([first_var, before_window])
                after.OnlyEnforceIf([first_var, before_window.Not()])


def _add_partial_register_dispatch_stalls(slothy):
    for node in slothy._get_nodes():
        output_indices, inout_indices = _partial_single_word_outputs(node.inst)
        for idx in output_indices:
            _add_dispatch_repair_for_consumers(slothy, node.dst_out[idx])
        for idx in inout_indices:
            _add_dispatch_repair_for_consumers(slothy, node.dst_in_out[idx])


def add_further_constraints(slothy):
    if slothy.config.constraints.functional_only:
        return
    _add_partial_register_dispatch_stalls(slothy)


def has_min_max_objective(config):
    _ = config
    return False


def get_min_max_objective(slothy):
    _ = slothy
