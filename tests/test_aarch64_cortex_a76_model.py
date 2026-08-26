#
# Copyright (c) 2026 Slothy contributors
# SPDX-License-Identifier: MIT
#

import logging
import unittest

from slothy import Slothy
from slothy.core.config import Config
from slothy.core.core import SlothyBase
from slothy.helper import SourceLine
from slothy.targets import query
from slothy.targets.aarch64 import aarch64_neon as aarch64
from slothy.targets.aarch64 import cortex_a76


def parse(text):
    return aarch64.Instruction.parser(SourceLine(text))[0]


class CortexA76ModelTest(unittest.TestCase):
    def test_target_is_registered(self):
        self.assertIs(
            query.Archery.get_target("Arm_Cortex_A76"),
            cortex_a76,
        )

    def test_asimd_store_uses_l_and_v_resources(self):
        instruction = parse("str q0, [x0]")
        usages = cortex_a76.get_resource_usages(instruction)
        self.assertEqual(len(usages), 4)
        for usage in usages:
            self.assertTrue(set(usage) & set(cortex_a76.ExecutionUnit.L()))
            self.assertTrue(set(usage) & set(cortex_a76.ExecutionUnit.V()))

    def test_ld4_st4_q_form_costs(self):
        load = parse("ld4 {v0.16b, v1.16b, v2.16b, v3.16b}, [x0]")
        store = parse("st4 {v4.16b, v5.16b, v6.16b, v7.16b}, [x1]")

        self.assertEqual(cortex_a76.get_dispatch_uops(load), 2)
        self.assertEqual(cortex_a76.get_inverse_throughput(load), 5)
        self.assertEqual(
            cortex_a76.get_resource_usages(load),
            [
                {
                    cortex_a76.ExecutionUnit.LSU0: 2,
                    cortex_a76.ExecutionUnit.LSU1: 2,
                    cortex_a76.ExecutionUnit.VEC0: 5,
                    cortex_a76.ExecutionUnit.VEC1: 5,
                }
            ],
        )
        self.assertEqual(cortex_a76.get_dispatch_uops(store), 2)
        self.assertEqual(cortex_a76.get_inverse_throughput(store), 6)
        self.assertEqual(
            cortex_a76.get_resource_usages(store),
            [
                {
                    cortex_a76.ExecutionUnit.LSU0: 2,
                    cortex_a76.ExecutionUnit.LSU1: 2,
                    cortex_a76.ExecutionUnit.VEC0: 6,
                    cortex_a76.ExecutionUnit.VEC1: 6,
                }
            ],
        )

    def test_ld4_st4_costs_cover_d_and_q_d_forms(self):
        cases = (
            (
                "ld4 {v0.8b, v1.8b, v2.8b, v3.8b}, [x0]",
                8,
                4,
                2,
            ),
            (
                "ld4 {v0.2d, v1.2d, v2.2d, v3.2d}, [x0]",
                10,
                5,
                2,
            ),
            (
                "st4 {v0.8b, v1.8b, v2.8b, v3.8b}, [x0]",
                7,
                3,
                2,
            ),
            (
                "st4 {v0.2d, v1.2d, v2.2d, v3.2d}, [x0]",
                6,
                4,
                2,
            ),
        )
        for text, latency, inverse_throughput, uops in cases:
            with self.subTest(instruction=text):
                instruction = parse(text)
                self.assertEqual(
                    cortex_a76.get_latency(instruction, 0, instruction),
                    latency,
                )
                self.assertEqual(
                    cortex_a76.get_inverse_throughput(instruction),
                    inverse_throughput,
                )
                self.assertEqual(
                    cortex_a76.get_dispatch_uops(instruction),
                    uops,
                )
                usage = cortex_a76.get_resource_usages(instruction)
                self.assertTrue(
                    all(
                        set(alternative) & set(cortex_a76.ExecutionUnit.L())
                        for alternative in usage
                    )
                )
                self.assertTrue(
                    all(
                        set(alternative) & set(cortex_a76.ExecutionUnit.V())
                        for alternative in usage
                    )
                )

    def test_ld4_writeback_adds_integer_uop(self):
        instruction = parse("ld4 {v0.16b, v1.16b, v2.16b, v3.16b}, [x0], #64")
        usages = cortex_a76.get_resource_usages(instruction)
        self.assertEqual(cortex_a76.get_dispatch_uops(instruction), 3)
        self.assertEqual(len(usages), 3)
        for usage in usages:
            self.assertTrue(set(usage) & set(cortex_a76.ExecutionUnit.I()))
            self.assertTrue(set(usage) & set(cortex_a76.ExecutionUnit.L()))
            self.assertTrue(set(usage) & set(cortex_a76.ExecutionUnit.V()))

    def test_vector_load_forwarding(self):
        vector_load = parse("ldr q0, [x1]")
        vector_consumer = parse("add v2.4s, v0.4s, v3.4s")
        integer_load = parse("ldr x0, [x1]")
        integer_consumer = parse("add x2, x0, x3")

        self.assertEqual(
            cortex_a76.get_latency(vector_load, 0, vector_consumer),
            5,
        )
        self.assertEqual(
            cortex_a76.get_latency(integer_load, 0, integer_consumer),
            4,
        )

    def test_partial_register_hazard_classification(self):
        partial_write = parse("uaddlv s0, v1.8h")
        multi_q_consumer = parse("add v2.4s, v0.4s, v3.4s")
        single_q_consumer = parse("cnt v2.16b, v0.16b")
        store_consumer = parse("st4 {v0.4s, v1.4s, v2.4s, v3.4s}, [x0]")

        self.assertTrue(cortex_a76.is_partial_register_write(partial_write))
        self.assertTrue(
            cortex_a76.is_partial_register_hazard_consumer(multi_q_consumer)
        )
        self.assertFalse(
            cortex_a76.is_partial_register_hazard_consumer(single_q_consumer)
        )
        self.assertFalse(cortex_a76.is_partial_register_hazard_consumer(store_consumer))

    def test_partial_register_repair_reserves_three_dispatch_cycles(self):
        config = Config(aarch64, cortex_a76, logging.getLogger("a76-repair-test"))
        config.outputs = {"v2", "x4", "x7", "x10", "x13"}
        config.inputs_are_outputs = False
        config.selftest = False
        config.variable_size = True
        config.constraints.stalls_allowed = 12
        config.constraints.allow_renaming = False

        source = SourceLine.read_multiline("""
            uaddlv s0, v1.8h
            add x4, x5, x6
            add x7, x8, x9
            add x10, x11, x12
            add x13, x14, x15
            add v2.4s, v0.4s, v3.4s
            """)
        core = SlothyBase(
            aarch64,
            cortex_a76,
            logger=logging.getLogger("a76-repair-test"),
            config=config,
        )
        self.assertTrue(core.optimize(source))

        consumer = next(
            node
            for node in core._model.tree.nodes
            if isinstance(node.inst, aarch64.vadd)
        )
        occupied_cycles = {
            node.real_pos_cycle
            for node in core._model.tree.nodes
            if node is not consumer
        }
        self.assertTrue(
            set(
                range(
                    consumer.real_pos_cycle - 3,
                    consumer.real_pos_cycle,
                )
            ).isdisjoint(occupied_cycles)
        )

    def test_optimize_loop_inserts_32_byte_alignment(self):
        slothy = Slothy(
            aarch64,
            cortex_a76,
            logger=logging.getLogger("a76-alignment-test"),
        )
        slothy.load_source_raw("""
            loop:
                add x0, x0, x1
                subs x2, x2, #1
                b.ne loop
            """)
        slothy.config.selftest = False
        slothy.config.inputs_are_outputs = True
        slothy.config.constraints.functional_only = True
        slothy.config.constraints.allow_renaming = False
        slothy.optimize_loop("loop")

        self.assertIn(
            ".p2align 5\nloop:",
            slothy.get_source_as_string(),
        )


class CortexA76SmullForwardingTest(unittest.TestCase):
    """Check the Pi 5 Cortex-A76 r4p1 measurements encoded by the model."""

    FORMS = (
        (
            "smull v0.4s, v1.4h, v2.4h",
            "smlal v0.4s, v3.4h, v4.4h",
            "smlal v5.4s, v3.4h, v4.4h",
        ),
        (
            "smull v0.2d, v1.2s, v2.2s",
            "smlal v0.2d, v3.2s, v4.2s",
            "smlal v5.2d, v3.2s, v4.2s",
        ),
    )

    def test_signed_widening_multiply_timing(self):
        for smull_text, dependent_smlal_text, independent_smlal_text in self.FORMS:
            with self.subTest(smull=smull_text):
                smull = parse(smull_text)
                dependent_smlal = parse(dependent_smlal_text)
                independent_smlal = parse(independent_smlal_text)

                self.assertEqual(cortex_a76.get_inverse_throughput(smull), 1)
                self.assertEqual(cortex_a76.get_inverse_throughput(dependent_smlal), 1)
                self.assertEqual(
                    cortex_a76.get_latency(smull, 0, independent_smlal),
                    4,
                )
                self.assertEqual(
                    cortex_a76.get_latency(smull, 0, dependent_smlal),
                    3,
                )
                self.assertEqual(
                    cortex_a76.get_latency(dependent_smlal, 0, dependent_smlal),
                    1,
                )

    def test_signed_widening_multiply_uses_v0(self):
        smull = parse(self.FORMS[0][0])
        smlal = parse(self.FORMS[0][1])
        expected = [[cortex_a76.ExecutionUnit.VEC0]]

        self.assertEqual(cortex_a76.get_units(smull), expected)
        self.assertEqual(cortex_a76.get_units(smlal), expected)


if __name__ == "__main__":
    unittest.main()
