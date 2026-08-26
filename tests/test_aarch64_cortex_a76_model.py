# SPDX-License-Identifier: MIT
"""Focused tests for measured Cortex-A76 SMULL/SMLAL behavior."""

import unittest

from slothy.helper import SourceLine
from slothy.targets import query
from slothy.targets.aarch64 import aarch64_neon
from slothy.targets.aarch64 import cortex_a76


def parse_instruction(text):
    """Parse one AArch64 instruction."""
    return aarch64_neon.Instruction.parser(SourceLine(text))[0]


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

    def test_target_is_registered(self):
        self.assertIs(
            query.Archery.get_target("Arm_Cortex_A76"),
            cortex_a76,
        )

    def test_signed_widening_multiply_timing(self):
        for smull_text, dependent_smlal_text, independent_smlal_text in self.FORMS:
            with self.subTest(smull=smull_text):
                smull = parse_instruction(smull_text)
                dependent_smlal = parse_instruction(dependent_smlal_text)
                independent_smlal = parse_instruction(independent_smlal_text)

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
        smull = parse_instruction(self.FORMS[0][0])
        smlal = parse_instruction(self.FORMS[0][1])
        expected = [cortex_a76.ExecutionUnit.VEC0]

        self.assertEqual(cortex_a76.get_units(smull), expected)
        self.assertEqual(cortex_a76.get_units(smlal), expected)


if __name__ == "__main__":
    unittest.main()
