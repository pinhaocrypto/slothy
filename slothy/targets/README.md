# SLOTHY Architecture and Microarchitecture models

This directory contains experimental architecture and microarchitecture models.

Currently, the following architectures have experimental support:

* Armv8.1-M (mostly focused on Helium SIMD instructions)
* AArch64 (scalar and Neon, yet incomplete)

The following microarchitectures have experimental support:

* Cortex-M55 (Armv8.1-M+Helium)
* Cortex-M85 (Armv8.1-M+Helium)
* Cortex-A55 (AArch64): Largely complete and accurate
* Cortex-A72 (AArch64)
* Cortex-A76 (AArch64): Experimental; signed SMULL/SMLAL timing validated on
  Raspberry Pi 5 Cortex-A76 r4p1
