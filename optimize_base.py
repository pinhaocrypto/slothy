import argparse
import re
import time
from pathlib import Path

from slothy import Slothy
from slothy.targets.aarch64 import aarch64_neon as AArch64_Neon
from slothy.targets.aarch64 import cortex_a55 as Target_CortexA55
from slothy.targets.aarch64 import cortex_a72_frontend as Target_CortexA72
from slothy.targets.aarch64 import (
    apple_m1_firestorm_experimental as Target_AppleM1Firestorm,
)
from slothy.targets.aarch64 import (
    apple_m1_icestorm_experimental as Target_AppleM1Icestorm,
)

ROOT = Path(__file__).resolve().parent

TARGETS = {
    "a55": Target_CortexA55,
    "a72": Target_CortexA72,
    "m1-firestorm": Target_AppleM1Firestorm,
    "m1-icestorm": Target_AppleM1Icestorm,
}

DEFAULT_LOOPS = {
    "base.s": ["_looptop", "_looptop_add", "_looptop_baseinv_1"],
    "base_gt.s": ["Lgt_basemul_loop", "Lgt_basemul_add_loop"],
}

VECTOR_RANGE_RE = re.compile(
    r"\{\s*v(?P<start>[0-9]+)\.(?P<dt>[0-9]+[bhsdBHSD])\s*-\s*"
    r"v(?P<end>[0-9]+)\.(?P=dt)\s*\}"
)


def expand_vector_ranges(line: str) -> str:
    def repl(match: re.Match[str]) -> str:
        start = int(match.group("start"))
        end = int(match.group("end"))
        if end < start:
            raise ValueError(f"Invalid vector register range in line: {line}")
        dt = match.group("dt")
        return "{" + ", ".join(f"v{i}.{dt}" for i in range(start, end + 1)) + "}"

    return VECTOR_RANGE_RE.sub(repl, line)


def normalize_source(source: str) -> str:
    normalized = []
    for line in source.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("#"):
            indent = line[: len(line) - len(stripped)]
            line = f"{indent}//{stripped[1:]}"
        normalized.append(expand_vector_ranges(line))
    return "\n".join(normalized) + "\n"


def output_path_for(source_path: Path, suffix: str) -> Path:
    return source_path.with_name(f"{source_path.stem}.{suffix}{source_path.suffix}")


def new_slothy(target, timeout: int, reserved_regs: list[str]) -> Slothy:
    slothy = Slothy(AArch64_Neon, target)
    slothy.config.with_llvm_mca = False
    slothy.config.selftest = False
    slothy.config.timeout = timeout
    slothy.config.variable_size = True
    slothy.config.inputs_are_outputs = True
    slothy.config.reserved_regs = reserved_regs
    slothy.config.constraints.stalls_first_attempt = 128
    return slothy


def optimize_file(
    source_path: Path, loops: list[str], args: argparse.Namespace
) -> Path:
    target = TARGETS[args.target]
    reserved_regs = [f"x{i}" for i in range(9)] + ["x30", "sp"]
    slothy = new_slothy(target, args.timeout, reserved_regs)
    slothy.load_source_raw(normalize_source(source_path.read_text(encoding="utf8")))

    print(f"{source_path.name}: target={args.target}, timeout={args.timeout}s")
    start = time.perf_counter()
    for loop in loops:
        loop_start = time.perf_counter()
        print(f"  optimize_loop({loop})")
        slothy.optimize_loop(loop)
        print(f"    done in {time.perf_counter() - loop_start:.3f}s")

    output_path = output_path_for(source_path, args.suffix)
    slothy.write_source_to_file(str(output_path))
    print(f"  wrote {output_path.name} in {time.perf_counter() - start:.3f}s")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Optimize base.s and base_gt.s with Slothy."
    )
    parser.add_argument(
        "sources",
        nargs="*",
        type=Path,
        default=[ROOT / "base.s", ROOT / "base_gt.s"],
        help="Assembly source files to optimize.",
    )
    parser.add_argument(
        "--target",
        choices=TARGETS,
        default="a72",
        help="Microarchitecture model to optimize for.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="Per-solver-call timeout in seconds.",
    )
    parser.add_argument(
        "--suffix",
        default="opt",
        help="Output filename suffix, e.g. base.opt.s.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for source_path in args.sources:
        source_path = source_path.resolve()
        loops = DEFAULT_LOOPS.get(source_path.name)
        if loops is None:
            known = ", ".join(sorted(DEFAULT_LOOPS))
            raise ValueError(
                f"No loop list for {source_path.name}; known files: {known}"
            )
        optimize_file(source_path, loops, args)


if __name__ == "__main__":
    main()
