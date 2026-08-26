import re
import time
from pathlib import Path

from slothy import Slothy
from slothy.targets.aarch64 import aarch64_neon as AArch64_Neon
from slothy.targets.aarch64 import cortex_a72_frontend as Target_CortexA72

ROOT = Path(__file__).resolve().parent
SOURCE_PATH = ROOT / "my_32ntt.s"
ALLOC_PATH = ROOT / "my_32ntt.alloc.s"
OUTPUT_PATH = ROOT / "my_32ntt.opt.s"

REGION_START_RE = re.compile(r"^\s*([A-Za-z0-9_.$]+)_slothy_start:\s*$")
REGION_END_RE = re.compile(r"^\s*([A-Za-z0-9_.$]+)_slothy_end:\s*$")

RESERVED_GPRS = [
    "x18",
    "x19",
    "x20",
    "x21",
    "x22",
    "x23",
    "x24",
    "x25",
    "x26",
    "x27",
    "x28",
    "x29",
    "x30",
]

RESERVED_VECTOR_REGS = [
    "v0",
]

WORK_REGS = [f"w{i:02d}" for i in range(32)]

REGION_OUTPUTS = {
    "_ntt32_stage12_slothy_start": [f"{reg}_s2" for reg in WORK_REGS],
    "_ntt32_stage3_slothy_start": [f"{reg}_s3" for reg in WORK_REGS],
    "_ntt32_stage4_slothy_start": [f"{reg}_s4" for reg in WORK_REGS],
    "_ntt32_stage5_store_slothy_start": [],
}


def find_regions(source_path: Path) -> list[tuple[str, str]]:
    regions: list[tuple[str, str]] = []
    pending_start: tuple[str, str] | None = None

    for line_number, line in enumerate(source_path.read_text(encoding="utf8").splitlines(), 1):
        start_match = REGION_START_RE.match(line)
        if start_match:
            if pending_start is not None:
                _, label = pending_start
                raise ValueError(f"Found {label}_slothy_start before its end label")
            label = start_match.group(1)
            pending_start = (str(line_number), label)
            continue

        end_match = REGION_END_RE.match(line)
        if end_match:
            label = end_match.group(1)
            if pending_start is None:
                raise ValueError(f"Found {label}_slothy_end without a start label")
            _, start_label = pending_start
            if start_label != label:
                raise ValueError(
                    f"Mismatched Slothy labels: {start_label}_slothy_start "
                    f"ends at {label}_slothy_end"
                )
            regions.append((f"{label}_slothy_start", f"{label}_slothy_end"))
            pending_start = None

    if pending_start is not None:
        _, label = pending_start
        raise ValueError(f"Found {label}_slothy_start without an end label")

    if not regions:
        raise ValueError(f"No Slothy regions found in {source_path}")

    return regions


def new_slothy() -> Slothy:
    slothy = Slothy(AArch64_Neon, Target_CortexA72)
    slothy.config.logger = slothy.logger
    slothy.config.with_llvm_mca = False
    slothy.config.selftest = False
    slothy.config.reserved_regs = (
        list(slothy.config.reserved_regs) + RESERVED_GPRS + RESERVED_VECTOR_REGS
    )
    return slothy


def configure_register_allocation(slothy: Slothy) -> None:
    slothy.config.allow_useless_instructions = True
    slothy.config.inputs_are_outputs = False
    slothy.config.constraints.functional_only = True
    slothy.config.constraints.allow_reordering = True
    slothy.config.constraints.allow_spills = True


def configure_window_optimization(slothy: Slothy) -> None:
    slothy.config.allow_useless_instructions = True
    slothy.config.constraints.allow_spills = True
    slothy.config.constraints.allow_reordering = True
    slothy.config.constraints.functional_only = False
    slothy.config.variable_size = True
    slothy.config.constraints.stalls_first_attempt = 32
    slothy.config.split_heuristic = True
    slothy.config.split_heuristic_stepsize = 0.025
    slothy.config.split_heuristic_factor = 7.5
    slothy.config.split_heuristic_repeat = 1
    slothy.config.split_heuristic_estimate_performance = False


def optimize_regions(slothy: Slothy, regions: list[tuple[str, str]]) -> None:
    for start, end in regions:
        slothy.config.outputs = REGION_OUTPUTS.get(start, [])
        print(f"  {start} -> {end}")
        slothy.optimize(start=start, end=end)


def run_register_allocation(regions: list[tuple[str, str]]) -> None:
    slothy = new_slothy()
    slothy.load_source_from_file(str(SOURCE_PATH))
    configure_register_allocation(slothy)
    optimize_regions(slothy, regions)
    slothy.write_source_to_file(str(ALLOC_PATH))


def run_window_optimization(regions: list[tuple[str, str]]) -> None:
    slothy = new_slothy()
    slothy.load_source_from_file(str(ALLOC_PATH))
    configure_window_optimization(slothy)
    optimize_regions(slothy, regions)
    slothy.write_source_to_file(str(OUTPUT_PATH))


def run_timed(step_name: str, fn) -> None:
    print(step_name)
    start = time.perf_counter()
    fn()
    elapsed = time.perf_counter() - start
    print(f"{step_name}: {elapsed:.6f}s")


def main() -> None:
    regions = find_regions(SOURCE_PATH)
    print(f"Found {len(regions)} Slothy regions in {SOURCE_PATH}")

    # This script intentionally loads the full source for every pass. If a
    # future workflow extracts one label body into a standalone temporary file,
    # also carry the preceding "adr tw_ptr, ..." line or initialize x12 in the
    # wrapper before entering the region.
    run_timed("pass1_register_allocation", lambda: run_register_allocation(regions))
    run_timed("pass2_window_opt", lambda: run_window_optimization(regions))
    print(ALLOC_PATH)
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
