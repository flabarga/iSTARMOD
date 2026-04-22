import argparse
import sys
from pathlib import Path

from iStarmod_automat_reloaded import starmod


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="istarmod",
        description="Run iSTARMOD from a .sm configuration file."
    )

    parser.add_argument(
        "config",
        help="Path to the iSTARMOD .sm configuration file"
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debugging mode"
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="Enable plotting of results"
    )

    args = parser.parse_args()
    config_path = Path(args.config)

    if not config_path.exists():
        print(f"Error: configuration file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        starmod(
            str(config_path),
            debugging=args.debug,
            plot=args.plot
        )
    except Exception as exc:
        print(f"iSTARMOD failed: {exc}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
