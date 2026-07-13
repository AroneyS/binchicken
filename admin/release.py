#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    print(f"+ {' '.join(cmd)}")
    return subprocess.run(cmd, check=check)


def capture(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True).strip()


def die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    sys.exit(1)


def validate_version(version: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?", version):
        die(
            f"invalid version {version!r}; expected something like "
            "1.2.3, 1.2.3-beta.1, or 1.2.3+build.1"
        )


def check_clean_git() -> None:
    status = capture(["git", "status", "--porcelain"])
    if status:
        die("git working tree is not clean; commit or stash changes first")


def get_current_version() -> str:
    init_path = Path("binchicken/__init__.py")
    for line in init_path.read_text().splitlines():
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]
    die(f"could not find __version__ in {init_path}")


def update_init_version(version: str) -> None:
    init_path = Path("binchicken/__init__.py")
    text = init_path.read_text()
    new_text, count = re.subn(
        r'(?m)^__version__\s*=\s*["\'][^"\']+["\']',
        f'__version__ = "{version}"',
        text,
        count=1,
    )
    if count != 1:
        die(f"could not find exactly one __version__ line in {init_path}")
    init_path.write_text(new_text)
    print(f"Updated {init_path}: __version__ = {version!r}")


def changelog_contains_version(version: str, path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text()
    return bool(re.search(rf"(?m)^## \[{re.escape(version)}\]", text))


def update_changelog(version: str, path: Path) -> None:
    if not path.exists():
        die(f"{path} not found — create it with an '## [Unreleased]' section first")

    text = path.read_text()

    if changelog_contains_version(version, path):
        print(f"{path} already contains a section for {version}; skipping.")
        return

    unreleased_match = re.search(r"(?m)^## \[Unreleased\]\s*$", text)
    if not unreleased_match:
        die(f"could not find '## [Unreleased]' section in {path}")

    next_heading_match = re.search(r"(?m)^## \[", text[unreleased_match.end():])
    if not next_heading_match:
        die(f"could not find next '## [' release section after Unreleased in {path}")

    unreleased_start = unreleased_match.end()
    next_heading_start = unreleased_start + next_heading_match.start()
    unreleased_body = text[unreleased_start:next_heading_start].strip()

    if not unreleased_body:
        die("CHANGELOG.md [Unreleased] section is empty; add entries before releasing")

    date_str = datetime.today().strftime("%Y-%m-%d")
    new_section = f"\n\n## [{version}] - {date_str}\n\n{unreleased_body}\n\n"

    before = text[:unreleased_start].rstrip()
    after = text[next_heading_start:].lstrip()
    path.write_text(before + new_section + after)
    print(f"Moved [Unreleased] entries to [{version}] - {date_str} in {path}")


def update_citation(version: str) -> None:
    citation_path = Path("CITATION.cff")
    if not citation_path.exists():
        print("No CITATION.cff found; skipping.")
        return
    lines = []
    r_version = re.compile(r"( *version: )")
    r_date = re.compile(r"( *date-released: )")
    for line in citation_path.read_text().splitlines(keepends=True):
        if m := r_version.match(line):
            line = m.group(1) + version + "\n"
        elif m := r_date.match(line):
            line = m.group(1) + datetime.today().strftime("%Y-%m-%d") + "\n"
        lines.append(line)
    citation_path.write_text("".join(lines))
    print(f"Updated {citation_path}: version = {version!r}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a binchicken release."
    )
    parser.add_argument(
        "--version",
        help="Release version, e.g. 0.14.0 (default: read from binchicken/__init__.py)",
    )
    parser.add_argument(
        "--no-commit",
        action="store_true",
        help="Modify files but do not commit, tag, or push",
    )
    parser.add_argument(
        "--tag-prefix",
        default="v",
        help='Git tag prefix; default is "v", producing tags like v0.14.0',
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow running with an already-dirty git working tree",
    )

    args = parser.parse_args()

    version = args.version or get_current_version()
    validate_version(version)
    tag = f"{args.tag_prefix}{version}"

    if not args.allow_dirty:
        check_clean_git()

    existing_tags = capture(["git", "tag", "--list", tag])
    if existing_tags:
        die(f"git tag {tag!r} already exists")

    yes_no = input(
        f"\nHave you run the full test suite (pixi run run-tests-at-cmr) and confirmed it passes? [y/n] "
    ).strip().lower()
    if yes_no not in {"y", "yes"}:
        die("run the full test suite first")

    # Update version in __init__.py if a new version was specified
    current_version = get_current_version()
    if version != current_version:
        update_init_version(version)

    # Move Unreleased → versioned section in CHANGELOG.md
    update_changelog(version, Path("CHANGELOG.md"))

    # Sync dependencies from pixi to pyproject / requirements.txt
    print("\nSyncing dependencies from pixi.lock/pixi.toml ...")
    run(["pixi", "run", "-e", "dev", "build_dep_defs"])

    # Update CITATION.cff
    update_citation(version)

    # Build docs
    print("\nBuilding docs ...")
    run(["pixi", "run", "-e", "dev", "python3", "admin/build_docs.py", "--version", version])

    # Check still clean (docs may have changed)
    status = capture(["git", "status", "--porcelain"])
    changed_files = [
        f for f in [
            "binchicken/__init__.py",
            "CHANGELOG.md",
            "CITATION.cff",
            "binchicken.yml",
            "admin/requirements.txt",
        ]
        if Path(f).exists()
    ]

    if args.no_commit:
        print("\nStopped before commit/tag because --no-commit was supplied.")
        print("Review changes with: git diff")
        return

    print(f"\nTagging the release as {tag}")
    run(["git", "add"] + changed_files + ["docs/"])
    run(["git", "commit", "-m", f"Release {tag}"])
    run(["git", "tag", tag])
    print(f"\nNow run: git push && git push --tags")
    print("GitHub Actions will create the GitHub release (from CHANGELOG) and publish to PyPI.")
    print("Then update the bioconda recipe.")
    print("Once the new version is on conda, run the 'binchicken-installation' GitHub workflow.")


if __name__ == "__main__":
    main()
