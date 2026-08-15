import re
import tomllib
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).parents[2]
MONOIMAGE_PACKAGES = {"biopython", "blast", "deacon", "meryl", "polars", "python"}
MODULE_DEPENDENCY = re.compile(
    r"^\s*-\s+(?:[^:]+::)?(?P<name>[a-z0-9_-]+)=(?P<version>[^=\s]+)",
    re.MULTILINE,
)
MANIFEST_VERSION = re.compile(
    r"^manifest\s*\{.*?^\s+version\s*=\s*['\"](?P<version>[^'\"]+)['\"]",
    re.MULTILINE | re.DOTALL,
)


def _workspace_pins() -> dict[str, str]:
    with (REPOSITORY_ROOT / "pyproject.toml").open("rb") as pyproject_file:
        pyproject = tomllib.load(pyproject_file)

    pins = {
        name: version.removeprefix("==").removesuffix(".*")
        for name, version in pyproject["tool"]["pixi"]["dependencies"].items()
        if name in MONOIMAGE_PACKAGES
    }
    for requirement in pyproject["project"]["dependencies"]:
        name, separator, version = requirement.partition("==")
        if name in MONOIMAGE_PACKAGES:
            assert separator == "==", f"{name} must use an exact workspace pin"
            pins[name] = version
    return pins


def _module_pins() -> dict[str, set[str]]:
    pins: dict[str, set[str]] = {}
    for environment_file in (REPOSITORY_ROOT / "modules").glob("**/environment.yml"):
        for match in MODULE_DEPENDENCY.finditer(environment_file.read_text()):
            name = match.group("name")
            if name in MONOIMAGE_PACKAGES:
                pins.setdefault(name, set()).add(match.group("version"))
    return pins


def test_module_dependency_versions_match_monoimage_workspace() -> None:
    workspace_pins = _workspace_pins()
    module_pins = _module_pins()

    assert module_pins.keys() == workspace_pins.keys()
    assert module_pins == {
        name: {version} for name, version in workspace_pins.items()
    }


def test_project_versions_match() -> None:
    with (REPOSITORY_ROOT / "pyproject.toml").open("rb") as pyproject_file:
        project_version = tomllib.load(pyproject_file)["project"]["version"]

    manifest_match = MANIFEST_VERSION.search(
        (REPOSITORY_ROOT / "nextflow.config").read_text(),
    )
    assert manifest_match is not None

    with (REPOSITORY_ROOT / "uv.lock").open("rb") as lock_file:
        packages = tomllib.load(lock_file)["package"]
    locked_version = next(
        package["version"]
        for package in packages
        if package["name"] == "baits"
        and package.get("source", {}).get("virtual") == "."
    )

    assert {
        project_version,
        manifest_match.group("version"),
        locked_version,
    } == {project_version}
