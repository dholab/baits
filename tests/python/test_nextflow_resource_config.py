import re
import subprocess
import sys
from pathlib import Path

NEXTFLOW = Path(sys.executable).with_name("nextflow")
REPOSITORY_ROOT = Path(__file__).parents[2]
PROCESS_LABEL = re.compile(r"^\s*label\s+['\"](?P<label>process_[^'\"]+)['\"]", re.MULTILINE)


def resolved_config(*global_options: str) -> dict[str, str]:
    completed = subprocess.run(
        [NEXTFLOW, *global_options, "config", "-flat"],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        key: value
        for line in completed.stdout.splitlines()
        if " = " in line
        for key, value in [line.split(" = ", maxsplit=1)]
    }


def declared_process_labels() -> set[str]:
    return {
        match.group("label")
        for module in (REPOSITORY_ROOT / "modules").glob("**/*.nf")
        for match in PROCESS_LABEL.finditer(module.read_text())
    }


def test_every_process_resource_label_has_effective_portable_defaults() -> None:
    config = resolved_config()
    expected = {
        "process_low": {"cpus": "1", "memory": "'2 GB'", "time": "'2h'"},
        "process_medium": {"cpus": "2", "memory": "'4 GB'", "time": "'8h'"},
        "process_high": {"cpus": "2", "memory": "'4 GB'", "time": "'24h'"},
    }

    assert declared_process_labels() == expected.keys()
    assert {
        label: {
            resource: config[f"process.'withLabel:{label}'.{resource}"]
            for resource in resources
        }
        for label, resources in expected.items()
    } == expected
    assert config["process.'withName:BLAST_MAKEBLASTDB'.cpus"] == "1"


def test_site_config_can_override_portable_resource_defaults(tmp_path: Path) -> None:
    site_config = tmp_path / "site.config"
    site_config.write_text(
        "process {\n"
        "    withLabel: process_medium {\n"
        "        cpus = 7\n"
        "        memory = '11 GB'\n"
        "        time = '1h'\n"
        "    }\n"
        "}\n",
    )

    config = resolved_config("-c", str(site_config))

    assert config["process.'withLabel:process_medium'.cpus"] == "7"
    assert config["process.'withLabel:process_medium'.memory"] == "'11 GB'"
    assert config["process.'withLabel:process_medium'.time"] == "'1h'"
