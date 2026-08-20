import re
import subprocess
import sys
from pathlib import Path

NEXTFLOW = Path(sys.executable).with_name("nextflow")
REPOSITORY_ROOT = Path(__file__).parents[2]
PROCESS_LABEL = re.compile(r"^\s*label\s+['\"](?P<label>process_[^'\"]+)['\"]", re.MULTILINE)


def resolved_config(*global_options: str, profile: str | None = None) -> dict[str, str]:
    command = [NEXTFLOW, *global_options, "config"]
    if profile is not None:
        command.extend(("-profile", profile))
    command.append("-flat")
    completed = subprocess.run(
        command,
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


def selected_resources(
    config: dict[str, str],
    expected: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    return {
        label: {
            resource: config[f"process.'withLabel:{label}'.{resource}"]
            for resource in resources
        }
        for label, resources in expected.items()
    }


def test_every_process_resource_label_has_effective_production_defaults() -> None:
    config = resolved_config()
    expected = {
        "process_low": {
            "cpus": "1",
            "memory": "'4 GB'",
            "disk": "'20 GB'",
            "time": "'4h'",
        },
        "process_medium": {
            "cpus": "6",
            "memory": "'16 GB'",
            "disk": "'100 GB'",
            "time": "'12h'",
        },
        "process_high": {
            "cpus": "20",
            "memory": "'64 GB'",
            "disk": "'300 GB'",
            "time": "'24h'",
        },
        "process_extra_high": {
            "cpus": "20",
            "memory": "'64 GB'",
            "disk": "'500 GB'",
            "time": "'48h'",
        },
    }

    assert declared_process_labels() == expected.keys()
    assert selected_resources(config, expected) == expected
    assert config["process.'withName:BLAST_MAKEBLASTDB'.cpus"] == "1"


def test_test_profile_fits_the_ci_executor() -> None:
    config = resolved_config(profile="test")
    expected = {
        "process_low": {
            "cpus": "1",
            "memory": "'2 GB'",
            "disk": "'5 GB'",
            "time": "'2h'",
        },
        "process_medium": {
            "cpus": "2",
            "memory": "'4 GB'",
            "disk": "'5 GB'",
            "time": "'8h'",
        },
        "process_high": {
            "cpus": "2",
            "memory": "'4 GB'",
            "disk": "'5 GB'",
            "time": "'12h'",
        },
        "process_extra_high": {
            "cpus": "2",
            "memory": "'4 GB'",
            "disk": "'5 GB'",
            "time": "'24h'",
        },
    }

    assert selected_resources(config, expected) == expected
    assert config["process.executor"] == "'local'"
    assert config["executor.cpus"] == "2"
    assert config["executor.memory"] == "'4 GB'"


def test_site_config_can_override_production_resource_defaults(tmp_path: Path) -> None:
    site_config = tmp_path / "site.config"
    site_config.write_text(
        "process {\n"
        "    withLabel: process_extra_high {\n"
        "        cpus = 7\n"
        "        memory = '11 GB'\n"
        "        disk = '17 GB'\n"
        "        time = '1h'\n"
        "    }\n"
        "}\n",
    )

    config = resolved_config("-c", str(site_config))

    assert config["process.'withLabel:process_extra_high'.cpus"] == "7"
    assert config["process.'withLabel:process_extra_high'.memory"] == "'11 GB'"
    assert config["process.'withLabel:process_extra_high'.disk"] == "'17 GB'"
    assert config["process.'withLabel:process_extra_high'.time"] == "'1h'"


def test_taxonomic_reference_searches_use_the_extra_high_tier() -> None:
    taxonomic_screening = (
        REPOSITORY_ROOT / "modules/local/blastn_taxonomic_screening/main.nf"
    ).read_text()
    read_classification = (
        REPOSITORY_ROOT / "modules/local/blastn_read_classification/main.nf"
    ).read_text()
    source_discovery = (
        REPOSITORY_ROOT / "modules/local/blastn_source_sequences/main.nf"
    ).read_text()

    assert "label 'process_extra_high'" in taxonomic_screening
    assert "label 'process_extra_high'" in read_classification
    assert "label 'process_medium'" in source_discovery
