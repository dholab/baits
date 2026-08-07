#!/usr/bin/env python3
"""Write deterministic per-design provenance tables."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence

INPUT_SCHEMA = {
    "input_role": pl.String,
    "input_id": pl.String,
    "attribute": pl.String,
    "value": pl.String,
}
PARAMETER_SCHEMA = {"parameter": pl.String, "value": pl.String}
SOFTWARE_SCHEMA = {"component": pl.String, "version": pl.String}


class ProvenanceInputError(ValueError):
    """Raised when provenance input cannot be parsed."""



@dataclass(frozen=True)
class InputFile:
    """A staged filesystem object contributing input provenance."""

    role: str
    identifier: str
    kind: str
    path: Path


@dataclass(frozen=True)
class ProvenanceRequest:
    """The complete facts and staged objects for one provenance bundle."""

    input_facts_json: str
    input_file_roles: Sequence[str]
    input_file_ids: Sequence[str]
    input_file_kinds: Sequence[str]
    input_file_paths: Sequence[Path]
    parameters_json: str
    software_versions_json: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-facts-base64", required=True)
    parser.add_argument("--input-file-roles-base64", required=True)
    parser.add_argument("--input-file-ids-base64", required=True)
    parser.add_argument("--input-file-kinds-base64", required=True)
    input_paths = parser.add_mutually_exclusive_group(required=True)
    input_paths.add_argument("--input-files", type=Path, nargs="+")
    input_paths.add_argument("--input-root", type=Path)
    parser.add_argument("--parameters-base64", required=True)
    parser.add_argument("--software-versions-base64", required=True)
    parser.add_argument("--inputs-out", type=Path, required=True)
    parser.add_argument("--parameters-out", type=Path, required=True)
    parser.add_argument("--software-versions-out", type=Path, required=True)
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_tsv_value(value: str) -> bool:
    return value.strip() != "" and not any(character in value for character in "\t\n\r")


def decode_base64(value: str, flag: str) -> str:
    try:
        return base64.b64decode(value, validate=True).decode()
    except (ValueError, UnicodeDecodeError) as error:
        message = f"{flag} must contain Base64-encoded UTF-8"
        raise ProvenanceInputError(message) from error


def decode_string_list(value: str, flag: str) -> tuple[str, ...]:
    raw = decode_base64(value, flag)
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError as error:
        message = f"{flag} must encode a JSON list: {error}"
        raise ProvenanceInputError(message) from error
    if not isinstance(decoded, list) or any(not isinstance(item, str) for item in decoded):
        message = f"{flag} must encode a JSON list of strings"
        raise ProvenanceInputError(message)
    return tuple(decoded)


def parse_json_rows(raw: str, fields: tuple[str, ...], flag: str) -> tuple[tuple[str, ...], ...]:
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError as error:
        message = f"{flag} must be valid JSON: {error}"
        raise ProvenanceInputError(message) from error
    if not isinstance(decoded, list):
        message = f"{flag} must be a JSON list"
        raise ProvenanceInputError(message)
    rows: list[tuple[str, ...]] = []
    for index, record in enumerate(decoded):
        if not isinstance(record, dict) or set(record) != set(fields):
            message = f"{flag} item {index} must have exactly these fields: {', '.join(fields)}"
            raise ProvenanceInputError(message)
        values = tuple(record[field] for field in fields)
        if any(not isinstance(value, str) or not is_tsv_value(value) for value in values):
            message = f"{flag} item {index} must contain only nonblank single-line TSV values"
            raise ProvenanceInputError(message)
        rows.append(values)
    return tuple(rows)


def input_files(
    roles: Sequence[str],
    identifiers: Sequence[str],
    kinds: Sequence[str],
    paths: Sequence[Path],
) -> tuple[InputFile, ...]:
    try:
        descriptors = tuple(
            InputFile(role, identifier, kind, path)
            for role, identifier, kind, path in zip(roles, identifiers, kinds, paths, strict=True)
        )
    except ValueError as error:
        message = "--input-file-roles, --input-file-ids, --input-file-kinds, and --input-files must have equal lengths"
        raise ProvenanceInputError(message) from error
    if not descriptors:
        message = "input-file descriptors must be nonempty"
        raise ProvenanceInputError(message)
    for descriptor in descriptors:
        if not is_tsv_value(descriptor.role) or not is_tsv_value(descriptor.identifier):
            message = "--input-file-roles and --input-file-ids must contain only nonblank single-line TSV values"
            raise ProvenanceInputError(message)
        if descriptor.kind not in {"file", "directory"}:
            message = "--input-file-kinds must contain only file or directory"
            raise ProvenanceInputError(message)
        if descriptor.kind == "file" and not descriptor.path.is_file():
            message = f"--input-files path must be a regular file: {descriptor.path}"
            raise ProvenanceInputError(message)
        if descriptor.kind == "directory" and not descriptor.path.is_dir():
            message = f"--input-files path must be a directory: {descriptor.path}"
            raise ProvenanceInputError(message)
        if not is_tsv_value(descriptor.path.name):
            message = f"--input-files path must have a nonblank single-line TSV basename: {descriptor.path}"
            raise ProvenanceInputError(message)
    return descriptors


def staged_input_paths(root: Path) -> tuple[Path, ...]:
    try:
        buckets = tuple(sorted(path for path in root.iterdir() if path.is_dir()))
        bucket_paths = tuple(tuple(bucket.iterdir()) for bucket in buckets)
    except OSError as error:
        message = f"--input-root cannot be read: {root}"
        raise ProvenanceInputError(message) from error
    if not bucket_paths or any(len(paths) != 1 for paths in bucket_paths):
        message = f"--input-root must contain one staged path in each input directory: {root}"
        raise ProvenanceInputError(message)
    return tuple(paths[0] for paths in bucket_paths)


def reject_duplicate_keys(frame: pl.LazyFrame, keys: tuple[str, ...], label: str) -> None:
    duplicates = frame.group_by(*keys).len().filter(pl.col("len") > 1).limit(1).collect()
    if duplicates.height:
        message = f"{label} contains a duplicate key"
        raise ProvenanceInputError(message)


def construct_provenance(request: ProvenanceRequest) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    input_facts = parse_json_rows(
        request.input_facts_json,
        ("input_role", "input_id", "attribute", "value"),
        "--input-facts-json",
    )
    parameters = parse_json_rows(request.parameters_json, ("parameter", "value"), "--parameters-json")
    versions = parse_json_rows(
        request.software_versions_json,
        ("component", "version"),
        "--software-versions-json",
    )
    descriptors = input_files(
        request.input_file_roles,
        request.input_file_ids,
        request.input_file_kinds,
        request.input_file_paths,
    )
    file_rows = tuple(
        row
        for descriptor in descriptors
        for row in (
            ((descriptor.role, descriptor.identifier, "directory", descriptor.path.name),)
            if descriptor.kind == "directory"
            else (
                (descriptor.role, descriptor.identifier, "filename", descriptor.path.name),
                (descriptor.role, descriptor.identifier, "sha256", sha256_file(descriptor.path)),
            )
        )
    )
    input_frame = pl.LazyFrame(input_facts + file_rows, schema=INPUT_SCHEMA, orient="row")
    parameter_frame = pl.LazyFrame(parameters, schema=PARAMETER_SCHEMA, orient="row")
    version_frame = pl.LazyFrame(
        (*versions, ("python", platform.python_version())),
        schema=SOFTWARE_SCHEMA,
        orient="row",
    ).unique()
    reject_duplicate_keys(input_frame, ("input_role", "input_id", "attribute"), "input facts")
    reject_duplicate_keys(parameter_frame, ("parameter",), "parameters")
    reject_duplicate_keys(version_frame, ("component",), "software versions")
    return (
        input_frame.sort("input_role", "input_id", "attribute").collect(),
        parameter_frame.sort("parameter").collect(),
        version_frame.sort("component").collect(),
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    input_paths = staged_input_paths(args.input_root) if args.input_root is not None else args.input_files
    inputs, parameters, versions = construct_provenance(ProvenanceRequest(
        decode_base64(args.input_facts_base64, "--input-facts-base64"),
        decode_string_list(args.input_file_roles_base64, "--input-file-roles-base64"),
        decode_string_list(args.input_file_ids_base64, "--input-file-ids-base64"),
        decode_string_list(args.input_file_kinds_base64, "--input-file-kinds-base64"),
        input_paths,
        decode_base64(args.parameters_base64, "--parameters-base64"),
        decode_base64(args.software_versions_base64, "--software-versions-base64"),
    ))
    inputs.write_csv(args.inputs_out, separator="\t", quote_style="never")
    parameters.write_csv(args.parameters_out, separator="\t", quote_style="never")
    versions.write_csv(args.software_versions_out, separator="\t", quote_style="never")


if __name__ == "__main__":
    main()
