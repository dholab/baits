from pathlib import Path

import apply_complexity_filter as complexity
import build_candidate_tables as candidates
import polars as pl
import pytest

KMER_A = "ACGTTGCATGTCAGTACGATCGTAGCTAGCA"
KMER_B = "GTTCCGAATCGTACGATCAGGTCATGGCTAA"
KMER_B_RC = "TTAGCCATGACCTGATCGTACGATTCGGAAC"


def sources(rows: dict[str, list[str]]) -> pl.LazyFrame:
    return pl.LazyFrame(rows, schema={"source_sequence_id": pl.String, "sequence": pl.String})


def groups(rows: dict[str, list[str]]) -> pl.LazyFrame:
    return pl.LazyFrame(rows, schema=candidates.QUERY_GROUP_SCHEMA)


def counts(rows: dict[str, list[object]]) -> pl.LazyFrame:
    return pl.LazyFrame(rows, schema=candidates.MERYL_COUNT_SCHEMA)


def manifest(rows: list[dict[str, str]]) -> pl.LazyFrame:
    return pl.LazyFrame(rows, schema=complexity.MANIFEST_SCHEMA)


def valid_manifest_rows() -> list[dict[str, str]]:
    return [
        dict(zip(complexity.MANIFEST_FIELDS, values, strict=True))
        for values in (
            ("candidate_kmer_000001", "", KMER_A, "1", "0", "PASS", "none", "NOT_RUN", "", ""),
            ("candidate_kmer_000002", "", KMER_B, "2", "4", "REJECT_INTERFERENCE_BACKGROUND", "background_occurrence", "NOT_APPLICABLE", "", ""),
        )
    ]


def write_build_inputs(tmp_path: Path, sequence: str, *, source_counts: str, background: str = "") -> dict[str, Path]:
    tmp_path.mkdir(exist_ok=True)
    paths = {name: tmp_path / f"{name}.{'fasta' if name == 'source' else 'tsv'}" for name in ("source", "groups", "source_counts", "background")}
    paths["source"].write_text(f">source\n{sequence}\n")
    paths["groups"].write_text("source_sequence_id\tquery_group\nsource\t\n")
    paths["source_counts"].write_text(source_counts)
    paths["background"].write_text(background)
    return paths


def test_python_enumeration_canonicalizes_and_preserves_all_overlapping_zero_based_occurrences() -> None:
    evidence = candidates.construct_candidate_evidence(
        sources({"source_sequence_id": ["a", "b"], "sequence": ["ACGTAC", "ACGT"]}),
        groups({"source_sequence_id": ["a", "b"], "query_group": ["alpha", "beta"]}),
        counts({"kmer": ["ACGT", "CGTA", "GTAC"], "count": [2, 1, 1]}),
        counts({"kmer": [], "count": []}),
        4,
    )

    assert evidence.manifest.collect().filter(pl.col("kmer") == "ACGT").select("source_copy_count").rows() == [(2,)]
    assert evidence.occurrences.collect().select("source_sequence_id", "start", "query_group").sort("source_sequence_id", "start").rows() == [
        ("a", 0, "alpha"),
        ("a", 1, "alpha"),
        ("a", 2, "alpha"),
        ("b", 0, "beta"),
    ]


def test_candidate_evidence_is_stable_under_source_and_query_group_permutation() -> None:
    results = []
    for ids in (("a", "b"), ("b", "a")):
        evidence = candidates.construct_candidate_evidence(
            sources({"source_sequence_id": list(ids), "sequence": [KMER_A if key == "a" else KMER_B for key in ids]}),
            groups({"source_sequence_id": list(reversed(ids)), "query_group": ["beta" if key == "b" else "alpha" for key in reversed(ids)]}),
            counts({"kmer": [KMER_A, KMER_B], "count": [1, 1]}),
            counts({"kmer": [], "count": []}),
            len(KMER_A),
        )
        results.append((evidence.manifest.collect(), evidence.occurrences.collect()))
    assert results[0][0].equals(results[1][0])
    assert results[0][1].equals(results[1][1])


def test_source_reader_allows_iupac_but_enumeration_skips_affected_windows_and_rejects_malformed_alphabet(tmp_path: Path) -> None:
    paths = write_build_inputs(tmp_path, "NN" + KMER_A, source_counts=f"{KMER_A}\t1\n")
    source = candidates.read_source_sequences(paths["source"])
    evidence = candidates.construct_candidate_evidence(source, candidates.scan_query_groups(paths["groups"], source), candidates.scan_meryl_counts(paths["source_counts"], len(KMER_A), "source"), counts({"kmer": [], "count": []}), len(KMER_A))
    assert evidence.occurrences.collect().select("start").rows() == [(2,)]

    paths["source"].write_text(">source\nACGTZ\n")
    with pytest.raises(candidates.SourceSequenceError, match="malformed DNA/IUPAC"):
        candidates.read_source_sequences(paths["source"])
    paths["source"].write_text(">source\nacgt\n")
    with pytest.raises(candidates.SourceSequenceError, match="uppercase bases"):
        candidates.read_source_sequences(paths["source"])


def test_candidate_enumeration_does_not_underflow_for_source_shorter_than_kmer() -> None:
    evidence = candidates.construct_candidate_evidence(
        sources({"source_sequence_id": ["short"], "sequence": ["ACGT"]}),
        groups({"source_sequence_id": ["short"], "query_group": [""]}),
        counts({"kmer": [], "count": []}),
        counts({"kmer": [], "count": []}),
        31,
    )

    assert evidence.manifest.collect().is_empty()
    assert evidence.occurrences.collect().is_empty()


def test_query_groups_preserve_blank_labels_and_require_exact_schema_and_id_relation(tmp_path: Path) -> None:
    source = sources({"source_sequence_id": ["a", "b"], "sequence": [KMER_A, KMER_B]})
    path = tmp_path / "groups.tsv"
    path.write_text("source_sequence_id\tquery_group\na\t\nb\tgroup\n")
    assert candidates.scan_query_groups(path, source).collect().rows() == [("a", ""), ("b", "group")]
    for content, error in (("query_group\tsource_sequence_id\n\ta\n", "columns must be exactly"), ("source_sequence_id\tquery_group\na\tg\n", "Missing Source Sequence Query Group ID"), ("source_sequence_id\tquery_group\na\tg\nb\tg\nc\tg\n", "Unknown Source Sequence Query Group ID")):
        path.write_text(content)
        with pytest.raises(candidates.QueryGroupError, match=error):
            candidates.scan_query_groups(path, source)


@pytest.mark.parametrize(("content", "error"), [("", None), ("AAAA\t1\textra\n", "Could not parse"), ("AAA\t1\n", "malformed k-mer"), ("AAAA\t0\n", "nonpositive count")])
def test_meryl_count_parser_requires_strict_width_valid_kmers_and_positive_counts(tmp_path: Path, content: str, error: str | None) -> None:
    path = tmp_path / "meryl.tsv"
    path.write_text(content)
    if error is None:
        assert candidates.scan_meryl_counts(path, 4, "Meryl evidence").collect().is_empty()
    else:
        with pytest.raises(candidates.MerylCountError, match=error):
            candidates.scan_meryl_counts(path, 4, "Meryl evidence")


def test_meryl_count_parser_canonicalizes_orientation_and_rejects_canonical_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "meryl.tsv"
    path.write_text(f"{KMER_B_RC}\t2\n")
    assert candidates.scan_meryl_counts(path, len(KMER_B), "Meryl evidence").collect().rows() == [(KMER_B, 2)]
    path.write_text(f"{KMER_B}\t1\n{KMER_B_RC}\t1\n")
    with pytest.raises(candidates.MerylCountError, match="Duplicate canonical k-mer"):
        candidates.scan_meryl_counts(path, len(KMER_B), "Meryl evidence")


def test_meryl_source_keys_and_positional_counts_must_match_python_candidates_exactly() -> None:
    source = sources({"source_sequence_id": ["a"], "sequence": [KMER_A + "A"]})
    query_groups = groups({"source_sequence_id": ["a"], "query_group": ["group"]})
    with pytest.raises(candidates.CandidateEvidenceError, match="keys and positional occurrence counts"):
        candidates.construct_candidate_evidence(source, query_groups, counts({"kmer": [KMER_A], "count": [1]}), counts({"kmer": [], "count": []}), len(KMER_A))


def test_background_evidence_is_candidate_subset_and_anti_join_uses_background_first_multiplicity() -> None:
    source = sources({"source_sequence_id": ["a"], "sequence": [KMER_A + "A"]})
    query_groups = groups({"source_sequence_id": ["a"], "query_group": ["group"]})
    evidence = candidates.construct_candidate_evidence(source, query_groups, counts({"kmer": [KMER_A, candidates.project_canonical(KMER_A[1:] + "A")], "count": [1, 1]}), counts({"kmer": [KMER_A], "count": [7]}), len(KMER_A))
    assert evidence.manifest.collect().select("kmer", "background_occurrences", "status").rows() == [
        (KMER_A, 7, "REJECT_INTERFERENCE_BACKGROUND"),
        (candidates.project_canonical(KMER_A[1:] + "A"), 0, "PASS"),
    ]
    assert evidence.complexity_candidates.collect().rows() == [("candidate_kmer_000002", candidates.project_canonical(KMER_A[1:] + "A"))]
    with pytest.raises(candidates.CandidateEvidenceError, match="non-Candidate"):
        candidates.construct_candidate_evidence(source, query_groups, counts({"kmer": [KMER_A, candidates.project_canonical(KMER_A[1:] + "A")], "count": [1, 1]}), counts({"kmer": [KMER_B], "count": [1]}), len(KMER_A))


@pytest.mark.parametrize("mutate", [lambda row: row.__setitem__("status", "PASS"), lambda row: row.__setitem__("bait_id", "bait_000001")])
def test_pre_complexity_manifest_requires_exact_states_and_empty_evidence(tmp_path: Path, mutate: object) -> None:
    rows = valid_manifest_rows()
    mutate(rows[1])  # type: ignore[operator]
    path = tmp_path / "manifest.tsv"
    pl.DataFrame(rows, schema=complexity.MANIFEST_SCHEMA).write_csv(path, separator="\t")
    with pytest.raises(complexity.ComplexityFilterError, match="pre-complexity"):
        complexity.scan_manifest(path)


def test_pre_complexity_manifest_requires_sorted_kmers_and_an_eligible_candidate(tmp_path: Path) -> None:
    rows = valid_manifest_rows()
    rows[0]["kmer"], rows[1]["kmer"] = rows[1]["kmer"], rows[0]["kmer"]
    path = tmp_path / "unsorted.tsv"
    pl.DataFrame(rows, schema=complexity.MANIFEST_SCHEMA).write_csv(path, separator="\t")
    with pytest.raises(complexity.ComplexityFilterError, match="not sorted by k-mer"):
        complexity.scan_manifest(path)

    rejected = valid_manifest_rows()[1:]
    rejected[0]["candidate_kmer_id"] = "candidate_kmer_000001"
    path = tmp_path / "no-eligible.tsv"
    pl.DataFrame(rejected, schema=complexity.MANIFEST_SCHEMA).write_csv(path, separator="\t")
    with pytest.raises(complexity.ComplexityFilterError, match="at least one eligible"):
        complexity.scan_manifest(path)


@pytest.mark.parametrize(("content", "error"), [("ACGT\n", "before the first header"), (">empty\n", "empty record"), (">bad\nACGN\n", "non-ACGT"), (">short\nACG\n", "wrong-length"), (">a\nACGT\n>b\nACGT\n", "duplicate canonical")])
def test_deacon_passing_fasta_requires_well_formed_unique_eligible_kmers(tmp_path: Path, content: str, error: str) -> None:
    path = tmp_path / "passing.fasta"
    path.write_text(content)
    with pytest.raises(complexity.ComplexityFilterError, match=error):
        complexity.read_passing_kmers(path, 4)


def test_complexity_filter_matches_a_deacon_reverse_complement_to_its_candidate_kmer(tmp_path: Path) -> None:
    manifest_in = tmp_path / "candidate_kmers.tsv"
    rows = valid_manifest_rows()
    rows[0]["kmer"] = KMER_B
    pl.DataFrame(rows[:1], schema=complexity.MANIFEST_SCHEMA).write_csv(manifest_in, separator="\t")
    passing_kmers = tmp_path / "passing_kmers.fasta"
    passing_kmers.write_text(f">deacon_record\n{KMER_B_RC}\n")
    manifest_out = tmp_path / "filtered_candidate_kmers.tsv"
    baits_out = tmp_path / "locally_filtered_baits.fasta"

    complexity.main([
        "--design-id", "design",
        "--source-sequence-origin", "curated_input",
        "--manifest-in", str(manifest_in),
        "--passing-kmers", str(passing_kmers),
        "--manifest-out", str(manifest_out),
        "--baits-out", str(baits_out),
        "--filtering-status-out", str(tmp_path / "filtering_status.tsv"),
        "--bait-set-status-out", str(tmp_path / "bait_set_status.tsv"),
    ])

    assert pl.read_csv(manifest_out, separator="\t").select("kmer", "status").rows() == [(KMER_B, "PASS")]
    assert baits_out.read_text() == f">bait_000001\n{KMER_B}\n"


def test_complexity_filter_rejects_unknown_pass_key_and_never_reactivates_background_rejection(tmp_path: Path) -> None:
    path = tmp_path / "manifest.tsv"
    pl.DataFrame(valid_manifest_rows(), schema=complexity.MANIFEST_SCHEMA).write_csv(path, separator="\t")
    base = complexity.scan_manifest(path)
    passing = tmp_path / "passing.fasta"
    passing.write_text(f">unknown\n{KMER_B}\n")
    with pytest.raises(complexity.ComplexityFilterError, match="not eligible"):
        complexity.main(["--design-id", "design", "--source-sequence-origin", "curated_input", "--manifest-in", str(path), "--passing-kmers", str(passing), "--manifest-out", str(tmp_path / "out.tsv"), "--baits-out", str(tmp_path / "baits.fasta"), "--filtering-status-out", str(tmp_path / "status.tsv"), "--bait-set-status-out", str(tmp_path / "bait-status.tsv")])
    result = complexity.apply_complexity_results(base, pl.LazyFrame({"kmer": []}, schema={"kmer": pl.String})).collect()
    assert result.select("status", "rejection_reason", "taxonomic_screening_status").rows() == [
        ("REJECT_LOW_COMPLEXITY", "low_complexity", "NOT_APPLICABLE"),
        ("REJECT_INTERFERENCE_BACKGROUND", "background_occurrence", "NOT_APPLICABLE"),
    ]


def test_complexity_filter_assigns_bait_ids_deterministically_by_kmer(tmp_path: Path) -> None:
    rows = valid_manifest_rows()
    rows[1] = dict(zip(complexity.MANIFEST_FIELDS, ("candidate_kmer_000002", "", KMER_B, "1", "0", "PASS", "none", "NOT_RUN", "", ""), strict=True))
    path = tmp_path / "manifest.tsv"
    pl.DataFrame(rows, schema=complexity.MANIFEST_SCHEMA).write_csv(path, separator="\t")
    result = complexity.apply_complexity_results(complexity.scan_manifest(path), pl.LazyFrame({"kmer": [KMER_A, KMER_B]})).collect()
    assert result.sort("kmer").select("kmer", "bait_id").rows() == [(KMER_A, "bait_000001"), (KMER_B, "bait_000002")]


def test_build_command_writes_terminal_outputs_for_no_candidates_and_all_background(tmp_path: Path) -> None:
    for name, sequence, source_counts, background, terminal in (("none", "ACGT", "", "", "candidate_kmer_enumeration"), ("background", KMER_A, f"{KMER_A}\t1\n", f"{KMER_A}\t3\n", "explicit_background_cancellation")):
        paths = write_build_inputs(tmp_path / name, sequence, source_counts=source_counts, background=background)
        out = {key: tmp_path / name / f"{key}.{'fasta' if key == 'fasta' else 'tsv'}" for key in ("manifest", "occurrences", "fasta", "status", "terminal")}
        candidates.main(["--source-sequences", str(paths["source"]), "--source-sequence-query-groups", str(paths["groups"]), "--meryl-source-counts", str(paths["source_counts"]), "--background-intersection-counts", str(paths["background"]), "--kmer-size", str(len(KMER_A)), "--design-id", "design", "--manifest-out", str(out["manifest"]), "--occurrences-out", str(out["occurrences"]), "--complexity-candidates-out", str(out["fasta"]), "--filtering-status-out", str(out["status"]), "--terminal-manifest-out", str(out["terminal"])])
        assert not out["fasta"].exists()
        assert pl.read_csv(out["status"], separator="\t").filter(pl.col("metric") == "terminal_stage").item(0, "value") == terminal
        assert pl.read_csv(out["terminal"], separator="\t").height == (0 if name == "none" else 1)


def test_complexity_command_writes_terminal_status_and_no_fasta_when_no_baits(tmp_path: Path) -> None:
    input_path = tmp_path / "manifest.tsv"
    pl.DataFrame(valid_manifest_rows()[:1], schema=complexity.MANIFEST_SCHEMA).write_csv(input_path, separator="\t")
    passing = tmp_path / "passing.fasta"
    passing.write_text("")
    output, baits, status, terminal_manifest = (
        tmp_path / name
        for name in ("out.tsv", "baits.fasta", "status.tsv", "terminal_candidate_kmers.tsv")
    )
    complexity.main(["--design-id", "design", "--source-sequence-origin", "curated_input", "--manifest-in", str(input_path), "--passing-kmers", str(passing), "--manifest-out", str(output), "--baits-out", str(baits), "--filtering-status-out", str(status), "--bait-set-status-out", str(tmp_path / "bait-status.tsv"), "--terminal-manifest-out", str(terminal_manifest)])
    assert not baits.exists()
    assert pl.read_csv(status, separator="\t").filter(pl.col("metric") == "terminal_stage").item(0, "value") == "low_complexity_filtering"
    assert pl.read_csv(output, separator="\t").item(0, "status") == "REJECT_LOW_COMPLEXITY"
    assert terminal_manifest.read_text() == output.read_text()
