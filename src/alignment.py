from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass(frozen=True)
class SequenceRecord:
    identifier: str
    label: str
    description: str
    sequence: str


def sanitize_sequence(sequence: str) -> str:
    allowed = {"A", "C", "G", "T", "N", "-"}
    cleaned = "".join(char for char in sequence.upper() if not char.isspace())
    invalid = {char for char in cleaned if char not in allowed}
    if invalid:
        raise ValueError(f"Sequence contains unsupported symbols: {', '.join(sorted(invalid))}")
    return cleaned


def parse_fasta(text: str) -> List[SequenceRecord]:
    """Parse FASTA text into strongly typed sequence records."""
    records: List[SequenceRecord] = []
    current_header = ""
    current_lines: List[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                records.append(_record_from_header(current_header, "".join(current_lines)))
            current_header = line[1:]
            current_lines = []
        else:
            current_lines.append(line)

    if current_header:
        records.append(_record_from_header(current_header, "".join(current_lines)))

    return records


def _record_from_header(header: str, sequence: str) -> SequenceRecord:
    parts = header.split(maxsplit=1)
    identifier = parts[0].strip()
    description = parts[1].strip() if len(parts) > 1 else identifier
    label = identifier
    return SequenceRecord(
        identifier=identifier,
        label=label,
        description=description,
        sequence=sanitize_sequence(sequence),
    )


def format_fasta(records: List[SequenceRecord], line_length: int = 80) -> str:
    lines: List[str] = []
    for record in records:
        lines.append(f">{record.identifier} {record.description}".strip())
        for index in range(0, len(record.sequence), line_length):
            lines.append(record.sequence[index : index + line_length])
    return "\n".join(lines)


def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -1,
) -> Tuple[str, str]:
    rows = len(seq_a) + 1
    cols = len(seq_b) + 1
    scores = [[0 for _ in range(cols)] for _ in range(rows)]
    trace = [["" for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        scores[i][0] = i * gap_penalty
        trace[i][0] = "up"
    for j in range(1, cols):
        scores[0][j] = j * gap_penalty
        trace[0][j] = "left"

    for i in range(1, rows):
        for j in range(1, cols):
            diag = scores[i - 1][j - 1] + (
                match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_score
            )
            up = scores[i - 1][j] + gap_penalty
            left = scores[i][j - 1] + gap_penalty
            best = max(diag, up, left)
            scores[i][j] = best
            if best == diag:
                trace[i][j] = "diag"
            elif best == up:
                trace[i][j] = "up"
            else:
                trace[i][j] = "left"

    aligned_a: List[str] = []
    aligned_b: List[str] = []
    i = len(seq_a)
    j = len(seq_b)

    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == "diag":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def build_multiple_alignment(records: List[SequenceRecord]) -> List[SequenceRecord]:
    if len(records) < 2:
        raise ValueError("At least 2 sequences are required for alignment.")

    # Use the longest sequence as the temporary alignment anchor.
    center_index = max(range(len(records)), key=lambda idx: len(records[idx].sequence))
    center_record = records[center_index]
    aligned_records = {
        center_record.identifier: center_record.sequence,
    }
    current_center = center_record.sequence

    for idx, record in enumerate(records):
        if idx == center_index:
            continue
        pair_center, pair_other = needleman_wunsch(center_record.sequence, record.sequence)
        current_center, aligned_records = _merge_alignment(
            current_center=current_center,
            existing_alignment=aligned_records,
            pair_center=pair_center,
            pair_other=pair_other,
            pair_identifier=record.identifier,
        )

    result: List[SequenceRecord] = []
    for record in records:
        result.append(
            SequenceRecord(
                identifier=record.identifier,
                label=record.label,
                description=record.description,
                sequence=aligned_records[record.identifier],
            )
        )
    return result


def _merge_alignment(
    current_center: str,
    existing_alignment: dict[str, str],
    pair_center: str,
    pair_other: str,
    pair_identifier: str,
) -> Tuple[str, dict[str, str]]:
    # Merge one new pairwise alignment into the current multiple alignment while
    # preserving gap positions that were already introduced earlier.
    merged_existing = {key: [] for key in existing_alignment}
    merged_pair: List[str] = []
    merged_center: List[str] = []

    i = 0
    j = 0
    while i < len(current_center) or j < len(pair_center):
        current_char = current_center[i] if i < len(current_center) else None
        pair_char = pair_center[j] if j < len(pair_center) else None

        if current_char == pair_char and current_char is not None:
            merged_center.append(current_char)
            for key, sequence in existing_alignment.items():
                merged_existing[key].append(sequence[i])
            merged_pair.append(pair_other[j])
            i += 1
            j += 1
            continue

        if current_char == "-":
            merged_center.append(current_char)
            for key, sequence in existing_alignment.items():
                merged_existing[key].append(sequence[i])
            merged_pair.append("-")
            i += 1
            continue

        if pair_char == "-":
            merged_center.append(pair_char)
            for key in merged_existing:
                merged_existing[key].append("-")
            merged_pair.append(pair_other[j])
            j += 1
            continue

        if current_char is None and pair_char is not None:
            merged_center.append(pair_char)
            for key in merged_existing:
                merged_existing[key].append("-")
            merged_pair.append(pair_other[j])
            j += 1
            continue

        if pair_char is None and current_char is not None:
            merged_center.append(current_char)
            for key, sequence in existing_alignment.items():
                merged_existing[key].append(sequence[i])
            merged_pair.append("-")
            i += 1
            continue

        raise ValueError("Unable to merge alignments because center sequences diverged unexpectedly.")

    updated_alignment = {key: "".join(chars) for key, chars in merged_existing.items()}
    updated_alignment[pair_identifier] = "".join(merged_pair)
    return "".join(merged_center), updated_alignment
