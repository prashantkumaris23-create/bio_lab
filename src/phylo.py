from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from src.alignment import SequenceRecord


@dataclass
class TreeNode:
    name: str
    left: Optional["TreeNode"] = None
    right: Optional["TreeNode"] = None
    left_length: float = 0.0
    right_length: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return self.left is None and self.right is None


def similarity_matrix(records: List[SequenceRecord]) -> List[List[float]]:
    matrix: List[List[float]] = []
    for record_a in records:
        row: List[float] = []
        for record_b in records:
            row.append(pairwise_identity(record_a.sequence, record_b.sequence))
        matrix.append(row)
    return matrix


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    # Ignore columns that are gaps in both sequences so the identity score only
    # reflects informative aligned positions.
    matches = 0
    compared = 0
    for char_a, char_b in zip(seq_a, seq_b):
        if char_a == "-" and char_b == "-":
            continue
        compared += 1
        if char_a == char_b and char_a != "-":
            matches += 1
    return matches / compared if compared else 1.0


def distance_matrix_from_similarity(similarity: List[List[float]]) -> List[List[float]]:
    return [[max(0.0, 1.0 - value) for value in row] for row in similarity]


def build_neighbor_joining_tree(
    records: List[SequenceRecord],
    distance_matrix: List[List[float]],
) -> TreeNode:
    labels = [record.label for record in records]
    nodes = {label: TreeNode(name=label) for label in labels}
    distances = _distance_dict(labels, distance_matrix)
    active = list(labels)
    next_internal_id = 1

    while len(active) > 2:
        size = len(active)
        row_sums = {
            label: sum(distances[label][other] for other in active if other != label)
            for label in active
        }

        best_pair = None
        best_q = None
        for i in range(len(active)):
            for j in range(i + 1, len(active)):
                a = active[i]
                b = active[j]
                # Neighbor Joining uses the Q-matrix to decide which pair should
                # be merged next.
                q_value = (size - 2) * distances[a][b] - row_sums[a] - row_sums[b]
                if best_q is None or q_value < best_q:
                    best_q = q_value
                    best_pair = (a, b)

        if best_pair is None:
            raise ValueError("Neighbor Joining failed to select a node pair.")

        a, b = best_pair
        distance_ab = distances[a][b]
        limb_a = 0.5 * distance_ab + (row_sums[a] - row_sums[b]) / (2 * (size - 2))
        limb_b = distance_ab - limb_a
        limb_a = max(limb_a, 0.0)
        limb_b = max(limb_b, 0.0)

        parent_name = f"Inner{next_internal_id}"
        next_internal_id += 1
        nodes[parent_name] = TreeNode(
            name=parent_name,
            left=nodes[a],
            right=nodes[b],
            left_length=limb_a,
            right_length=limb_b,
        )

        distances[parent_name] = {}
        for other in active:
            if other in {a, b}:
                continue
            new_distance = 0.5 * (
                distances[a][other] + distances[b][other] - distance_ab
            )
            distances[parent_name][other] = new_distance
            distances[other][parent_name] = new_distance

        for label in [a, b]:
            distances.pop(label, None)
        for row in distances.values():
            row.pop(a, None)
            row.pop(b, None)

        active = [label for label in active if label not in {a, b}]
        active.append(parent_name)

    left_label, right_label = active
    final_distance = distances[left_label][right_label]
    return TreeNode(
        name="Root",
        left=nodes[left_label],
        right=nodes[right_label],
        left_length=max(final_distance / 2, 0.0),
        right_length=max(final_distance / 2, 0.0),
    )


def _distance_dict(labels: List[str], matrix: List[List[float]]) -> Dict[str, Dict[str, float]]:
    distances: Dict[str, Dict[str, float]] = {label: {} for label in labels}
    for i, label_a in enumerate(labels):
        for j, label_b in enumerate(labels):
            if i == j:
                continue
            distances[label_a][label_b] = matrix[i][j]
    return distances


def newick_from_tree(node: TreeNode) -> str:
    if node.is_leaf:
        return _sanitize_newick_label(node.name)

    left_text = newick_from_tree(node.left) + f":{node.left_length:.5f}"
    right_text = newick_from_tree(node.right) + f":{node.right_length:.5f}"
    return f"({left_text},{right_text})"


def matrix_to_tsv(matrix: List[List[float]], labels: Optional[List[str]] = None) -> str:
    resolved_labels = labels or [f"S{i + 1}" for i in range(len(matrix))]
    header = ["label"] + resolved_labels
    lines = ["\t".join(header)]
    for index, row in enumerate(matrix):
        lines.append(
            "\t".join([resolved_labels[index]] + [f"{value:.6f}" for value in row])
        )
    return "\n".join(lines)


def average_similarity(matrix: List[List[float]]) -> float:
    values: List[float] = []
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if j > i:
                values.append(value)
    return sum(values) / len(values) if values else 1.0


def _sanitize_newick_label(label: str) -> str:
    safe = label.replace(" ", "_")
    for char in "(),:;":
        safe = safe.replace(char, "_")
    return safe
