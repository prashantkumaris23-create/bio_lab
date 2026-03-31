from __future__ import annotations

from typing import List, Tuple

from src.phylo import TreeNode


def render_heatmap_html(matrix: List[List[float]], labels: List[str] | None = None) -> str:
    labels = labels or [f"S{i + 1}" for i in range(len(matrix))]
    header_cells = "".join(_header_cell(label) for label in [""] + labels)
    rows = [f"<tr>{header_cells}</tr>"]
    for row_label, row in zip(labels, matrix):
        row_html = [_header_cell(row_label)]
        for value in row:
            row_html.append(
                f'<td style="padding:0.55rem 0.65rem; text-align:center;'
                f' background:{_color_for_value(value)}; border:1px solid rgba(125,163,201,0.10);'
                f' color:#eff7ff; font-weight:600;">{value:.3f}</td>'
            )
        rows.append("<tr>" + "".join(row_html) + "</tr>")
    return (
        '<div style="overflow-x:auto; background:linear-gradient(180deg, rgba(15,24,43,0.96) 0%, rgba(10,18,32,0.98) 100%); border-radius:16px;'
        ' box-shadow:0 12px 30px rgba(0,0,0,0.22); border:1px solid rgba(125,163,201,0.12);">'
        '<table style="border-collapse:collapse; width:100%; min-width:380px;">'
        + "".join(rows)
        + "</table></div>"
    )


def _header_cell(label: str) -> str:
    return (
        '<th style="padding:0.6rem 0.7rem; background:#10233d; color:#f3fbff;'
        f' border:1px solid rgba(255,255,255,0.06);">{label}</th>'
    )


def _color_for_value(value: float) -> str:
    value = max(0.0, min(1.0, value))
    red = int(22 + 20 * value)
    green = int(33 + 155 * value)
    blue = int(48 + 110 * value)
    return f"rgb({red},{green},{blue})"


def render_tree_svg(tree: TreeNode, width: int = 920, height: int = 420) -> str:
    leaf_order = _collect_leaves(tree)
    if not leaf_order:
        return "<div>No tree available</div>"

    max_distance = max(_leaf_distances(tree)) or 1.0
    left_pad = 54
    right_pad = 260
    top_pad = 44
    bottom_pad = 40
    usable_width = width - left_pad - right_pad
    usable_height = height - top_pad - bottom_pad
    y_step = usable_height / max(1, len(leaf_order) - 1)

    leaf_y = {
        leaf.name: top_pad + index * y_step
        for index, leaf in enumerate(leaf_order)
    }

    lines: List[str] = []
    labels: List[str] = []
    branch_lengths: List[str] = []
    ornaments: List[str] = []

    def visit(node: TreeNode, x_position: float) -> Tuple[float, float]:
        if node.is_leaf:
            y_position = leaf_y[node.name]
            tip_x = left_pad + usable_width + 18
            label_x = tip_x + 18
            safe_name = _xml_escape(node.name)
            label_width = max(110, min(210, 10 * len(node.name) + 28))
            label_rect_y = y_position - 15
            labels.append(
                f'<rect x="{label_x - 10:.1f}" y="{label_rect_y:.1f}" width="{label_width:.1f}" height="28"'
                ' rx="14" ry="14" fill="#0f1a2d" stroke="rgba(115,240,214,0.16)" />'
            )
            labels.append(
                f'<text x="{label_x:.1f}" y="{y_position + 5:.1f}" font-size="13.5"'
                f' font-family="Segoe UI, sans-serif" font-weight="600" fill="#eef7ff">{safe_name}</text>'
            )
            lines.append(
                f'<line x1="{x_position:.1f}" y1="{y_position:.1f}" '
                f'x2="{tip_x:.1f}" y2="{y_position:.1f}" '
                'stroke="#79a2c1" stroke-width="2.1" stroke-linecap="round" />'
            )
            ornaments.append(
                f'<circle cx="{tip_x:.1f}" cy="{y_position:.1f}" r="4.6" fill="#ff8a3d" stroke="#091121" stroke-width="2" />'
            )
            return x_position, y_position

        left_x = x_position + (node.left_length / max_distance) * usable_width
        right_x = x_position + (node.right_length / max_distance) * usable_width
        left_child_x, left_y_pos = visit(node.left, left_x)
        right_child_x, right_y_pos = visit(node.right, right_x)

        lines.append(
            f'<line x1="{x_position:.1f}" y1="{left_y_pos:.1f}" '
            f'x2="{x_position:.1f}" y2="{right_y_pos:.1f}" '
            'stroke="#33d2c0" stroke-width="3" stroke-linecap="round" />'
        )
        lines.append(
            f'<line x1="{x_position:.1f}" y1="{left_y_pos:.1f}" '
            f'x2="{left_child_x:.1f}" y2="{left_y_pos:.1f}" '
            'stroke="#33d2c0" stroke-width="3" stroke-linecap="round" />'
        )
        lines.append(
            f'<line x1="{x_position:.1f}" y1="{right_y_pos:.1f}" '
            f'x2="{right_child_x:.1f}" y2="{right_y_pos:.1f}" '
            'stroke="#33d2c0" stroke-width="3" stroke-linecap="round" />'
        )
        node_y = (left_y_pos + right_y_pos) / 2
        ornaments.append(
            f'<circle cx="{x_position:.1f}" cy="{node_y:.1f}" r="4.2" fill="#7ff3e0" stroke="#091121" stroke-width="1.8" />'
        )
        if node.left_length > 0:
            branch_lengths.append(
                f'<text x="{(x_position + left_child_x) / 2:.1f}" y="{left_y_pos - 10:.1f}" font-size="11.5"'
                f' font-family="Segoe UI, sans-serif" fill="#91a9c0">{node.left_length:.3f}</text>'
            )
        if node.right_length > 0:
            branch_lengths.append(
                f'<text x="{(x_position + right_child_x) / 2:.1f}" y="{right_y_pos - 10:.1f}" font-size="11.5"'
                f' font-family="Segoe UI, sans-serif" fill="#91a9c0">{node.right_length:.3f}</text>'
            )
        return x_position, node_y

    visit(tree, left_pad)
    scale_x1 = left_pad
    scale_x2 = left_pad + min(usable_width * 0.18, 120)
    scale_y = top_pad - 18
    background = (
        f'<rect x="0" y="0" width="{width}" height="{height}" rx="20" ry="20" '
        'fill="url(#panelGradient)" stroke="rgba(125,163,201,0.10)" />'
    )
    return (
        f'<div style="overflow-x:auto;"><svg width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
        '<defs>'
        '<linearGradient id="panelGradient" x1="0" y1="0" x2="1" y2="1">'
        '<stop offset="0%" stop-color="#0d1729" />'
        '<stop offset="100%" stop-color="#0a1220" />'
        '</linearGradient>'
        '</defs>'
        + background
        + f'<text x="{left_pad:.1f}" y="{top_pad - 24:.1f}" font-size="12" font-family="Segoe UI, sans-serif" fill="#91a9c0">Branch length scale</text>'
        + f'<line x1="{scale_x1:.1f}" y1="{scale_y:.1f}" x2="{scale_x2:.1f}" y2="{scale_y:.1f}" stroke="#ff8a3d" stroke-width="3" stroke-linecap="round" />'
        + f'<text x="{scale_x2 + 8:.1f}" y="{scale_y + 4:.1f}" font-size="11.5" font-family="Segoe UI, sans-serif" fill="#91a9c0">{max_distance * (scale_x2 - scale_x1) / usable_width:.3f}</text>'
        + "".join(lines)
        + "".join(branch_lengths)
        + "".join(ornaments)
        + "".join(labels)
        + "</svg></div>"
    )


def _collect_leaves(node: TreeNode) -> List[TreeNode]:
    if node.is_leaf:
        return [node]
    leaves: List[TreeNode] = []
    leaves.extend(_collect_leaves(node.left))
    leaves.extend(_collect_leaves(node.right))
    return leaves


def _leaf_distances(node: TreeNode, current_distance: float = 0.0) -> List[float]:
    if node.is_leaf:
        return [current_distance]
    distances: List[float] = []
    distances.extend(_leaf_distances(node.left, current_distance + node.left_length))
    distances.extend(_leaf_distances(node.right, current_distance + node.right_length))
    return distances


def _xml_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )
