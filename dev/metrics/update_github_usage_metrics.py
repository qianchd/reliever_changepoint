#!/usr/bin/env python3
"""Accumulate GitHub repository traffic and update the README summary."""

from __future__ import annotations

import csv
import json
import os
import sys
import urllib.error
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path


API_ROOT = "https://api.github.com"
TRAFFIC_FIELDS = [
    "date_utc",
    "repository",
    "views",
    "unique_visitors",
    "clones",
    "unique_cloners",
]
USAGE_SNAPSHOT_FIELDS = [
    "snapshot_date_utc",
    "repository",
    "traffic_start_utc",
    "traffic_end_utc",
    "stars",
    "forks",
    "release_asset_downloads",
    "views",
    "daily_unique_visitors",
    "clones",
    "daily_unique_cloners",
]


def repository_name() -> str:
    return os.environ.get("GITHUB_REPOSITORY", "qianchd/reliever_changepoint")


def github_token() -> str | None:
    return os.environ.get("GH_TRAFFIC_TOKEN") or os.environ.get("GITHUB_TOKEN")


def api_get(path: str, token: str) -> object | None:
    request = urllib.request.Request(
        API_ROOT + path,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "reliever-usage-metrics",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace").strip()
        permission = exc.headers.get("X-Accepted-GitHub-Permissions", "")
        message = f"GitHub API request failed for {path}: HTTP {exc.code}"
        if permission:
            message += f"; accepted permissions: {permission}"
        if detail:
            message += f"; response: {detail}"
        print(f"error: {message}", file=sys.stderr)
    except urllib.error.URLError as exc:
        print(f"error: GitHub API request failed for {path}: {exc}", file=sys.stderr)
    return None


def release_downloads(repo: str, token: str) -> int | None:
    total = 0
    page = 1
    while True:
        releases = api_get(f"/repos/{repo}/releases?per_page=100&page={page}", token)
        if releases is None:
            return None
        if not isinstance(releases, list) or not releases:
            return total
        for release in releases:
            for asset in release.get("assets", []):
                total += int(asset.get("download_count", 0))
        if len(releases) < 100:
            return total
        page += 1


def parse_daily_traffic(payload: object, series_name: str) -> dict[str, tuple[int, int]]:
    if not isinstance(payload, dict) or not isinstance(payload.get(series_name), list):
        raise ValueError(f"GitHub traffic response has no '{series_name}' daily series")

    daily: dict[str, tuple[int, int]] = {}
    for item in payload[series_name]:
        timestamp = str(item.get("timestamp", ""))
        day = timestamp[:10]
        try:
            date.fromisoformat(day)
        except ValueError as exc:
            raise ValueError(f"invalid traffic timestamp: {timestamp!r}") from exc
        daily[day] = (int(item.get("count", 0)), int(item.get("uniques", 0)))
    return daily


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_daily_traffic(path: Path) -> list[dict[str, str]]:
    return load_csv_rows(path)


def load_usage_snapshots(path: Path) -> list[dict[str, str]]:
    return load_csv_rows(path)


def merge_daily_traffic(
    rows: list[dict[str, str]],
    repo: str,
    views: dict[str, tuple[int, int]],
    clones: dict[str, tuple[int, int]],
    cutoff: date,
) -> list[dict[str, str]]:
    by_key = {
        (row["repository"], row["date_utc"]): dict(row)
        for row in rows
    }
    for day in sorted(set(views) | set(clones)):
        if date.fromisoformat(day) > cutoff:
            continue
        previous = by_key.get((repo, day), {})
        view_count, visitor_count = views.get(
            day,
            (
                int(previous.get("views", 0)),
                int(previous.get("unique_visitors", 0)),
            ),
        )
        clone_count, cloner_count = clones.get(
            day,
            (
                int(previous.get("clones", 0)),
                int(previous.get("unique_cloners", 0)),
            ),
        )
        by_key[(repo, day)] = {
            "date_utc": day,
            "repository": repo,
            "views": str(view_count),
            "unique_visitors": str(visitor_count),
            "clones": str(clone_count),
            "unique_cloners": str(cloner_count),
        }
    return [by_key[key] for key in sorted(by_key)]


def write_csv_rows(
    path: Path, rows: list[dict[str, str]], fieldnames: list[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_daily_traffic(path: Path, rows: list[dict[str, str]]) -> None:
    write_csv_rows(path, rows, TRAFFIC_FIELDS)


def write_usage_snapshots(path: Path, rows: list[dict[str, str]]) -> None:
    write_csv_rows(path, rows, USAGE_SNAPSHOT_FIELDS)


def summarize_traffic(
    rows: list[dict[str, str]], repo: str, cutoff: date
) -> dict[str, int | str | None]:
    selected = [
        row for row in rows
        if row["repository"] == repo and date.fromisoformat(row["date_utc"]) <= cutoff
    ]
    if not selected:
        return {
            "traffic_start": None,
            "traffic_end": None,
            "views": None,
            "unique_visitors": None,
            "clones": None,
            "unique_cloners": None,
        }
    return {
        "traffic_start": min(row["date_utc"] for row in selected),
        "traffic_end": max(row["date_utc"] for row in selected),
        "views": sum(int(row["views"]) for row in selected),
        "unique_visitors": sum(int(row["unique_visitors"]) for row in selected),
        "clones": sum(int(row["clones"]) for row in selected),
        "unique_cloners": sum(int(row["unique_cloners"]) for row in selected),
    }


def merge_usage_snapshot(
    rows: list[dict[str, str]],
    metrics: dict[str, int | str | None],
) -> list[dict[str, str]]:
    snapshot = {
        field: "" if metrics.get(field) is None else str(metrics.get(field))
        for field in USAGE_SNAPSHOT_FIELDS
    }
    snapshot_date = snapshot["snapshot_date_utc"]
    repository = snapshot["repository"]
    date.fromisoformat(snapshot_date)
    if not repository:
        raise ValueError("usage snapshot repository must not be empty")

    by_key = {
        (row["repository"], row["snapshot_date_utc"]): dict(row)
        for row in rows
    }
    by_key[(repository, snapshot_date)] = snapshot
    return [by_key[key] for key in sorted(by_key)]


def main() -> int:
    repo = repository_name()
    token = github_token()
    if not token:
        print(
            "error: GH_TRAFFIC_TOKEN is required. Add a fine-grained GitHub token "
            "with repository Administration read permission.",
            file=sys.stderr,
        )
        return 2

    repo_info = api_get(f"/repos/{repo}", token)
    views_payload = api_get(f"/repos/{repo}/traffic/views?per=day", token)
    clones_payload = api_get(f"/repos/{repo}/traffic/clones?per=day", token)
    downloads = release_downloads(repo, token)
    if not isinstance(repo_info, dict) or views_payload is None or clones_payload is None:
        print("error: metrics were not updated because required API data is unavailable", file=sys.stderr)
        return 2

    try:
        views = parse_daily_traffic(views_payload, "views")
        clones = parse_daily_traffic(clones_payload, "clones")
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    today = datetime.now(timezone.utc).date()
    cutoff = today - timedelta(days=3)
    root = Path(__file__).resolve().parents[2]
    data_dir = root / "dev" / "metrics" / "data"
    traffic_path = data_dir / "github_traffic_daily.csv"
    snapshot_path = data_dir / "github_usage_snapshots.csv"
    traffic_rows = merge_daily_traffic(
        load_daily_traffic(traffic_path), repo, views, clones, cutoff
    )
    traffic_summary = summarize_traffic(traffic_rows, repo, cutoff)
    metrics: dict[str, int | str | None] = {
        "snapshot_date_utc": today.isoformat(),
        "repository": repo,
        "stars": int(repo_info.get("stargazers_count", 0)),
        "forks": int(repo_info.get("forks_count", 0)),
        "release_asset_downloads": downloads,
        "traffic_start_utc": traffic_summary["traffic_start"],
        "traffic_end_utc": traffic_summary["traffic_end"],
        "views": traffic_summary["views"],
        "daily_unique_visitors": traffic_summary["unique_visitors"],
        "clones": traffic_summary["clones"],
        "daily_unique_cloners": traffic_summary["unique_cloners"],
    }

    write_daily_traffic(traffic_path, traffic_rows)
    snapshot_rows = merge_usage_snapshot(
        load_usage_snapshots(snapshot_path), metrics
    )
    write_usage_snapshots(snapshot_path, snapshot_rows)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
