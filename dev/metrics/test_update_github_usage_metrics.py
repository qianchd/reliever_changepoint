import tempfile
import unittest
from datetime import date
from pathlib import Path

from update_github_usage_metrics import (
    load_daily_traffic,
    load_usage_snapshots,
    merge_daily_traffic,
    merge_usage_snapshot,
    parse_daily_traffic,
    summarize_traffic,
    write_daily_traffic,
    write_usage_snapshots,
)


class UsageMetricsTest(unittest.TestCase):
    def test_parse_daily_traffic(self):
        payload = {
            "views": [
                {"timestamp": "2026-07-12T00:00:00Z", "count": 8, "uniques": 3}
            ]
        }
        self.assertEqual(parse_daily_traffic(payload, "views"), {"2026-07-12": (8, 3)})

    def test_merge_is_cutoff_aware_and_idempotent(self):
        existing = [{
            "date_utc": "2026-07-10",
            "repository": "qianchd/reliever_changepoint",
            "views": "2",
            "unique_visitors": "1",
            "clones": "0",
            "unique_cloners": "0",
        }]
        views = {
            "2026-07-10": (5, 2),
            "2026-07-12": (8, 3),
            "2026-07-13": (99, 9),
        }
        clones = {
            "2026-07-10": (1, 1),
            "2026-07-12": (2, 1),
            "2026-07-13": (99, 9),
        }
        cutoff = date(2026, 7, 12)
        merged = merge_daily_traffic(
            existing, "qianchd/reliever_changepoint", views, clones, cutoff
        )
        merged_again = merge_daily_traffic(
            merged, "qianchd/reliever_changepoint", views, clones, cutoff
        )

        self.assertEqual(merged, merged_again)
        self.assertEqual([row["date_utc"] for row in merged], ["2026-07-10", "2026-07-12"])
        self.assertEqual(merged[0]["views"], "5")
        summary = summarize_traffic(
            merged, "qianchd/reliever_changepoint", cutoff
        )
        self.assertEqual(summary["views"], 13)
        self.assertEqual(summary["clones"], 3)

        views_only = merge_daily_traffic(
            merged,
            "qianchd/reliever_changepoint",
            {"2026-07-12": (9, 4)},
            {},
            cutoff,
        )
        self.assertEqual(views_only[1]["views"], "9")
        self.assertEqual(views_only[1]["clones"], "2")

    def test_daily_csv_round_trip(self):
        rows = [{
            "date_utc": "2026-07-12",
            "repository": "qianchd/reliever_changepoint",
            "views": "8",
            "unique_visitors": "3",
            "clones": "2",
            "unique_cloners": "1",
        }]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "traffic.csv"
            write_daily_traffic(path, rows)
            self.assertEqual(load_daily_traffic(path), rows)

    def test_snapshot_csv_is_idempotent_per_repository_and_date(self):
        first = {
            "snapshot_date_utc": "2026-07-15",
            "repository": "qianchd/reliever_changepoint",
            "traffic_start_utc": "2026-07-01",
            "traffic_end_utc": "2026-07-12",
            "stars": "1",
            "forks": "2",
            "release_asset_downloads": "3",
            "views": "10",
            "daily_unique_visitors": "5",
            "clones": "4",
            "daily_unique_cloners": "2",
        }
        replacement = dict(first, stars="4", views="12")
        later = dict(
            replacement,
            snapshot_date_utc="2026-07-23",
            traffic_end_utc="2026-07-20",
            views="20",
        )

        rows = merge_usage_snapshot([], first)
        rows = merge_usage_snapshot(rows, replacement)
        rows = merge_usage_snapshot(rows, later)

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["stars"], "4")
        self.assertEqual(rows[0]["views"], "12")
        self.assertEqual(rows[1]["snapshot_date_utc"], "2026-07-23")

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshots.csv"
            write_usage_snapshots(path, rows)
            self.assertEqual(load_usage_snapshots(path), rows)


if __name__ == "__main__":
    unittest.main()
