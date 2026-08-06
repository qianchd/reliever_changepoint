# GitHub usage metrics

This folder stores the updater used by
`.github/workflows/update-usage-metrics.yml`. Collected records live under
`data/`, separately from the script and its tests:

- `data/github_traffic_daily.csv` accumulates the delayed daily view, visitor,
  clone, and cloner series returned by the GitHub traffic API.
- `data/github_usage_snapshots.csv` records one repository-level summary per
  collection date, including stars, forks, release-asset downloads, and the
  cumulative traffic values available on that date.

These files record user-facing access metrics only:

- GitHub stars
- GitHub forks
- GitHub Release asset downloads
- repository views and unique visitors over GitHub's rolling 14-day window
- repository clones and unique cloners over GitHub's rolling 14-day window

The workflow does not write these metrics into the package `README.md`, and it
does not collect development-maintenance metrics such as issues, pull requests,
CI runs, commits, or contributors.

Traffic endpoints require authenticated repository traffic access. Configure a
repository secret named `GH_TRAFFIC_TOKEN`; a fine-grained token scoped to this
repository with `Administration: Read` is sufficient for the GitHub traffic
API. The built-in `GITHUB_TOKEN` is still used by checkout for committing the
data files.
