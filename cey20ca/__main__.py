"""Command line: python -m cey20ca {download,reference,metrics}."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import compute, download, reference


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m cey20ca", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("download", help="download the Zenodo files and verify MD5 and SHA-256")
    d.add_argument("--data", type=Path, default=Path("data"))
    r = sub.add_parser("reference", help="download ERA5-Land from the CDS and pack the 20 reference years")
    r.add_argument("--data", type=Path, default=Path("data"))
    r.add_argument("--cache", type=Path, default=Path("reference/packed"))
    r.add_argument("--downloads", type=Path, default=Path("reference/cds"))
    r.add_argument("--workers", type=int, default=4, help="parallel CDS requests")
    m = sub.add_parser("metrics", help="recompute point_metrics.csv and compare with the published file")
    m.add_argument("--data", type=Path, default=Path("data"))
    m.add_argument("--cache", type=Path, default=Path("reference/packed"))
    m.add_argument("--output", type=Path, default=Path("point_metrics_recomputed.csv"))
    args = parser.parse_args(argv)

    if args.command == "download":
        download.fetch(args.data)
        print(f"verified {len(download.verify(args.data))} files against manifest.json")
        return 0
    nc_path = args.data / "cey20_v2_california_points.nc"
    if args.command == "reference":
        product = compute.load_product(nc_path)
        reference.build_cache(product["latitude"], product["longitude"], args.cache, args.downloads,
                              workers=args.workers)
        return 0
    rows = compute.point_metrics(nc_path, args.cache)
    compute.write_csv(args.output, rows)
    problems = compute.compare(args.data / "point_metrics.csv", rows)
    if problems:
        print(f"{len(problems)} values differ from the published point_metrics.csv; per column:")
        for line in compute.summarize(args.data / "point_metrics.csv", rows):
            print(f"  {line}")
        print(f"recomputed file written to {args.output}")
        return 1
    print(f"point_metrics.csv reproduced exactly ({len(rows)} points); written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
