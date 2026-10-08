#!/usr/bin/env python3
"""
Build a 1957-2002 archive of ERA-40 u and v at 300 hPa (N80 grid, 6-hourly)
from NCAR GDEX dataset d117001 without keeping the 1.3-1.5 GB monthly files.

What is in d117001: the ERA-40 pressure-level analyses as ECMWF archived them,
i.e. spectral (T159) fields.  The horizontal wind is stored as vorticity (vo)
and divergence (d), NOT as u and v.  So this script

  1. pulls the vo and d messages at 300 hPa out of each monthly file
     (2 fields x 4 times/day, ~26k spectral coefficients each) -> vod300_YYYYMM.grib
  2. converts them to u and v on the 320 x 160 regular Gaussian N80 grid
     (1.125 deg) with  cdo dv2uv,linear  -> uv300_YYYYMM.grib
     (this is the same transform MARS applies when you ask ECMWF for u/v)

Step 2 needs the `cdo` executable (in the conda env file); without it the
script stops after step 1 and you can run the cdo command yourself later.

Two download modes, chosen per month:

  full   download e4oper.an.pl.YYYYMM, copy out the wanted messages with
         ecCodes, delete the big file.
  range  fetch only the byte ranges that hold those messages, using HTTP
         Range requests.  Needs a "template" (built automatically from the
         first month done in full mode) and assumes every month has the same
         message order and the same per-timestep byte layout.  Every fetched
         message is validated (GRIB magic, shortName, level, validity time);
         if anything is off, the month is redone in full mode.

  auto   (default) first month full, then range with fallback to full.

Months whose output already exists and validates are skipped, so the script
can be re-run after an interruption.  If a downloaded file turns out not to
contain the expected fields, the script prints an inventory of what it does
contain and keeps the file.

Winter (DJF) version: --year Y --nwinters n fetches the months
    Y-12, (Y+1)-01, (Y+1)-02, (Y+1)-12, (Y+2)-01, (Y+2)-02, ..., (Y+n)-02
i.e. n consecutive December-January-February seasons starting with Dec Y.
Months outside the ERA-40 record (Sep 1957 - Aug 2002) are skipped.

Usage
    python fetch_era40_uv300_djf.py --year 1960 --nwinters 5 --out-dir era40_uv300
    python fetch_era40_uv300_djf.py --year 1957 --nwinters 45 --out-dir era40_uv300 --workers 3
    python fetch_era40_uv300_djf.py --year 1998 --nwinters 2 --downloader wget

Afterwards, e.g.
    cat era40_uv300/uv300_*.grib > uv300_1957-2002.grib
    cdo -f nc4 -z zip copy uv300_1957-2002.grib uv300_1957-2002.nc

Requirements: python-eccodes, requests, cdo   (all in the conda env file)
Optional:     wget, or the Pelican client (https://pelicanplatform.org) for
              --downloader wget|pelican
"""
import argparse
import calendar
import datetime as dt
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import requests
import eccodes as ec

BASE = "https://osdf-director.osg-htc.org/ncar/gdex/d117001"
FNAME = "e4oper.an.pl.{ym}"
WANT = {"vo", "d"}         # shortNames to keep (spectral vorticity + divergence)
UV = {"u", "v"}            # what cdo dv2uv produces
LEVEL = 300                # hPa
LEVTYPE = "isobaricInhPa"
STEPS_PER_DAY = 4
CHUNK = 1 << 20            # 1 MiB streaming chunk


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def months(start, end):
    y, m = int(start[:4]), int(start[4:])
    while True:
        ym = f"{y}{m:02d}"
        if ym > end:
            return
        yield ym
        m += 1
        if m == 13:
            y, m = y + 1, 1


ERA40_FIRST, ERA40_LAST = "195709", "200208"


def djf_months(year, nwinters):
    """Dec of `year` through Feb of `year + nwinters`, DJF seasons only."""
    out = []
    for k in range(nwinters):
        out += [f"{year + k}12", f"{year + k + 1}01", f"{year + k + 1}02"]
    skipped = [ym for ym in out if not ERA40_FIRST <= ym <= ERA40_LAST]
    if skipped:
        print("outside the ERA-40 record, skipped:", " ".join(skipped))
    return [ym for ym in out if ERA40_FIRST <= ym <= ERA40_LAST]


def ndays(ym):
    return calendar.monthrange(int(ym[:4]), int(ym[4:]))[1]


def expected_times(ym):
    """List of validity datetimes for a month, in file order."""
    d0 = dt.datetime(int(ym[:4]), int(ym[4:]), 1)
    return [d0 + dt.timedelta(hours=6 * i) for i in range(STEPS_PER_DAY * ndays(ym))]


def msg_meta(gid):
    d, t = ec.codes_get(gid, "validityDate"), ec.codes_get(gid, "validityTime")
    return dict(
        shortName=ec.codes_get(gid, "shortName"),
        typeOfLevel=ec.codes_get(gid, "typeOfLevel"),
        level=ec.codes_get(gid, "level"),
        valid=dt.datetime(d // 10000, d // 100 % 100, d % 100, t // 100, t % 100),
    )


def wanted(meta, names=None):
    return (meta["shortName"] in (names or WANT) and meta["level"] == LEVEL
            and meta["typeOfLevel"] == LEVTYPE)


def inventory(path, nmax=200000):
    """Count (shortName, typeOfLevel, gridType) combos in a GRIB file."""
    from collections import Counter
    c = Counter()
    with open(path, "rb") as f:
        for _ in range(nmax):
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            c[(ec.codes_get(gid, "shortName"), ec.codes_get(gid, "typeOfLevel"),
               ec.codes_get(gid, "gridType"))] += 1
            ec.codes_release(gid)
    return c


def validate_output(path, ym, verbose=False, names=None):
    """Check an output file holds exactly `names` at 300 hPa for every 6-h step."""
    names = names or WANT
    if not os.path.exists(path):
        return False
    seen = set()
    try:
        with open(path, "rb") as f:
            while True:
                gid = ec.codes_grib_new_from_file(f)
                if gid is None:
                    break
                m = msg_meta(gid)
                ec.codes_release(gid)
                if not wanted(m, names):
                    if verbose:
                        print(f"  {ym}: unexpected message {m}")
                    return False
                seen.add((m["shortName"], m["valid"]))
    except Exception as e:  # truncated / corrupt
        if verbose:
            print(f"  {ym}: cannot read output ({e})")
        return False
    need = {(v, t) for v in names for t in expected_times(ym)}
    ok = seen == need
    if verbose and not ok:
        print(f"  {ym}: have {len(seen)} messages, need {len(need)}")
    return ok


# ---------------------------------------------------------------------------
# full mode
# ---------------------------------------------------------------------------
class NoProgress(Exception):
    """Server answered but sent no (new) bytes."""


def resolve(session, url, log=None):
    """Follow the OSDF director redirect once and return the list of cache URLs
    to try, best first.  The director answers with a 307 whose Location is the
    preferred cache and whose Link header lists alternatives (rel="duplicate",
    pri=N).  A plain web server (no redirect) just returns [url]."""
    r = session.get(url, allow_redirects=False, stream=True, timeout=60)
    r.close()
    if r.status_code in (301, 302, 303, 307, 308) and "Location" in r.headers:
        urls = [r.headers["Location"]]
        for part in r.headers.get("Link", "").split(","):
            m = re.match(r'\s*<([^>]+)>(.*)', part)
            if m and 'rel="duplicate"' in m.group(2):
                pri = re.search(r"pri=(\d+)", m.group(2))
                urls.append((int(pri.group(1)) if pri else 99, m.group(1)))
        alts = [u for _, u in sorted(x for x in urls[1:] if isinstance(x, tuple))]
        urls = [urls[0]] + [u for u in alts if u != urls[0]]
        if log:
            log(f"  director offers {len(urls)} cache(s): " + ", ".join(
                re.sub(r"https?://([^/]+).*", r"\1", u) for u in urls))
        return urls
    r.raise_for_status()
    return [url]


def _stream_to_file(session, url, dest, log):
    """One streaming attempt with resume.  Raises NoProgress if nothing new
    arrived, re-raises connection errors otherwise (caller decides)."""
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    with session.get(url, headers=headers, stream=True, timeout=(30, 300)) as r:
        if r.status_code == 416:            # already complete
            return have
        r.raise_for_status()
        if have and r.status_code != 206:   # server ignored Range: start over
            have = 0
        total = int(r.headers.get("Content-Length", 0)) + have
        done, t0 = have, time.time()
        try:
            with open(dest, "ab" if have else "wb") as f:
                for chunk in r.iter_content(CHUNK):
                    f.write(chunk)
                    done += len(chunk)
                    if time.time() - t0 > 30:
                        log(f"  {done / 1e9:.2f}/{total / 1e9:.2f} GB")
                        t0 = time.time()
        except (requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout) as e:
            if done == have:
                raise NoProgress(str(e)[:120])
            raise
    return done


def download(session, url, dest, log, downloader="python", max_rounds=6):
    """Download url -> dest robustly.  Returns final size.

    downloader = "python"  : requests, with cache fail-over and resume
                 "wget"    : shell out to wget -c (resumes, retries itself)
                 "pelican" : shell out to the Pelican/OSDF client
    """
    if downloader == "wget":
        cmd = ["wget", "-c", "--tries=20", "--waitretry=60", "--retry-connrefused",
               "--read-timeout=300", "-q", "--show-progress", "-O", dest, url]
        subprocess.run(cmd, check=True)
        return os.path.getsize(dest)
    if downloader == "pelican":
        osdf = re.sub(r"^https?://[^/]+", "osdf://", url)
        subprocess.run(["pelican", "object", "get", osdf, dest], check=True)
        return os.path.getsize(dest)

    caches = resolve(session, url, log)
    for rnd in range(max_rounds):
        for cache in caches:
            host = re.sub(r"https?://([^/]+).*", r"\1", cache)
            before = os.path.getsize(dest) if os.path.exists(dest) else 0
            try:
                return _stream_to_file(session, cache, dest, log)
            except requests.exceptions.HTTPError as e:
                if e.response is not None and e.response.status_code in (403, 404):
                    raise RuntimeError(f"HTTP {e.response.status_code} for {url}: file not on server")
                log(f"  {host}: HTTP {e.response.status_code if e.response is not None else '?'}; trying next cache")
            except (NoProgress, requests.exceptions.ChunkedEncodingError,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.ReadTimeout) as e:
                after = os.path.getsize(dest) if os.path.exists(dest) else 0
                if after > before:
                    log(f"  {host}: interrupted after {after / 1e9:.2f} GB; resuming")
                else:
                    log(f"  {host}: no data ({str(e)[:90]}); trying next cache")
        wait = 60 * (rnd + 1)
        log(f"  all caches failed; waiting {wait}s for the OSDF cache to warm up (round {rnd + 1}/{max_rounds})")
        time.sleep(wait)
        caches = resolve(session, url)      # re-ask the director; it may route elsewhere
    raise RuntimeError("could not download after repeated attempts")


def extract_full(src, out, ym, template_path=None):
    """Copy wanted messages src -> out. Optionally record a byte-layout template."""
    layout = []            # (offset, length, shortName, time index) for wanted msgs
    first_valid = None
    with open(src, "rb") as fi, open(out, "wb") as fo:
        while True:
            gid = ec.codes_grib_new_from_file(fi)
            if gid is None:
                break
            m = msg_meta(gid)
            if first_valid is None:
                first_valid = m["valid"]
            if wanted(m):
                ec.codes_write(gid, fo)
                off = int(ec.codes_get(gid, "offset"))
                ln = int(ec.codes_get(gid, "totalLength"))
                tidx = int((m["valid"] - first_valid).total_seconds() // (6 * 3600))
                layout.append((off, ln, m["shortName"], tidx))
            ec.codes_release(gid)
    if not layout:
        inv = inventory(src)
        lines = "\n".join(f"    {n:>6} x {name:<6} {lt:<16} {gt}" for (name, lt, gt), n in inv.most_common(30))
        raise RuntimeError(f"no {sorted(WANT)} messages at {LEVEL} hPa in {os.path.basename(src)}; "
                           f"file kept. It contains:\n{lines}")
    if template_path:
        _write_template(layout, os.path.getsize(src), ym, template_path)


def _write_template(layout, fsize, ym, path):
    """Derive the per-timestep byte block from the full-file layout."""
    nt = STEPS_PER_DAY * ndays(ym)
    if fsize % nt:
        print(f"  template: file size {fsize} not divisible by {nt} steps; range mode disabled")
        return
    block = fsize // nt
    per_step = {}
    for off, ln, name, tidx in layout:
        rel = off - tidx * block
        if not 0 <= rel < block:
            print("  template: message outside its time block; range mode disabled")
            return
        per_step.setdefault(tidx, []).append((rel, ln, name))
    steps = list(per_step.values())
    if any(s != steps[0] for s in steps):
        print("  template: layout differs between time steps; range mode disabled")
        return
    json.dump({"block": block, "msgs": steps[0], "from": ym}, open(path, "w"), indent=1)
    print(f"  template written: {len(steps[0])} msgs/step, block {block} bytes")


# ---------------------------------------------------------------------------
# range mode
# ---------------------------------------------------------------------------
class CacheUnavailable(Exception):
    """No cache would serve the request (network trouble, not a layout problem)."""


def fetch_range(session, caches, start, end, rounds=3):
    """GET bytes start-end, trying each cache in turn; raises CacheUnavailable."""
    for rnd in range(rounds):
        for url in caches:
            try:
                r = session.get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=(30, 120))
                if r.status_code != 206:
                    raise RuntimeError(f"range request not honoured (HTTP {r.status_code})")
                if len(r.content) != end - start + 1:
                    raise requests.exceptions.ConnectionError("short body")
                return r.content
            except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError,
                    requests.exceptions.ReadTimeout):
                continue
        time.sleep(30 * (rnd + 1))
    raise CacheUnavailable(f"no cache served bytes {start}-{end}")


def extract_range(session, url, out, ym, tpl, log):
    """Fetch only the wanted messages. Raises on any inconsistency."""
    nt = STEPS_PER_DAY * ndays(ym)
    block, msgs = tpl["block"], sorted(tpl["msgs"])
    caches = resolve(session, url)           # talk to the caches directly
    size = None
    for c in caches:
        try:
            size = int(session.head(c, allow_redirects=True, timeout=60).headers["Content-Length"])
            break
        except (requests.exceptions.RequestException, KeyError, ValueError):
            continue
    if size is None:
        raise CacheUnavailable("no cache answered HEAD")
    if size != nt * block:
        raise RuntimeError(f"size {size} != {nt} x {block}")
    times = expected_times(ym)
    # merge the per-step ranges into as few HTTP requests as possible
    lo, hi = msgs[0][0], max(r + ln for r, ln, _ in msgs)
    with open(out, "wb") as fo:
        for t in range(nt):
            buf = fetch_range(session, caches, t * block + lo, t * block + hi - 1)
            for rel, ln, name in msgs:
                raw = buf[rel - lo: rel - lo + ln]
                if raw[:4] != b"GRIB" or raw[-4:] != b"7777":
                    raise RuntimeError(f"step {t}: bad GRIB framing")
                gid = ec.codes_new_from_message(raw)
                m = msg_meta(gid)
                ec.codes_release(gid)
                if m["shortName"] != name or not wanted(m) or m["valid"] != times[t]:
                    raise RuntimeError(f"step {t}: got {m}, expected {name} @ {times[t]}")
                fo.write(raw)
            if t % 40 == 0:
                log(f"  step {t}/{nt}")


# ---------------------------------------------------------------------------
# vorticity/divergence -> u/v
# ---------------------------------------------------------------------------
def have_cdo():
    from shutil import which
    return which("cdo") is not None


def vod_to_uv(vod, uv, ym, log):
    """Spectral vo,d -> u,v on the linear (N80, 320x160) Gaussian grid."""
    cmd = ["cdo", "-s", "-f", "grb", "dv2uv,linear", vod, uv + ".part"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"cdo failed: {r.stderr.strip()[:300]}")
    os.replace(uv + ".part", uv)
    if not validate_output(uv, ym, verbose=True, names=UV):
        raise RuntimeError("u/v output failed validation")
    log("u/v on N80 grid written")


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def process_month(ym, args, tpl):
    out = os.path.join(args.out_dir, f"vod300_{ym}.grib")      # spectral vo, d
    uv = os.path.join(args.out_dir, f"uv300_{ym}.grib")        # gridded u, v
    tmp = os.path.join(args.out_dir, FNAME.format(ym=ym))
    url = f"{BASE}/{FNAME.format(ym=ym)}"
    log = lambda s: print(f"[{ym}] {s}", flush=True)

    def finish(how):
        if args.no_uv:
            return how
        if validate_output(uv, ym, names=UV):
            return how
        if not have_cdo():
            log("cdo not found: vod300 file written, u/v conversion skipped")
            return how
        vod_to_uv(out, uv, ym, log)
        return how

    if validate_output(out, ym):
        log("already downloaded")
        return finish("skip")

    session = requests.Session()
    mode = args.mode
    if mode == "auto":
        mode = "range" if tpl else "full"

    if mode == "range":
        try:
            log("range fetch")
            extract_range(session, url, out + ".part", ym, tpl, log)
            os.replace(out + ".part", out)
            if validate_output(out, ym, verbose=True):
                log("ok (range)")
                return finish("range")
            raise RuntimeError("validation failed")
        except Exception as e:
            log(f"range mode failed ({e}); falling back to full download")
            for p in (out, out + ".part"):
                if os.path.exists(p):
                    os.remove(p)

    log(f"full download ({args.downloader})")
    for attempt in range(1, 4):
        try:
            download(session, url, tmp, log, downloader=args.downloader)
            break
        except RuntimeError as e:          # 404 etc.: no point retrying
            log(f"download failed: {e}")
            return "failed"
        except Exception as e:
            log(f"download attempt {attempt} failed: {str(e)[:160]}")
            time.sleep(120 * attempt)
    else:
        return "failed"
    tpl_path = os.path.join(args.out_dir, "template.json")
    try:
        extract_full(tmp, out + ".part", ym, None if os.path.exists(tpl_path) else tpl_path)
    except RuntimeError as e:
        log(str(e))
        return "failed"
    os.replace(out + ".part", out)
    if not validate_output(out, ym, verbose=True):
        log("output failed validation; keeping the full file for inspection")
        return "failed"
    if not args.keep:
        os.remove(tmp)
    log("ok (full)")
    return finish("full")


def main():
    global BASE, WANT
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default="era40_uv300")
    p.add_argument("--year", type=int, required=True, help="Y: December of this year starts the first winter")
    p.add_argument("--nwinters", type=int, required=True, help="n: number of DJF seasons (last month is Feb Y+n)")
    p.add_argument("--mode", choices=["auto", "full", "range"], default="auto")
    p.add_argument("--workers", type=int, default=1, help="months processed in parallel")
    p.add_argument("--keep", action="store_true", help="do not delete the full monthly files")
    p.add_argument("--no-uv", action="store_true", help="stop after extracting vo/d (skip the cdo step)")
    p.add_argument("--fields", default="vo,d",
                   help="comma-separated shortNames to extract (default vo,d; use u,v if the "
                        "files turn out to hold gridded winds)")
    p.add_argument("--downloader", choices=["python", "wget", "pelican"], default="python",
                   help="how to fetch the full monthly files (default: python requests with "
                        "OSDF cache fail-over; wget and pelican must be on PATH)")
    p.add_argument("--base-url", default=BASE, help=argparse.SUPPRESS)
    args = p.parse_args()
    BASE = args.base_url
    WANT = set(args.fields.split(","))
    if WANT != {"vo", "d"}:
        args.no_uv = True          # cdo step only makes sense for vo/d
    os.makedirs(args.out_dir, exist_ok=True)

    todo = djf_months(args.year, args.nwinters)
    print(f"{len(todo)} months to fetch: {todo[0]} .. {todo[-1]}" if todo else "nothing to fetch")
    if not todo:
        return
    tpl_path = os.path.join(args.out_dir, "template.json")
    tpl = json.load(open(tpl_path)) if os.path.exists(tpl_path) else None

    # months are done one at a time until a template exists (needs one full
    # download); after that the rest can run in parallel in range mode
    results = {}
    while todo and tpl is None and args.mode != "full":
        results[todo[0]] = process_month(todo[0], args, None)
        todo = todo[1:]
        tpl = json.load(open(tpl_path)) if os.path.exists(tpl_path) else None
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for ym, res in zip(todo, ex.map(lambda ym: process_month(ym, args, tpl), todo)):
            results[ym] = res

    counts = {k: sum(1 for v in results.values() if v == k) for k in ("skip", "range", "full", "failed")}
    print("\nsummary:", counts)
    failed = [ym for ym, v in results.items() if v == "failed"]
    if failed:
        print("failed months (re-run to retry):", " ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()