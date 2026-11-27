#!/usr/bin/env bash
# Smoke test of a built `briareus` (build.jl) as a user without Julia would run it: the
# launcher alone, with no Julia on PATH and an empty depot. It writes a starter design
# with `briareus run`, runs the headless commands on it and checks their files, and
# prints the wall time of each command.
#
#     dev/release/smoke.sh build/Briareus [work dir]
set -euo pipefail
app="$(cd "$1" && pwd)"
work="${2:-$(mktemp -d)}"
mkdir -p "$work"
exe="$app/bin/briareus"
[ -x "$exe" ] || exe="$app/bin/briareus.exe"
here="$(cd "$(dirname "$0")" && pwd)"
export JULIA_DEPOT_PATH="$work/depot"
export PATH="$app/bin:/usr/bin:/bin"
# QML creates a Qt application as the system image loads, which stops the process on a
# Linux machine without a display unless Qt draws offscreen.
if [ "$(uname)" = Linux ] && [ -z "${DISPLAY:-}" ]; then export QT_QPA_PLATFORM=offscreen; fi
timed() { local t0=$SECONDS s=0; "$@" || s=$?; echo "  $(( SECONDS - t0 )) s: ${*:2}"; return $s; }
timed "$exe" --help > /dev/null
timed "$exe" run "$here/fixture.jl" "$work/design.briar"
timed "$exe" talis "$work/design.briar" --out "$work/talis.csv" || [ $? -eq 2 ]
timed "$exe" review-1w "$work/design.briar" --out "$work/review"
timed "$exe" roe "$work/design.briar" --out "$work/roe"
for f in design.briar talis.csv review/review-1w.html roe/roe.html roe/ports.svg; do
    [ -s "$work/$f" ] || { echo "missing: $f"; exit 1; }
done
echo "smoke test passed in $work"
