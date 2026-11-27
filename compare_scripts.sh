#!/bin/bash
# compare_scripts.sh -- compare the scripts in a laptop folder (searched recursively)
# with those in the local copy of the repository, by file name and by content.
#
# Usage:
#   bash compare_scripts.sh <laptop-folder> <repo-folder> [YYYY-MM-DD]
#
# The optional date (default 2026-07-01) limits the "only on laptop" list to files
# modified on or after that date, so that very old files do not flood the report.
# Compared file types: .py .for .f90 .sh .tex .bib.  Folders named .git, old,
# __pycache__ are skipped.  Works with the bash 3.2 and the md5 tool of macOS.

# Use the system tools only: no aliases or shell functions from the user's settings,
# and on macOS no GNU versions (Homebrew) that shadow the macOS ones in PATH.
unalias -a 2>/dev/null
while read -r _ _ F; do unset -f "$F"; done < <(declare -F)
[ "$(/usr/bin/uname)" = Darwin ] && PATH=/usr/bin:/bin:/usr/sbin:/sbin

LAP="$1"; REPO="$2"; SINCE="${3:-2026-07-01}"
if [ ! -d "$LAP" ] || [ ! -d "$REPO" ]; then
  echo "Usage: bash compare_scripts.sh <laptop-folder> <repo-folder> [YYYY-MM-DD]"; exit 1
fi
if command -v md5 >/dev/null 2>&1; then H() { md5 -q "$1"; }; else H() { md5sum "$1" | cut -d' ' -f1; }; fi
if stat -f %Sm /dev/null >/dev/null 2>&1; then
  D() { stat -f '%Sm' -t '%Y-%m-%d' "$1"; }           # macOS
else
  D() { date -r "$1" +%Y-%m-%d; }                      # Linux
fi

FIND='-name .git -prune -o -name old -prune -o -name __pycache__ -prune -o -type f
      ( -name *.py -o -name *.for -o -name *.f90 -o -name *.sh -o -name *.tex -o -name *.bib )
      -print'
TMP=$(mktemp -d)
set -f
find "$REPO" $FIND > "$TMP/repo"
find "$LAP"  $FIND > "$TMP/lap"
set +f

SAME=""; DIFF=""; ONLY=""
while IFS= read -r L; do
  B=$(basename "$L")
  R=$(grep "/$B\$" "$TMP/repo" | head -1)
  if [ -z "$R" ]; then
    [ ! "$(D "$L")" \< "$SINCE" ] && ONLY="$ONLY
  $(D "$L")  $L"
  elif [ "$(H "$L")" = "$(H "$R")" ]; then
    SAME="$SAME
  $B"
  else
    DIFF="$DIFF
  laptop $(D "$L")  $L
  repo            ${R#$REPO/}"
  fi
done < "$TMP/lap"

echo "=== Identical in both (by content):"; echo "$SAME" | sort -u | sed '/^$/d'
echo; echo "=== Same name, DIFFERENT content (check which is newer):"; echo "$DIFF" | sed '/^$/d'
echo; echo "=== Only on the laptop, modified since $SINCE:"; echo "$ONLY" | sort | sed '/^$/d'

echo; echo "=== Same file name in several laptop places (possible duplicate versions):"
sed 's|.*/||' "$TMP/lap" | sort | uniq -d | while read -r B; do
  grep "/$B\$" "$TMP/lap" | while read -r L; do echo "  $(D "$L")  $L"; done
done
rm -rf "$TMP"
