#!/usr/bin/env bash
# The two length limits in .claude/CLAUDE.md ("Writing in this guide"), so the
# comments and NEWS do not regrow:
#   - every NEWS.md release entry is at most NEWS_MAX non-blank lines;
#   - in every header under inst/include/phylloptim/ and every R file under R/,
#     whole-line comments are at most COMMENT_MAX of the file's non-blank lines.
#     Roxygen lines (#') are documentation, not comments, and the generated
#     RcppR6 and RcppExports files are not ours to shorten.
# Plain awk, so cpp-tests.yml can run it on runners with no R.
#
# Usage, from the package root: tools/check_lengths.sh
set -u
NEWS_MAX=${NEWS_MAX:-12}
COMMENT_MAX=${COMMENT_MAX:-0.40}
status=0

awk -v max="$NEWS_MAX" '
  function flush() { if (head != "" && n > max) { printf "NEWS.md: \"%s\" is %d lines; the limit is %d\n", head, n, max; bad = 1 } }
  /^# / { flush(); head = substr($0, 3); n = 0; next }
  NF { n++ }
  END { flush(); exit bad }' NEWS.md || status=1

for f in inst/include/phylloptim/*.hpp R/*.R; do
  case "$f" in *RcppR6*|*RcppExports*) continue ;; esac
  awk -v max="$COMMENT_MAX" -v f="$f" '
    /^[ \t]*#'"'"'/ { next }
    NF { n++ }
    /^[ \t]*(\/\/|#)/ { c++ }
    END {
      if (n && c / n > max) {
        printf "%s: comments are %.0f%% of its lines; the limit is %.0f%%\n", f, 100 * c / n, 100 * max
        exit 1
      }
    }' "$f" || status=1
done

[ "$status" -eq 0 ] && echo "NEWS entries at most $NEWS_MAX lines; comments at most $(awk -v m="$COMMENT_MAX" 'BEGIN{printf "%.0f", 100*m}')% of every file"
exit $status
