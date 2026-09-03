#!/bin/bash
# PreToolUse guard for CLAUDE.md 6.2/6.3 (no deletion of tracked VTU/PRJ under Tests/Data, or of any .md) and an early
# warning for 7 (no Claude co-author trailer; the authoritative check is .git/hooks/commit-msg). Written 2026-09-03 by an
# adversarial verifier, re-tested independently against 162 + 21 cases; residual bypasses recorded in the audit README.
c=$(jq -r ".tool_input.command // empty")
set -f
hit=""
# split on shell separators; each segment is inspected on its own so that text in
# echo/grep/commit messages/heredoc bodies of OTHER segments cannot trigger it
segs=$(printf '%s\n' "$c" | tr ';&|()`' '\n\n\n\n\n\n')
while IFS= read -r seg; do
  seg=$(printf '%s' "$seg" | sed -E 's/^[[:space:]]*((sudo|command|time|nohup|exec|do|then|else)[[:space:]]+)*//')
  verb=""
  case "$seg" in
    rm|rm[[:space:]]*|\\rm[[:space:]]*|unlink[[:space:]]*|shred[[:space:]]*|trash[[:space:]]*) verb=rm ;;
    /*) printf '%s' "$seg" | grep -Eq '^/[^[:space:]]*/rm([[:space:]]|$)' && verb=rm ;;
    git[[:space:]]*) printf '%s' "$seg" | grep -Eq '^git([[:space:]]+-[cC][[:space:]]+[^[:space:]]+)*[[:space:]]+rm([[:space:]]|$)' && verb=rm ;;
    find[[:space:]]*) printf '%s' "$seg" | grep -Eq -- '-delete|-exec[[:space:]]+(/[^[:space:]]*/)?rm[[:space:]]' && verb=find ;;
    xargs[[:space:]]*) printf '%s' "$seg" | grep -Eq '^xargs([[:space:]]+-[^[:space:]]+)*[[:space:]]+(/[^[:space:]]*/)?rm([[:space:]]|$)' && verb=xargs ;;
  esac
  [ -z "$verb" ] && continue
  rec=""; printf '%s' "$seg" | grep -Eq '(^|[[:space:]])(-[a-zA-Z]*[rR][a-zA-Z]*|--recursive)([[:space:]]|$)' && rec=1
  if [ "$verb" = find ] || [ "$verb" = xargs ]; then
    # cannot see the file list: block if the segment or the whole command names the guarded classes
    printf '%s' "$c" | grep -Eq 'Tests/Data|\.md' && hit="$seg"
    continue
  fi
  for tok in $(printf '%s' "$seg" | tr '"'"'"'' '  '); do
    case "$tok" in
      *.md|*.MD|*.markdown) hit="$tok" ;;
      *build*/Tests/Data*|*/__pycache__*|*Tests/Data/*/out*|*Tests/Data/*/_out*|*Tests/Data/*/Output*|*Tests/Data/*/results*) ;;   # output trees are free to delete
      *Tests/Data/*.vtu|*Tests/Data/*.prj|*Tests/Data/*.VTU|*Tests/Data/*.PRJ) hit="$tok" ;;
      *Tests/Data/*|*Tests/Data) [ -n "$rec" ] && hit="$tok" ;;
      \$*) printf '%s' "$c" | grep -Eq 'Tests/Data|\.md' && hit="$tok (variable-indirect)" ;;
    esac
  done
done <<EOF
$segs
EOF
if [ -n "$hit" ]; then echo "Blocked by CLAUDE.md 6.2/6.3 (matched: $hit): never delete tracked VTU/PRJ files or any .md file; annotate as historical instead" >&2; exit 2; fi
if printf '%s' "$c" | grep -Eq 'git([[:space:]]+-[cC][[:space:]]+[^[:space:]]+)*[[:space:]]+commit' \
   && printf '%s' "$c" | grep -Eiq '(^|[[:space:]]-[a-z]*m[[:space:]]+["'"'"']?|\\n|--trailer[= ]+["'"'"']?)[[:space:]]*co-authored-by[=:][[:space:]\\]*claude'; then
  echo "Blocked by CLAUDE.md 7: no Co-Authored-By: Claude trailer unless Vinay explicitly asked for one" >&2; exit 2
fi
exit 0
