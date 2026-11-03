#!/usr/bin/env bash
# Keep run_all.sh going until all 9 runs finish.
# Every CHECK_SECONDS: if run_all.sh is not running and progress.log doesn't
# say the series completed, restart it (it resumes from checkpoints / skips
# finished steps). Gives up after MAX_SAME restarts that stop at the same
# step, so a deterministic error can't restart forever -- the log says where.
#
# Usage (from repo root):  nohup bash prompt_afterglow_runs/watchdog.sh > /dev/null 2>&1 &
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"
PA="$REPO/prompt_afterglow_runs"
LOG="$PA/progress.log"
WLOG="$PA/watchdog.log"
CHECK_SECONDS=600
MAX_SAME=3

wlog() { echo "$(date '+%Y-%m-%d %H:%M:%S')  watchdog: $*" | tee -a "$WLOG" >> "$LOG"; }

last_stop=""; same=0
wlog "started (pid $$), checking every ${CHECK_SECONDS}s"
while true; do
  if grep -q "=== all 9 runs complete ===" "$LOG" 2>/dev/null; then
    wlog "series complete -- exiting"; exit 0
  fi
  if ! pgrep -f "prompt_afterglow_runs/run_all.sh" > /dev/null; then
    stop="$(grep -E '\] (formula search|SuperLearner) (START|FAILED)' "$LOG" | tail -1 | cut -c22-)"
    if [ "$stop" = "$last_stop" ]; then same=$((same + 1)); else same=1; last_stop="$stop"; fi
    if [ "$same" -gt "$MAX_SAME" ]; then
      wlog "stopped at the same step $MAX_SAME times in a row ($stop) -- giving up, needs a look"; exit 1
    fi
    wlog "run_all.sh not running (last step: $stop) -- restarting (attempt $same/$MAX_SAME at this step)"
    nohup caffeinate -i bash "$PA/run_all.sh" >> "$PA/run_all.nohup.log" 2>&1 &
  fi
  sleep "$CHECK_SECONDS"
done
