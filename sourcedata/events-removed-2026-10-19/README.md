# Events removed on 2026-10-19 (version 1.0.3)

Rows removed from two `events.tsv` files because they lie outside their recording and belong to the neighbouring block.
Each file here has the header and the removed rows exactly as they were released (byte-for-byte lines), at the same
relative path as the BIDS file they came from.

- `sub-02/ieeg/sub-02_task-interruptedreach_run-01_events.tsv`: 119 `ttl_up` rows (CDIG_IN_1), onsets 767.39-1426.52 s,
  after the end of the run-01 data (677.102864 s). They are the run-02 `ttl_up` rows at onset + 677.102887 s.
- `sub-02/ieeg/sub-02_task-interruptedreach_run-02_events.tsv`: 99 `ttl_up` rows (CDIG_IN_1) with negative onsets
  (-582.32 to -99.58 s), before the start of the run-02 data. They are the run-01 `ttl_up` rows.

See README.md, section "Digital-input events of the other block (corrected 2026-10-19, version 1.0.3)".
