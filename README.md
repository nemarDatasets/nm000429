[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000429-blue)](https://doi.org/10.82901/nemar.nm000429)

# Human subthalamic nucleus microelectrode recordings during an interrupted reach task (London et al., 2021)

Intraoperative microelectrode recordings from the subthalamic nucleus (STN) of patients with Parkinson's disease during
deep brain stimulation surgery (NYU Langone), while they performed an extended reaching task on a touchscreen in which the
target sometimes changed during the movement ("interrupted reach"). Article: subjects were off anti-Parkinsonian medication;
one or two Neuroprobe Sonus tungsten microelectrodes (1 MOhm, Alpha Omega), 2 mm apart, driven on one microdrive;
recordings referenced to the outer cannula; raw signals sampled at 44 kHz.

## Source
- Dryad: Dennis London, Arash Fazl, Kalman Katlowitz, Marisol Soula, Michael Pourfar, Alon Mogilner, Roozbeh Kiani. Distinct population code for movement kinematics and changes of ongoing movements in human
  subthalamic nucleus. doi:10.5061/dryad.2jm63xsq2. License CC0 1.0.
- Article: eLife 10:e64893 (2021), doi:10.7554/eLife.64893. Code: https://github.com/dlondon12/InterruptedReachClustMass
- The 8 subject zips were downloaded through the Dryad API on Voyager and matched the Dryad SHA-256 digests.
  Not mirrored: CWT_mean.zip, Hilbert_beta_mean.zip, PPC_*.zip (time-frequency and phase-consistency derivatives) and
  processed_data.mat (spiking/behaviour summary); they remain on Dryad.

## Contents
- `sub-<nn>/ieeg/sub-<nn>_task-interruptedreach_run-<block>_ieeg.*`: 16 task blocks of 8 subjects
  (2.57 h). `run` = the release block number. Channels: `CRAW_xx` = microelectrode tip, `CMacro_RAW_xx` = macro
  contact of the same lead, wideband RAW streams at the acquisition rate (44 kHz; Sub04: 22 kHz, rate from the stored
  sample count and begin/end times). Values: int16 x BitResolution / Gain (both stored per channel) in µV, float32.
- `events.tsv`: digital port words (CPORT__1) and digital-input rising edges (CDIG_IN_x) on the RAW time axis; see
  task-interruptedreach_events.json for the time convention. Trial-level behaviour (touch trajectories, trial types) is in
  the behaviour .mat files in sourcedata.
- `sourcedata/dryad-2jm63xsq2/SubXX/BlockYY/`: every file of the subject zips: the complete Alpha Omega export (also the
  online LFP streams at 1.375 kHz, spike segments, digital channels), the behaviour/task .mat and the Plexon spike-time
  .txt files (sorted units). Changes: MAT-file header dates set to day 01; file names with patient initials and a date
  were renamed (4 files: Sub02_OR_HandR_STNL_Block02.mat, Sub09_HandR_STNL_Block02.mat, plexon_Sub02_OR_HandR_STNL_Block02_spikes1.txt, plexon_Sub02_OR_HandR_STNL_Block02_spikes2.txt).

## Timing caveat
Sampling rate = the nominal rate stored per channel (_KHz). In some exports the stored begin/end times span more time than
samples / rate, which indicates recording pauses concatenated in the export; event onsets there are exact only up to the
first pause: Sub02/Block01: 29792526 samples over TimeBegin-TimeEnd 1634.2 s (18231.1 Hz vs nominal 44000 Hz); Sub04/Block01: 8873514 samples over TimeBegin-TimeEnd 172.2 s (51522.8 Hz vs nominal 22000 Hz); Sub11/Block01: 29672592 samples over TimeBegin-TimeEnd 133.8 s (221785.2 Hz vs nominal 44000 Hz); Sub11/Block05: 26669440 samples over TimeBegin-TimeEnd 139.9 s (190625.3 Hz vs nominal 44000 Hz); Sub11/Block06: 28558928 samples over TimeBegin-TimeEnd 138.0 s (207009.3 Hz vs nominal 44000 Hz); Sub11/Block07: 27740369 samples over TimeBegin-TimeEnd 135.2 s (205178.1 Hz vs nominal 44000 Hz); Sub13/Block02: 30603712 samples over TimeBegin-TimeEnd 139.1 s (219961.2 Hz vs nominal 44000 Hz).

## Task events from the behaviour files (Oct 2026)
`events.tsv` rows with `source_channel` = `behaviour` hold the trial events of the behaviour .mat files (sourcedata), placed on
the RAW time axis (code/fix-20261030_trial_events.py, report in code/fix-20261030_trial_events_report.json). The release has no
alignment code for this step (the authors' GitHub repository holds the statistics only), so the alignment was measured from the
files: port word 21002 marks each trial start; the sequence of 21002 times matches the behaviour `start_t` sequence at a single
trial offset in every block (median residual 1-7 ms, next candidate offset 210 ms or more), and a straight line maps behaviour time
to the RAW axis (slope about 0.99995; residuals up to 2 ms, except the first trial of sub-08 run-02/03, sub-09 run-01 and sub-11
run-01: 13-15 ms). Event times are `start_t + event_time`; the event words 21004 (fp_on), 21006 (targ_on), 21300 (targ_chg), 21101 (response / window entries) and
21015 (feedback) lie at a median 0.2-1.2 ms from the mapped times (per row: column `port_check_ms`); fp_off has no word of
its own.
- Added for 13 runs (819 trials, 5,107 rows); behaviour trials outside the recording (before the first or after the last start word) are not listed.
  sub-08 run-02 / run-03 have 13 / 7 start words after the last behaviour trial; those trials have no behaviour rows.
- Runs with the timing caveat above (sub-04 run-01, sub-11 run-01/05/07, sub-13 run-02) list only trials that end before
  TimeEnd - TimeBegin; sub-11 run-06 has none.
- sub-02 run-01 / run-02 are not changed here: they match as well (129 / 159 trials, all paired), but another pending change
  edits those two files; they will be added on top of it.
- `condition` is `trial_type_names{trial_type}` as stored in the file ('1mov', '2mov', '1mov_2mov' occur); `result` and `rt` are
  the stored fields.

## Localisation
Recordings were from the STN (article); the side comes from the file names (STNL/STNR). No coordinates or depths are in
the release; electrodes.tsv lists the channels with x, y, z = n/a.

## Participants
Age/sex per patient are in the article's Supplementary file 1, which is not linked to subject numbers in the release; n/a here.

## Privacy
File names with initials and dates were renamed; MAT headers keep year and month only. String scan of the behaviour .mat
files: no names, dates or identifiers found.
