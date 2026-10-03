BT-DAT-Q022

**Built on** `PREREG.md` (sha256 `94db5ca1…d54f`, frozen 19:31:53Z before any download) and the interrupted session's `inputs/Q022_model.py` (sha256 `9d97bd4b…d635`, read-only, unchanged). The BodyTwin repo anchors named in `inputs/NIGHT_PREAMBLE.md` are not mounted here, so no node id is invented.

**What was found** — `DATA_SOURCES.json`, 10 public measured datasets, every licence and total size read off the repository page. Top four, all Open Data Commons / CC BY and downloadable without credentials: MIMIC-IV Waveform (ODbL, 12.8 GB) with ECG+Pleth+impedance on one time base; MIMIC-IV-ECG (ODbL, 90.4 GB, 12-lead matched to invasive pressure); HeartCycle (ODC-BY, 2.3 GB, simultaneous ICG+echo+ECG); SensSmartTech (CC BY 4.0, 775 MB, PCG+ECG+PPG+ACC over 58–173 bpm). PTB-XL (CC BY 4.0, 3.0 GB) constrains the electrical drive. MIMIC-IV-Echo is credentialed — not usable.

**Sample proof** — `samples/`, 6,736,362 B ≤ 50 MiB, sha256 per file in `results.json`. PTB-XL record 00001_hr parsed by `verify_sample.py`: 12 leads, 16-bit, gain 1000/mV, ADC zero 0, 500 Hz, 5000 samples = 120,000 B exactly. Coordinate frame **verified**: III=II−I, aVR=−(I+II)/2, aVL=I−II/2, aVF=II−I/2 all hold to ≤0.0028 relative rms, i.e. the 1 µV/LSB floor. Cohort from the downloaded CSV: 21,799 records, 52.08% male, median 62 y; NORM 9514, MI 5469, CD 4898 match the release exactly, STTC/HYP do not (shipped class map covers 44 codes) and are not claimed.

**Measured vs model** — median RR 0.652 s, HR 92.0 bpm, QRS 0.040 s. Model drive FWHM is 0.0589 s, i.e. **1.47× wider than measured**. C3 FAILED as preregistered: no public beat-resolved coronary flow dataset exists, so `gregg_gain`, `hypoxia_penalty`, `perfusion_gamma`, `flow_max`, `p_ext` stay literature-anchored. MIMIC-IV payloads are FLAC-in-WFDB; form verified, values UNKNOWN (no decoder, no pip).

**Next** — install a FLAC reader, then MIMIC-IV Waveform for the R-to-ABP lag, which is the only open test of the model's [0, 0.080] s and [0, 0.150] s windows.
