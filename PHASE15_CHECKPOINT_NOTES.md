# Phase15 exact-context recovery

Each completed unrestrained production segment and umbrella window now saves
an atomic archive containing an OpenMM binary checkpoint, raw samples and an
identity manifest. Source, scientific configuration, System, Integrator,
OpenMM version and platform must match before explicit `--resume-allostery`.
Binary restoration retains velocities and stochastic integrator state; no
portable-State fallback is presented as exact continuation. Sampling remains
833,000 steps per state. An interrupted unfinished window restarts from the
last completed checkpoint; this is not a per-step journal.

Eight regression tests passed, including real Reference-platform random-state
continuation and CPU-platform restoration, corruption rejection, atomic-write
failure and sample/step-count validation. These tests do not validate the full
free-energy calculation. The interrupted legacy trajectory has no such binary
checkpoint and cannot be resumed from its coordinate/console records. A new
full allostery run is still required. Existing completed spin outputs may only
be reused with separately recorded provenance.
