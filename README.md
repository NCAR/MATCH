# MATCH

Model for Atmospheric Transport and Chemistry

MATCH is a global offline chemical transport model used within the
NASA CERES project to produce aerosol datasets. It is driven by NCEP
reanalysis meteorology and assimilates MODIS aerosol optical depth
observations.

## Documentation

- [Architecture](docs/architecture.md) — Code organization, modules, and data flow
- [Building](docs/build.md) — Compiler setup, Makefile configuration, compilation
- [Running](RUN.md) — Run directory setup, namelist reference, emission tuning, CDAS↔CORe migration
- [AOD tuning plan](AOD_TUNE.md) — Per-month tuning of dust, sea salt, sulfate, and carbon emission scales
