# Gannet tests

Unit and integration tests written with the MATLAB unit testing framework
(`matlab.unittest`).

| Test | What it checks |
| --- | --- |
| `SpectralRegistrationTest` | Spectral registration recovers known frequency and phase offsets in synthetic transients |
| `DataReaderTest` | `GERead`, `PhilipsRead`, and `SiemensTWIXRead` load Big GABA MEGA-PRESS data with the expected dimensions and parameters |
| `EndToEndTest` | `GannetLoad` and `GannetFit` run on a Big GABA Philips dataset and reproduce reference GABA+ estimates |

## Test data

The reader and end-to-end tests use subject S01 (TE = 68 ms) from the Big GABA
GE (G6), Philips (P3), and Siemens (S5) sites:

```text
GE_data/G6_MP/S01/S01_GABA_68.7
Philips_data/P3_MP/S01/S01_GABA_68_act.sdat (+ .spar)
Philips_data/P3_MP/S01/S01_GABA_68_ref.sdat (+ .spar)
Siemens_data/S5_MP/S01/S01_GABA_68.dat
Siemens_data/S5_MP/S01/S01_GABA_68_H2O.dat
```

Put these in `tests/data/` (ignored by Git), or point the `GANNET_TEST_DATA`
environment variable to a folder with the same layout. Tests that need missing
data are skipped, unless `GANNET_TEST_DATA_REQUIRED` is set to `1`, in which
case they fail. GitHub Actions downloads the data from the `test-data-v1`
release.

## Running the tests

From the Gannet folder in MATLAB:

```matlab
setenv('GANNET_TEST_DATA', '/path/to/test_data') % optional
results = runtests('tests')
```
