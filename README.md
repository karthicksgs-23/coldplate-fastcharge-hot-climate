# Cold-plate liquid cooling of a prismatic battery module during DC fast charging in hot climates

Python code and results for a thermal study of a bottom-cooled module of twelve 50 Ah prismatic NMC cells on an
aluminium cold plate, charged from 10 % to 80 % SOC at 1C–3C with a 45 °C current-derating rule. The study covers:

- a layered thermal network (cell layers, TIM, plate nodes, quasi-steady coolant segments), solved by backward Euler;
- laminar thermal-entry heat transfer (Shah & London) shifted to the flat-channel asymptote, Gnielinski above Re 2300;
- single-pass versus U-turn (two-pass) channel layouts, flow 0.5–8 L/min, TIM conductivity 1–6 W/(m·K),
  coolant supply 20–44 °C;
- hourly daytime air temperature 2005–2023 (PVGIS-ERA5) at Jodhpur, Bhuj, Anantapur and Paradip, mapped onto
  radiator-only and chiller (30 °C / 25 °C) cooling, with a Carnot-fraction chiller COP;
- ±20 % one-at-a-time sensitivity and numerical verification (layers, time step, energy balance, exact steady solution).

The cell parameters are generic values for the VDA PHEV2 format, not a measured product.

## Layout

| File | Purpose |
|---|---|
| `lc_model.py` | module model: `simulate()`, `hydraulics()`, `seg_hA()`, `chiller_cop()`; parameters in `P` |
| `run_study.py` | design matrix (450 runs), climate grid, site mapping, sensitivity, hydraulics, verification -> `results/study.json`, `results/traces.npz` |
| `run_extra.py` | 1 K supply-temperature sweep for 2C/3C and three TIMs -> `results/ts_sweep.json` |
| `make_schematic.py` | Fig. 1 (layout and thermal network) |
| `make_figures.py` | Figs. 2–6 |
| `make_tables.py` | formatted tables -> `results/tables.xlsx` |
| `data/climate_*.npz` | hourly air temperature (`T`, °C), irradiance (`G`), year, month, UTC hour, site latitude/longitude |
| `results/`, `figures/` | outputs of the scripts as used in the study |

## Running

```bash
pip install -r requirements.txt
python run_study.py      # about 8 min on 8 cores
python run_extra.py
python make_schematic.py
python make_figures.py
python make_tables.py
```

The scripts run from the repository root. Small dense solves are fastest single-threaded, so `lc_model.py` sets
`OMP_NUM_THREADS=1` (and the OpenBLAS/MKL equivalents) before importing NumPy.

## Data source

Climate series: PVGIS-ERA5 hourly data from the European Commission Joint Research Centre PVGIS service
(Huld et al., Sol. Energy 86 (2012) 1803–1815; Hersbach et al., Q. J. R. Meteorol. Soc. 146 (2020) 1999–2049).

## Limitations

Joule heat only (no entropic term); constant coolant properties except viscosity; equal flow in every channel;
the heat-transfer and friction correlations are verified numerically but not validated against a measured cold plate.

## License

MIT (see `LICENSE`).
