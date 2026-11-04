# ecohydro_feedback

Code and figure notebooks for:

> L. Cultra et al., **Ecohydrological Controls on Moist Convection and Long-Term Rainfall Feedback**, arXiv:2604.08442 (2026).

## Files

| File | Contents |
|---|---|
| `ParamsECO1.py` | Model parameters: soil types (`soil_params`), stomatal/atmospheric, plant hydraulic, carbon assimilation and boundary-layer constants |
| `Functionsdiurnal1.py` | Plant hydraulics and surface fluxes: soil and plant conductances, stomatal conductance, leaf energy balance, transpiration, initial boundary-layer profiles |
| `Atmosphericdynamics1.py` | Atmospheric column: stochastic rainfall, LCL height, moist adiabat, CAPE and convective precipitation (`RAINC` is current) |
| `PrecipGraph.ipynb` | CAPE and precipitation over near-surface temperature and humidity (Figure 5) |
| `Figures6and7_Forest_Grass.ipynb` | Figures 6 and 7: soil-moisture time series, CAPE and soil moisture at crossing time, and event frequency for the forest and grass (37°N) cases |

The `.py` files import each other by name, so keep them in the same folder as the notebooks.

## Requirements

```bash
pip install numpy scipy matplotlib jupyterlab
```

## Running

```bash
jupyter lab
```

- `PrecipGraph.ipynb` runs from the scripts alone.
- `Figures6and7_Forest_Grass.ipynb` plots saved simulation output (`forest.npz`, `grass.npz`, containing `ss`, `ssCAC` and `CAC`) and does not rerun the simulations. These files are not yet included in this repository.

## Citation

```bibtex
@article{cultra2026ecohydrological,
  title   = {Ecohydrological Controls on Moist Convection and Long-Term Rainfall Feedback},
  author  = {Cultra, Libby and others},
  journal = {arXiv preprint arXiv:2604.08442},
  year    = {2026}
}
```
