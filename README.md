# Copepod Latitudinal Diversity Gradient

A data pipeline and statistical analysis that fetches global copepod occurrence records from OBIS, bins them by latitude, and fits a quadratic model to explore the **Latitudinal Diversity Gradient (LDG)**.

## Output

| File | Description |
|---|---|
| `copepod_quadratic_fit.png` | Side-by-side plot: quadratic fit + residual bar chart |
| `copepod_quadratic_fit.csv` | Per-bin table of observed S, predicted S, and residuals |

The fitted model takes the form:

$$S(\phi) = a\phi^2 + b\phi + c$$

where $\phi$ is latitude in degrees. The peak richness latitude is estimated at $\hat{\phi} = -b / 2a$.


## Requirements

```
Python >= 3.9
numpy
pandas
matplotlib
rose-pine-matplotlib
pyobis
```

Dependencies:

```bash
pip install numpy pandas matplotlib pyobis
# Rose Pine theme (optional, just for 'pretty-ness')
pip install rose-pine-matplotlib
```

##  Usage

```bash
python main.py
```


## Methodology

### Data Source
Species occurrence data is collected from OBIS using taxon ID `1080` (Copepoda), with a maximum of 50,000 records.

### Binning
Latitudes are discretised into **5° bins** spanning −90° to +90°. Species richness per bin is the count of unique species names within that band.

### Regression
A quadratic polynomial is fit using the **Normal Equations**:

$$\hat{\beta} = (X^TX)^{-1}X^Ty$$

Model diagnostics:

| Metric | Description |
|---|---|
| $R^2$ | Coefficient of determination |
| RMSE | Root mean squared error |
| SE | Standard errors of coefficients via the covariance matrix |


## Interpretation

Copepods — unlike terrestrial taxa — often show a **bimodal or reverse LDG**, with higher diversity at mid-to-high latitudes in productive upwelling zones. The quadratic fit captures the overall curvature, and the residual plot highlights where the model over- or under-predicts richness.


## Structure

```
.
├── copepod_lgd.py                 # Main analysis script
├── copepod_quadratic_fit.png      # Output figure
├── copepod_quadratic_fit.csv      # Output data table
└── README.md
```

## Acknowledgements

- Occurrence data sourced from [OBIS](https://obis.org/) via [`pyobis`](https://github.com/iobis/pyobis)
- Styled with [Rose Pine](https://rosepinetheme.com/) for Matplotlib
