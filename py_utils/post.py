import sys 
from pathlib import Path
from helpers import plot_latest_product
import matplotlib.pyplot as plt
import xarray as xr

data_dir = Path(sys.argv[1])
figs_dir = str(sys.argv[2])

extent = [float(x) for x in sys.argv[3:8]]

products = ["CAL", "DNI", "GHI", "CSR"]

for p in products: 
    plot_latest_product(p, data_dir/p, figs_dir, extent)