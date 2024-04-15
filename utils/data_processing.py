import pandas as pd
import numpy as np
from typing import Union
from sklearn.metrics import r2_score
import statsmodels.api as sm
from tabulate import tabulate
import bisect


def linear_regression(
    X: pd.DataFrame, y: pd.Series, print_whole_table: bool = False
) -> list:
    X_np = np.array(X, dtype=float)
    y_np = np.array(y, dtype=float)
    mod = sm.OLS(y_np, X_np)
    fii = mod.fit()

    if print_whole_table:
        print(fii.summary2())
    regression_coeficitents = fii.summary2().tables[1]["Coef."].tolist()

    data = []
    for i, j in zip(regression_coeficitents, X.columns.tolist()):
        data.append([j, round(i, 4)])

    print(tabulate(data, headers=["coeficient"], tablefmt="double_outline"))
    print(
        "Pearson correlation coeficient "
        f"= {round(np.sqrt(float(fii.summary2().tables[0][1][6])),3)}"
    )

    return regression_coeficitents


def correlation(
    x: Union[list, np.array, pd.Series],
    y: Union[list, np.array, pd.Series],
) -> None:
    pearson_corr = np.sqrt(r2_score(x, y))
    print(f"Pearson correlation coeficient = {round(pearson_corr,3)}")


def return_bin_column(
    df: pd.DataFrame,
    column_to_bin: str,
    bins: list,
    name_of_new_column: str = "binned_column",
) -> pd.DataFrame:
    df[name_of_new_column] = df[column_to_bin].apply(
        lambda i: bins[int(bisect.bisect_left(bins, i)) - 1]
    )
    return df


def num_bins(arr: list, tolerance: float) -> int:
    """Bins values by the bin size (tolerance)
    and return number of bins.

    Args:
        arr (numpy array): array of peaks wavenumbers
        tolerance (float): bin size

    Returns:
        int: Number of bins
    """
    import numpy as np

    if isinstance(arr, list):
        arr = np.array(arr)
    binned = arr[~(np.triu(np.abs(arr[:, None] - arr) <= tolerance, 1)).any(0)]
    return int(len(binned))


def get_super(x: str) -> str:
    normal = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    )
    super_s = (
        "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    )
    res = x.maketrans("".join(normal), "".join(super_s))
    return x.translate(res)


def process_nmr(
    df_local: np.array, min_ppm: float, max_ppm: float, bin_size: float
) -> int:
    # transform to np.array if not already
    ppm_local = np.array(df_local)

    # select from range of chemical shifts
    bool_store = (ppm_local > min_ppm) & (ppm_local <= max_ppm)
    return num_bins(ppm_local[bool_store], bin_size)


def process_ir(
    df_local: np.array,
    min_wavenumber: float,
    max_wavenumber: float,
    ir_intensity_threshold: float,
    bin_size: float,
) -> int:
    # split intensity and wavenumber
    wn_arr_local = df_local[0]
    in_arr_local = df_local[1]

    # select from range of wavenumbers
    bool1_local = (wn_arr_local > min_wavenumber) & (
        wn_arr_local < max_wavenumber
    )

    # select above certain treshold
    bool3_local = in_arr_local >= ir_intensity_threshold
    return num_bins(wn_arr_local[bool1_local & bool3_local], bin_size)
