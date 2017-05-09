import numpy as np
import sys
from datetime import datetime

def time_print(*s):
#    import sys
#    from datetime import datetime
    print ' '.join(map(str, s)), datetime.today()
    sys.stdout.flush()

def load_table_as_list(f, delimiter='\t', encoding=None):    
    if encoding is None:
        f = open(f)
    else:
        import codecs
        f = codecs.open(f, encoding=encoding)
    tmp = [x.strip().split(delimiter) for x in f.read().splitlines()]
    f.close()
    return tmp    

def pivot_2_square(mat):
    """ Turn table into a symmetric array.
    Different sets of genes index the rows and columns. Make them the same sets and in the same order."""

    rows, cols = mat.index, mat.columns
    total_genes = np.unique(np.append(rows.values, cols.values))
    rows_needed = np.setdiff1d(total_genes, mat.index.unique())
    cols_needed = np.setdiff1d(total_genes, mat.columns.unique())
    assert rows.size + rows_needed.size == cols.size + cols_needed.size

    arr = mat.values.copy()
    assert arr.dtype == np.float32
    tmp = np.empty((rows_needed.size, arr.shape[1]), arr.dtype)
    tmp.fill(np.nan)
    arr = np.vstack([arr, tmp])
    tmp = np.empty((arr.shape[0], cols_needed.size), arr.dtype)
    tmp.fill(np.nan)
    arr = np.hstack([arr, tmp])
    assert arr.flags['C_CONTIGUOUS']

    rows, cols = np.append(rows.values, rows_needed), np.append(cols.values, cols_needed)
    rows_argsort, cols_argsort = np.argsort(rows), np.argsort(cols)
    rows, cols = rows[rows_argsort], cols[cols_argsort]
    assert np.all(rows==cols)

    arr = arr[rows_argsort, :][:, cols_argsort]
    arr[np.isnan(arr)] = np.inf
    arr = np.minimum(arr, arr.T)
    arr[np.isinf(arr)] = np.nan

    rows_index = {b : a for a, b in enumerate(rows)}

    return arr, rows, rows_index 

