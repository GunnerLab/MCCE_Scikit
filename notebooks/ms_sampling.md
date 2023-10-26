---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:mce]
    language: python
    name: conda-env-mce-py
---

_Cell 1: from jupyterlab template_: run it

```python
import sys
from pathlib import Path
import time
import numpy as np
from pprint import pprint as pp
import matplotlib as mpl
from matplotlib import pyplot as plt
plt.ion()
#plt.style.use('seaborn-v0_8-muted')
#from IPython.display import HTML, Markdown #, IFrame

# To get multiple outputs into 1 cell w/o using print:
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# autoreload extension
from IPython import get_ipython

ipython = get_ipython()
if 'autoreload' not in ipython.extension_manager.loaded:
    %load_ext autoreload
%autoreload 2

# -----------------------------------------
# TWO USEFUL FUNCTIONS:

def add_to_sys_path(this_path, up=False):
    """
    Prepend this_path to sys.path.
    If up=True, path refers to parent folder (1 level up).
    """
    if up:
        newp = Path(this_path).parent
    else:
        newp = Path(this_path)

    src = newp.joinpath("src")
    if src.exists():
        newp = str(src)
    else:
        newp = str(newp)

    if newp not in sys.path:
        sys.path.insert(1, newp)
        print('Path added to sys.path: {}'.format(newp))

# Filtered dir() for method discovery:
def fdir(obj, start_with_str='_', exclude=True):
    return [d for d in dir(obj) if not d.startswith(start_with_str) == exclude]

```

_Cell 2: from jupyterlab template_: Uncomment (and amend it) to enable import of local modules.

```python
# Insert current src dir into sys.path so that modules in ../src can be imported:
# CHANGE THIS IF NEEDED:

add_to_sys_path(Path.cwd(), up=True)
```

---

---

# MCCE - MS Sampling (using test data in ../tests/data/)
---

```python
import base
import mcce_io as io
import ms_sampling as sampling
```

```python
DATA = Path.cwd().parent.joinpath("tests/data")
DATA
```

```python
!ls {DATA}
```

```python
# filepaths of inputs used by MS class:
h3_path = DATA.joinpath("head3.lst")
mcce_output_path = h3_path.parent
mcce_output_path

step2_path = mcce_output_path.joinpath("step2_out.pdb")
msout_path = mcce_output_path.joinpath("ms_out")
msout_path, msout_path.is_dir()
```

```python
# filepaths of outputs:

pH = 5.0
Eh= 0.0
msout_file = io.get_msout_filename(mcce_output_path, pH, Eh)
msout_file

msout_file_dir = msout_file.parent.joinpath(msout_file.stem)
msout_file_dir
```

```python
start_time = time.time()

io.split_msout_file(mcce_output_path, pH, Eh)

end_time = time.time()
print("io.divide_msout_file() took {:.2f} mins".format((end_time - start_time)/60))
```

```python
!ls {msout_file_dir}
```

```python
pdbs_dir = msout_file_dir.joinpath("pdbs_from_ms")
```

```python
!ls {pdbs_dir}
```

```python
io.clear_folder(pdbs_dir)
!ls {pdbs_dir}
```

# base.MC class

```python
print(base.MS.__doc__)
print(base.MS.__init__.__doc__)
```

```python
# create instance
start_time = time.time()

ms = base.MS(mcce_output_path, pH, Eh)

d = time.time() - start_time
print(f"Loading of base.MS instance took {d/60:.2f} mins or {d:.2f} seconds")
print(ms)
```

```python
# Public vars in MC:
fdir(ms)
```

```python
ms.method
```

# ms sampling

```python
fdir(sampling)
```

```python
n_sample_size = 5
ms_sort_by = "energy"
output_dir = msout_file_dir
```

```python
# create pdbs from samples ms
start_time = time.time()

sampling.pdbs_from_ms_samples(ms,
                              mcce_output_path,
                              n_sample_size,
                              ms_sort_by,
                              output_dir,
                              list_files=True)

d = time.time() - start_time
print(f"`sampling.pdbs_from_ms_samples` with sample size={n_sample_size:,} took {d/60:.2f} mins or {d:.2f} seconds")
```

```python
!head -n 20 ../tests/data/ms_out/pH5eH0ms/pdbs_from_ms/mc0_ms1.pdb
```

```python

```
