<!---
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
--->

# Run the first 2 cells:


_Cell 1: from jupyterlab template_: run it

```python jupyter={"source_hidden": true}
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
# Workflow to producing a collections of pdbs from sampled microstates
#### 5 steps to pdbs!

### 1. Necessary imports
```
from pathlib import Path
import numpy as np
import time         # only needed if you want to time a process

import base
import mcce_io as io
import ms_sampling as sampling
```
### 2. Path definition
 * mcce_output_path: path to  where a MCCE simulation was run. Must include step2_out.pdb, head3.lst, ms_out dir.

### 3. MS class instanciation:
 * The class needs values for pH and Eh in addition to the MCCE output path.
 * Call using variables:
```
pH, Eh = 5.0, 0.0   # can be 5, 0 (int) as well
ms = base.MS(mcce_output_path, pH, Eh)
```

### 4. Define arguments for sampling and pdb creation:
 - sample size
 - "sort by" key
 - output folder (optional: if not given, the output is given by ms.msout_file_dir)

#### Note: Rationale for using a folder created from the msout file, e.g. pH5eH0ms/:
The pdb file names only have the chosen MC and selected ms index as identifiers, hence,
a file must be open them to obtain the T, pH and Eh information, so keeping them in a folder
named after the msout file they come from is the simplest way to keep things tidy.

```
n_sample_size = 4
ms_sort_by = "energy"

# optional:
output_dir = some_folder_path  # defaults to ms.msout_file_dir if not given
```

### 5. Function call to create pdbs from sampled ms:
```
start_time = time.time()        # optional

sampling.pdbs_from_ms_samples(ms,
                              mcce_output_path,
                              n_sample_size,
                              ms_sort_by,
                              output_dir,
                              clear_pdbs_folder=True,  # default:True
                              list_files=True          # default:False
                            )

# next 2 lines: # optional
d = time.time() - start_time
print(f"`ms_sampling.pdbs_from_ms_samples` with sample size={n_sample_size:,} took {d/60:.2f} mins or {d:.2f} seconds")
```


---
# Example using repo data

```python
import base
import mcce_io as io
import ms_sampling as sampling
```

```python
mcce_output_path = Path.cwd().parent.joinpath("tests/data")
mcce_output_path

!ls {mcce_output_path}
```

# base.MC class

```python
print(base.MS.__doc__)
print(base.MS.__init__.__doc__)
```

```python
# create instance

pH, Eh = 5.0, 0.0

start_time = time.time()

ms = base.MS(mcce_output_path, pH, Eh)

d = time.time() - start_time
print(f"Loading of base.MS instance took {d/60:.2f} mins or {d:.2f} seconds")
print(ms)
```

```python
# Public vars in MC:  (uncomment to view)
#fdir(ms)
```

<!-- #raw -->
print(f"ms.counts: {ms.counts:,}") # == N_MONTERUNS! redundant
<!-- #endraw -->

# Call to `ms_sampling.pdbs_from_ms_samples`

```python
#fdir(sampling)  (uncomment to view)
```

```python
# create pdbs from samples ms

n_sample_size = 4
ms_sort_by = "energy"

start_time = time.time()

sampling.pdbs_from_ms_samples(ms,
                              mcce_output_path,
                              n_sample_size,
                              ms_sort_by,
                              clear_pdbs_folder=True,  # default:True
                              list_files=True          # default:False
                            )

d = time.time() - start_time
print(f"`ms_sampling.pdbs_from_ms_samples` with sample size={n_sample_size:,} took {d/60:.2f} mins or {d:.2f} seconds")
```

# Inspect a pdb head:

```python
!head -n 20 ../tests/data/ms_out/pH5eH0ms/pdbs_from_ms/mc0_ms1.pdb
```

```python

```
