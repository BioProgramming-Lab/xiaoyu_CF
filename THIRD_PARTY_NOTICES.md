# Third-Party Notices

This repository vendors local copies of the libraries needed by `xiaoyu_CF.py`
so users do not need to install the unmaintained `cytoflow` package separately.

## Cytoflow

The `cytoflow/` directory is copied from the local working installation at:

`/Users/wuxiaoyu/miniforge3/envs/cytoflow/lib/python3.8/site-packages/cytoflow/`

Cytoflow is copyright Massachusetts Institute of Technology 2015-2018 and Brian
Teague 2018-2022. Its source headers state that it is distributed under the GNU
General Public License, version 2 or later.

Local packaging compatibility changes in this repository:

- `cytoflow/__init__.py` tolerates Matplotlib versions without
  `matplotlib.text._log`.
- The compiled scale extension and its Python wrapper were removed so the
  package remains pure Python. Linear and log scales remain available.
- The vendored Cytoflow tree was reduced to the operations and views used by
  `xiaoyu_CF.py` and the included example workflow.

## fcsparser

The `fcsparser/` directory is copied from the local working installation at:

`/Users/wuxiaoyu/miniforge3/envs/cytoflow/lib/python3.8/site-packages/fcsparser/fcsparser/`

The upstream package metadata identifies `fcsparser` as MIT licensed.
