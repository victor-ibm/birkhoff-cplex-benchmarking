# Minimum Birkhoff Decomposition - CPLEX Benchmarking

This repository contains a notebook to reproduce the CPLEX benchmarking for the Minimum Birkhoff Decomposition problem in: https://arxiv.org/pdf/2504.03832

### Setup

```
conda create -n my_env python=3.8
conda activate my_env
pip install .
```

Download and install CPLEX: https://www.ibm.com/products/ilog-cplex-optimization-studio

### Run notebook

Download the QOBLIB dataset from: https://git.zib.de/qopt/qoblib-quantum-optimization-benchmarking-library. Change `PATH_TO_QOBLIB` in [mbd_cplex.ipynb](mbd_cplex.ipynb) to your local copy of QOBLIB. 

Select the matrix size (`n`) and density (`type`). By default, the results are stored in the `logs` folder. 

### IBM Public Repository Disclosure

All content in these repositories including code has been provided by IBM under the associated open source software license and IBM is under no obligation to provide enhancements, updates, or support. IBM developers produced this code as an open source project (not as an IBM product), and IBM makes no assertions as to the level of quality nor security, and will not be maintaining this code going forward.
