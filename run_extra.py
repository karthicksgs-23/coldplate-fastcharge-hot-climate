# -*- coding: utf-8 -*-
"""Fine supply-temperature sweep (1 K steps) for Fig. 4: C-rate 2 and 3, TIM 1/3/6 W/(m K), reference cold plate."""
import json
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from run_study import _run

RES = Path(__file__).parent / "results"

if __name__ == "__main__":
    Ts = np.arange(20, 45, 1.0)
    cases = [dict(C_rate=c, Q_lpm=2.0, layout="upass", tim_k=k, T_supply=float(t))
             for c, k, t in product([2, 3], [1.0, 3.0, 6.0], Ts)]
    with Pool(8) as pool:
        rs = pool.map(_run, cases)
    json.dump(rs, open(RES / "ts_sweep.json", "w"), indent=1, default=float)
    print("done", len(rs))
