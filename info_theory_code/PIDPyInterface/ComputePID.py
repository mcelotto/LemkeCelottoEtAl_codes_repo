#
#  This source code is part of:
#  NIT - Neuroscience Information Toolbox
#  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import sys
sys.path.append('$BROJA_2PID_PATH')
import numpy as np
import BROJA_2PID


class ComputeII:
    def __init__(self, p_src, n_S, n_R, n_C, n_singleton_dims):
        self.p_src = np.asarray(p_src)
        self.n_S = int(n_S)
        self.n_R = int(n_R)
        self.n_C = int(n_C)
        self.n_singleton_dims = int(n_singleton_dims)
        assert 0 <= self.n_singleton_dims <= 3, "Number of non-singleton dimensions must be between 0 and 3"

    def calculate(self):
        p_cr_s = dict()
        p_sr_c = dict()

        for c in range(0, self.n_C):
            for r in range(0, self.n_R):
                for s in range(0, self.n_S):
                    # matlab compacts singleton dimensions. handle this.
                    if self.n_singleton_dims == 0:
                        p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                        p_sr_c[(s,r,c)] = float(self.p_src[s][r][c])
                    elif self.n_singleton_dims == 1:
                        if self.n_C == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r])
                            p_sr_c[(s,r,c)] = float(self.p_src[s][r])
                        if self.n_R == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                            p_sr_c[(s,r,c)] = float(self.p_src[s][r][c])
                        if self.n_S == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                            p_sr_c[(s,r,c)] = float(self.p_src[s][r][c])
                    elif self.n_singleton_dims == 2:
                        if self.n_C != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                            p_sr_c[(s,r,c)] = float(self.p_src[s][r][c])
                        if self.n_R != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[r])
                            p_sr_c[(s,r,c)] = float(self.p_src[r])
                        if self.n_S != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s])
                            p_sr_c[(s,r,c)] = float(self.p_src[s])
                    else:
                        p_cr_s[(c,r,s)] = float(self.p_src)
                        p_sr_c[(s,r,c)] = float(self.p_src)

        pid_s = BROJA_2PID.pid(p_cr_s)
        pid_c = BROJA_2PID.pid(p_sr_c)
        return min(pid_s['SI'], pid_c['SI'])

class ComputePID:
    def __init__(self, p_src, n_S, n_R, n_C, n_singleton_dims):
        self.p_src = np.asarray(p_src)
        self.n_S = int(n_S)
        self.n_R = int(n_R)
        self.n_C = int(n_C)
        self.n_singleton_dims = int(n_singleton_dims)
        assert 0 <= self.n_singleton_dims <= 3, "Number of non-singleton dimensions must be between 0 and 3"

    def calculate(self):
        p_cr_s = dict()

        for c in range(0, self.n_C):
            for r in range(0, self.n_R):
                for s in range(0, self.n_S):
                    # matlab compacts singleton dimensions. handle this.
                    if self.n_singleton_dims == 0:
                        p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                    elif self.n_singleton_dims == 1:
                        if self.n_C == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r])
                        if self.n_R == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                        if self.n_S == 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                    elif self.n_singleton_dims == 2:
                        if self.n_C != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s][r][c])
                        if self.n_R != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[r])
                        if self.n_S != 1:
                            p_cr_s[(c,r,s)] = float(self.p_src[s])
                    else:
                        p_cr_s[(c,r,s)] = float(self.p_src)

        pid_s = BROJA_2PID.pid(p_cr_s)
        return [ pid_s['SI'], pid_s['UIY'], pid_s['UIZ'], pid_s['CI'] ]