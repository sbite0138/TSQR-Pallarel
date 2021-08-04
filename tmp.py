import numpy as np
from pprint import pprint
from copy import deepcopy
res = [
    [-0.86593792,  1.16065968,  1.12301562, ],
    [0.94169364, -1.32029803,  1.28081268, ],
    [0.80258140, -0.51878546, -0.49875779, ],
    [0.50753022,  0.11355116,  0.08305245, ],
    [0.05188427,  0.41321093, -0.10360333, ],
    [0.28370011,  0.12993623,  0.41804100, ],
    [0.40222766, -0.35296807,  0.33954812, ],
    [0.27418665, -0.05283901,  0.70895174, ],
    [0.10249297,  0.37543815, -0.20296889, ],
]
L = deepcopy(res)
U = deepcopy(res)
for i in range(len(res)):
    for j in range(len(res[0])):
        if (i == j):
            U[i][j] = 1.0
        elif (i > j):
            U[i][j] = 0.0
        else:
            L[i][j] = 0.0


L = np.matrix(L)
U = np.matrix(U[:len(res[0])])
# print(L)
# print(U)
print(L*U)
# Q = np.matrix(Q)
# pprint(Q.T[2:, 2:]*A[2:, 2:]*Q[2:, 2:])
# pprint((Q.T*A*Q)[2:, 2:])
