import numpy as np
a = [
    [0.315598, 0.284943, ],
    [0.240601, 0.484127, ],
]
R = [
    [-0.396851, -0.520117, ],
    [0.000000, 0.212250, ],
]
Y = [
    [1.000000, 0.000000, ],
    [0.337710, 1.000000, ],
]
T = [
    [1.795255, 0.000000, ],
    [0.000000, 0.000000, ],
]


a = np.matrix(a)
R = np.matrix(R)
Y = np.matrix(Y)
T = np.matrix(T)
tmp = -Y*T*Y.T
for i in range(len(tmp)):
    tmp[i, i] += 1.0
tmp = tmp*R
tmp = tmp-a
print("error:", np.linalg.norm(tmp, 2)/np.linalg.norm(a, 2))
