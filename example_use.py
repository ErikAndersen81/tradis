from tradis import tradis
import numpy as np

# Create two arrays each consisting of two spatial
# and one temporal dimension, i.e. the outer
# axis of an array contains 2D points (x,y) at
# time t, giving [..., [x_i,y_i,t_i], ...] where
# t_i <= t_(i+1).
p = np.random.rand(15, 3)
q = np.random.rand(12, 3)
# sort by t
p = p[p[:, 2].argsort()]
q = q[q[:, 2].argsort()]

# call tradis
print(tradis(p, q))
