import numpy as np

alpha = 1 # thermal diffusitivity
west = 0 # domains limits
east = 1
c = 1 / (3*alpha)

k = 2 # operators order of accuracy
m = 2*k+1 # min number of cells to attain the desired accuracy (why +1?)
dx = (east - west)/m

t = 1
dt = c * dx**2 # von neumann stability criteriot for explicit scheme 

L = None # let's define Laplacian

# |> IC --------------

U = np.zeros((m+2, 1))

# |> BC ----------------

# hold 100 at each end
U[0] = 100
U[-1] = 100

# |>  Grid -------------







