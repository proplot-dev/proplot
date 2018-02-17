from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot
import numpy as np

# Create "colormap object" from list of colors
# Here we create a random colormap. You might also want to load the colormap array
# with numpy "np.loadtxt" or pandas "pd.read_table"
N = 10
array = np.random.rand(N,3) # list/tuple of length-3 lists/tuples === N by 3 numpy array
# array = np.loadtxt("cmap.txt", delimiter=",") # loading from file
cmap = LinearSegmentedColormap.from_list('mymap', array) # this is a "colormap object"

# If you have a "startup" file that you import every time, add this line to "register"
# the colormap. Then you can use the colormap with "cmap='mymap'"
pyplot.register_cmap(cmap=cmap)

# Quick visualization
# Note the "cmap" argument can always be a string OR the object itself
n = 5
fig, ax = pyplot.subplots()
ax.set_title("Randomly generated colormap")
mappable = ax.contourf(range(n), range(n), np.random.rand(n,n), cmap='mymap')
# mappable = ax.contourf(range(n), range(n), np.random.rand(n,n), cmap=cmap)
fig.colorbar(mappable)
pyplot.show()
