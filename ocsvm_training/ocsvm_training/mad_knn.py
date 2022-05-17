import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation

class MAD_KNN():
	def __init__(self, data, kernel, verbose = True):
		self.graphs = data
		self.kernel = kernel
		self.verbose = verbose

	def find_vals(self, dist, n_end):
		vals = np.mean(dist[:,n_end-1:-1], axis = 1)
		vals[np.isnan(vals)] = 0

		return np.sort(vals)

	def fit(self, n = 7, plots = True, find_n = False):

		if find_n:
			n = max(10, int(0.03*len(self.graphs)))

		if not hasattr(self, 'dist'):
			self.dist = self.kernel.fit_transform(self.graphs)
			self.sorted_dist = np.sort(self.dist)
		self.vals = self.find_vals(self.sorted_dist, -n)

		mad = median_abs_deviation(self.vals, scale= 'normal')
		median = np.median(self.vals)
		
		tmin = median - 3*mad
		tmax = median + 3*mad

		mask = ((self.vals < tmin) | (self.vals > tmax))

		self.nu = len(self.graphs[mask])/len(self.graphs)

		if plots:
			fig, ax = plt.subplots()
			ax.scatter(np.arange(len(self.vals)), self.vals, color = 'blue', label = 'Inlier')
			ax.scatter(np.arange(len(self.vals))[mask], self.vals[mask], color = 'red', label = 'Outlier')
			ax.legend()
			plt.show()