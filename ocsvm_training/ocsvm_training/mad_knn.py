import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation

class MAD_KNN():
	'''
	Train an OCSVM model using the MAD-KNN method.
	The method estimates the fraction of outliers in the dataset based on the median absolute deviation (MAD) based on the KNN-average similarity.
	This method can be used with any python module for graph kernel estimation as long as the kernel have a fit_transform method returning a square numpy array.
	:param data: Input graphs
	:type data: list of graphs
	:param kernel: The kernel used to compute the graph similairty
	:type kernel: graph kernel
	'''
	def __init__(self, data, kernel):
		self.graphs = data
		self.kernel = kernel

	def find_vals(self, dist, n_end):
		vals = np.mean(dist[:,n_end-1:-1], axis = 1)
		vals[np.isnan(vals)] = 0

		return np.sort(vals)

	def fit(self, n = 7, plots = True, find_n = False):
		'''
		Find the nu value for the given dataset.
		The method can be called multiple times changing the parameters, without the need to recompute the similarity matrix.
		Multiple paramters can be tested.
		:param n: Number of nearest neighbours for the analysis
		:type n: int
		:param plots: show the plots of the obtained distribution
		:type plots: bool
		:param find_n: Defines K of KNN as the value corresponding to 3% of the data
		:type find_n: bool
		'''
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
			ax.set_xlabel('Rank')
			ax.set_ylabel('KNN-average similarity')
			plt.show()