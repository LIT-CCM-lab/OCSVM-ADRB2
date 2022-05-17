import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator

class QMS():
	def __init__(self, data, kernel, verbose = True):
		self.graphs = data
		self.kernel = kernel
		self.verbose = verbose

	def _find_vals(self, dist, n_end):
		
		return np.mean(dist[:,n_end-1:-1], axis = 1)


	def _find_knee(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, eta = 0, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):

		if find_n:
			n = max(10, int(0.03*len(self.graphs)))

		if not hasattr(self, 'dist'):
			self.dist = self.kernel.fit_transform(self.graphs)
			self.sorted_dist = np.sort(self.dist)
			self.vals = self._find_vals(self.sorted_dist, -n)
			self.vals[np.isnan(self.vals)] = 0
		
		knee_finder = KneeLocator(np.arange(len(self.vals)), np.sort(self.vals), curve = curve, direction = direction, S = S_knee, interp_method = interp_method, polynomial_degree= polynomial_degree)
		if knee_finder.knee_y is None:
			raise ValueError('The given dataset does not present a knee/elbow with the desired characteristics')
		else:
			self.knee = knee_finder.knee_y
		if plots:
			knee_finder.plot_knee_normalized()
			plt.show()
			fig, ax = plt.subplots()
			ax.scatter(np.arange(len(self.vals)), np.sort(self.vals), color = 'blue')
			#ax.set_ylim(bottom = 0, top = 1)
			ax.vlines(knee_finder.knee, ymin = np.min(self.vals), ymax = np.max(self.vals), linestyles = 'dashed')
			ax.hlines(self.knee, xmin = -1, xmax = len(self.vals), linestyles = 'dashed')
			plt.show()


class QMS2(QMS):
	def __init__(self, data, kernel, verbose = True):
		super().__init__(data, kernel, verbose)
		self.graphs = data
		self.kernel = kernel
		self.verbose = verbose


	def fit(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, eta = 0, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):

		self._find_knee(n = n, curve = curve, direction = direction, S_knee = S_knee, plots = plots, find_n = find_n, interp_method = interp_method, polynomial_degree = polynomial_degree)

		self.mask = self.vals > self.knee

		return self.graphs[self.mask]

class QMS1(QMS):
	def __init__(self, data, kernel, verbose = True):
		super().__init__(data, kernel, verbose)
		self.graphs = data
		self.kernel = kernel
		self.verbose = verbose


	def fit(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, eta = 0, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):

		self._find_knee(n = n, curve = curve, direction = direction, S_knee = S_knee, plots = plots, find_n = find_n, interp_method = interp_method, polynomial_degree = polynomial_degree)

		self.mask = self.vals > self.knee 

		self.fitted_nu = 1-((len(self.graphs[self.mask])*(1- eta))/len(self.graphs))