import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator

class QMS():
	'''
	Base class for the implementation of Quick Model Selection (QMS) methods for parameters estimation for OCSVM models.
	This method can be used with any python module for graph kernel estimation as long as the kernel have a fit_transform method returning a square numpy array.
	Similairty is defined as the average similairty to the KNN points.
	:param data: Input graphs
	:type data: list of graphs
	:param kernel: The kernel used to compute the graph similairty
	:type kernel: graph kernel
	'''
	def __init__(self, data, kernel):
		self.graphs = data
		self.kernel = kernel

	def _find_vals(self, dist, n_end):
		
		return np.mean(dist[:,n_end-1:-1], axis = 1)


	def _find_knee(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):
		'''
		Find the knee for the pairwise similairty distribution.
		The method can be called multiple times changing the parameters, without the need to recompute the similarity matrix.
		Further information on the parameters used for knee identification can be found in the kneed module documentations:
		https://kneed.readthedocs.io/en/stable/
		:param n: Number of nearest neighbours for the analysis
		:type n: int
		:param curve: kneed parameter
		:type curve: string
		:param direction: kneed parameter
		:type direction: string
		:param S_knee: Sensibility of the kneed algorithm
		:type S_knee: int
		:param plots: show the plots of the obtained distribution
		:type plots: bool
		:param find_n: Defines K of KNN as the value corresponding to 3% of the data
		:type find_n: bool
		:param interp_method: kneed parameter
		:type interp_method: string
		:param polynomial_degree: kneed parameter
		:tpye polynomial_degree: int
		'''

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
			knee_finder.plot_knee()
			plt.show()
			knee_finder.plot_knee_normalized()
			plt.show()
			


class QMS2(QMS):
	'''
	Implementation of the QMS2 method
	'''
	def __init__(self, data, kernel):
		super().__init__(data, kernel)
		self.graphs = data
		self.kernel = kernel


	def fit(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):
		'''
		Select the graphs to be used for model training.
		:param n: Number of nearest neighbours for the analysis
		:type n: int
		:param curve: kneed parameter
		:type curve: string
		:param direction: kneed parameter
		:type direction: string
		:param S_knee: Sensibility of the kneed algorithm
		:type S_knee: int
		:param plots: show the plots of the obtained distribution
		:type plots: bool
		:param find_n: Defines K of KNN as the value corresponding to 3% of the data
		:type find_n: bool
		:param interp_method: kneed parameter
		:type interp_method: string
		:param polynomial_degree: kneed parameter
		:tpye polynomial_degree: int
		:returns: graphs to be used for model training
		:rtype: list of graphs
		'''

		self._find_knee(n = n, curve = curve, direction = direction, S_knee = S_knee, plots = plots, find_n = find_n, interp_method = interp_method, polynomial_degree = polynomial_degree)

		self.mask = self.vals > self.knee

		return self.graphs[self.mask]

class QMS1(QMS):
	'''
	Implementation of the QMS1 method
	:param n: Number of nearest neighbours for the analysis
	:type n: int
	:param curve: kneed parameter
	:type curve: string
	:param direction: kneed parameter
	:type direction: string
	:param S_knee: Sensibility of the kneed algorithm
	:type S_knee: int
	:param eta: The eta value of QMS1, the values must be in [0,1)
	:type eta: float
	:param plots: show the plots of the obtained distribution
	:type plots: bool
	:param find_n: Defines K of KNN as the value corresponding to 3% of the data
	:type find_n: bool
	:param interp_method: kneed parameter
	:type interp_method: string
	:param polynomial_degree: kneed parameter
	:tpye polynomial_degree: int
	'''
	def __init__(self, data, kernel):
		super().__init__(data, kernel)
		self.graphs = data
		self.kernel = kernel


	def fit(self, n = 7, curve = 'convex', direction = 'decreasing', S_knee = 1, eta = 0, plots = True, find_n = False, interp_method = 'interp1d', polynomial_degree = 7):

		self._find_knee(n = n, curve = curve, direction = direction, S_knee = S_knee, plots = plots, find_n = find_n, interp_method = interp_method, polynomial_degree = polynomial_degree)

		self.mask = self.vals > self.knee 

		self.fitted_nu = 1-((len(self.graphs[self.mask])*(1- eta))/len(self.graphs))