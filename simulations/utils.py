"""
Utilities to create and analyze LAMMPS simulations from a "colloid" prepared by C++ code.
"""
import json
import scipy as sp
import numpy as np
from scipy.optimize import minimize
import copy

class Colloid:
	"""
	Colloid created by C++ code.
	"""
	def __init__(self):
		self.motif_coords = []
		self.motif_types = []
		self.boundary_coords = []
		self.boundary_types = []
		self.chirality = None
	
	def load(self, filename, chirality=True, sigma_b=0.1, sigma_m=1.0):
		"""
		Load colloid from a JSON file.

		This converts types in the JSON file into integers for LAMMPS
		to use and stores the mapping internally. Typically, the motif
		is specified with letter and boundary points labeled with 
		integers, where 0 is reserved for "stop codons."

		Parameters
		----------
		filename : str
			Name of file, e.g., colloid.json.
		chirality : bool
			Default (arbitrary) chirality associated with this file.
		sigma_b : float
			Diameter of boundary points.
		sigma_m : float
			Diameter of motif points.
		"""
		try:
			data = json.load(open(filename, 'r'))
		except IOError:
			raise Exception("Unable to load "+filename)

		try:
			self.motif_coords = data['Motif']['coords']
			self.motif_types = data['Motif']['types']
			self.boundary_coords = data['Properties']['boundary_coords']
			self.boundary_types = data['Properties']['boundary_ids']
		except:
			raise Exception("Unable to read from "+filename)

		# Map the types to integers
		assert(set(self.motif_types).intersection(set(self.boundary_types)) == set()), "Motif's types overlap with the boundary types"
		self.forward = dict(enumerate(np.unique(self.boundary_types).tolist()+np.unique(self.motif_types).tolist(), start=1))

		self.reverse = {v:k for k,v in self.forward.items()}

		# Assign chirality
		self.chirality = chirality

		# Assume all boundary points have the same size
		self.sigma_b = {k:sigma_b for k in self.reverse.keys() if not (k in self.motif_types)}
		self.sigma_m = {k:sigma_m for k in self.reverse.keys() if (k in self.motif_types)}

	def scale(self, inplace=False):
		"""
		Scale coordinates so that the minimum motif pair distance is unity.

		If inplace, then internally resets the diameter (sigma) of the beads.
		to sigma_m = 1, and sigma_b to the minimum pairwise distance between
		boundary points after the scaling.
		"""
		coords_ = self.coords

		# Scale based on motif coordinates
		min_pdist = np.min(sp.spatial.distance.pdist(self.motif_coords))
		com = np.mean(self.coords, axis=0)

		# Center, scale, shift back
		bc = ((np.array(self.boundary_coords) - com)/min_pdist + com).tolist()
		mc = ((np.array(self.motif_coords) - com)/min_pdist + com).tolist()
		if inplace:
			self.boundary_coords = bc
			self.motif_coords = mc

			min_pdist = np.min(sp.spatial.distance.pdist(self.boundary_coords))
			self.sigma_b = {k:min_pdist for k in self.reverse.keys() if not (k in self.motif_types)}
			self.sigma_m = {k:1.0 for k in self.reverse.keys() if (k in self.motif_types)}

			return min_pdist
		else:
			return bc, mc

	def eps(self, type_1):
		"""
		Return the interaction energy associated with a given type."""
		return 1.0
	
	def sigma(self, type_1):
		"""Return the particle diameter associated with a given type."""
		t = self.file_type(type_1)
		if t in self.sigma_b.keys():
			return self.sigma_b[t]
		if t in self.sigma_m.keys():
			return self.sigma_m[t]

	def lammps_type(self, file_type):
		"""Return the integer assigned to LAMMPS for a given type."""
		return self.reverse[file_type]

	def file_type(self, lammps_type):
		"""Return the type from the JSON file associated with an integer LAMMPS type."""
		return self.forward[lammps_type]

	def involves_motif(self, type_1, type_2):
		"""Is one of these LAMMPS types a motif?"""
		if self.file_type(type_1) in self.motif_types:
			return True
		if self.file_type(type_2) in self.motif_types:
			return True
		return False 

	def is_stop_codon(self, type_1):
		"""Is one of these LAMMPS types a stop codon?"""
		if self.reverse[0] == type_1: # '0' is reserved in C++ code for stop codon
			return True
		else:
			return False

	@property
	def coords(self):
		return np.vstack((np.array(self.boundary_coords), np.array(self.motif_coords))).copy()

	@property
	def types(self):
		return np.array([self.lammps_type(x) for x in self.boundary_types]+[self.lammps_type(x) for x in self.motif_types])

	def save(self, filename):
		"""Save object to disk."""
		f = open(filename, 'wb')
		import pickle
		pickle.dump(self, f, protocol=4)

class Analysis:
	"""Tools to analyze simulations after they have been run."""

	@staticmethod
	def thermo(log_filename):
		"""
		Get the thermodynamic output from a simulation log file.

		Returns
		-------
		results : dict(str, ndarray)
			Dictionary of thermodynamic properties over time.
		"""
		f = open(log_filename, 'r')
		lines = f.readlines()
		f.close()

		start = []
		end = []
		keys = 'Step Atoms Temp Press PotEng KinEng E_pair TotEng Volume'
		for i,l in enumerate(lines):
			if keys in l:
				start.append(i+1)
			elif 'Loop time of' in l:
				end.append(i-1)

		# In case loop didn't finish, just use the chunks that did for simplicity
		length = np.min([len(start), len(end)])

		# The first one is the initialization of velocities, so skip
		total = None
		for i,j in zip(start[1:length], end[1:length]):
			chunk = np.array([x.strip().split() for x in lines[i:j+1]], dtype=np.float64)
			if total is None:
				total = chunk
			else: 
				total = np.vstack((total, chunk))

		results = {}
		for i, k in enumerate(keys.split()):
			results[k] = total[:,i]

		return results

	@staticmethod
	def performance(log_filename):
		"""Get the speed (timesteps per second) of the simulations."""
		f = open(log_filename, 'r')
		lines = f.readlines()
		f.close()

		steps_per_second = []
		for l in lines:
			if "Performance" in l:
				steps_per_second.append(float(l.strip().split()[-2]))

		return np.array(steps_per_second, dtype=np.float64)

	@staticmethod
	def sanity_checks(log_filename):
		"""Check the LAMMPS log.thermo file for any errors."""
		f = open(log_filename, 'r')
		lines = f.readlines()
		f.close()

		# 1. Check for dangerous builds
		for l in lines:
			if "Dangerous" in l:
				if int(l.strip().split('=')[-1]) > 0:
					raise Exception('Dangerous builds found in {}'.format(log_filename))

		# 2. Check max distance to body ownership
		for i,l in enumerate(lines):
			if "comm_cut" in l:
				ctr = i+1
				break
        
		for l in lines:
			if "max distance from body owner to body atom" in l:
				d = float(l.strip().split('=')[0])
				pseudo_diameter = float(lines[ctr].split()[4])
				if d > pseudo_diameter:
					raise Exception('Max distance from body owner to body atom ({}) exceeds pseudo_diameter ({})'.format(d, pseudo_diameter))

		f.close()
		return True

	@staticmethod
	def unwrap(coords, box):
		"""
		Unwrap atom coordinates.

		Parameters
		----------
		coords : ndarray
			Potentially wrapped coordinates.
		box : ndarray
			Box lengths.
		Returns
		-------
		unwrapped : ndarray
			Unwrapped coordinates, based on first atom.
		"""
		box, coords = np.array(box), np.array(coords)

		# Take first coordinate as reference
		return coords - box * np.round((coords - coords[0]) / box)

	@staticmethod
	def identify_chirality(template, coords):
		"""
		Check if the chirality of a nanoparticle is the same as the template.
		Coordinates must be specified in the same order in the template and
		test nanoparticle.  Coordinates must also be unwrapped to work correctly.

		Parameters
		----------
		template : ndarray
			Unwrapped coordinates of the template.
		coords : ndarray
			Unwrapped coordinates of the molecule/nanoparticle.

		Example
	  	-------
		>>> identify_chirality(unwrap(template, box), unwrap(coords, box))

		Returns
		-------
		x : bool
			If chiralities are the same or not.
		"""
		template, coords = np.array(template), np.array(coords)

		core_t = template[-1, :2]
		core_c = coords[-1, :2]

		# Take cross product and compare sign
		vec_t = np.empty((2, 2), dtype=np.float64)
		vec_c = np.empty((2, 2), dtype=np.float64)
		for i in range(2):  # Use first two particles to make vectors
			vec_t[i, :] = template[i, :2] - core_t
			vec_c[i, :] = coords[i, :2] - core_c

		# If same sign, same chirality
		return np.sign(np.cross(vec_t[0], vec_t[1])) == np.sign(
			np.cross(vec_c[0], vec_c[1])
		)

class LAMMPS:
	"""Tools to build LAMMPS simulations easily."""

	@staticmethod
	def force_shifted_tanh(r_min, r_cut, eps, sigma, kappa, bins=1000):
		"""
		Compute a force-shifted tanh() potential in LAMMPS RSQ layout.
		"""
		rvals = np.sqrt(np.linspace(np.min([r_min, r_cut])**2, r_cut**2, bins))
		def u_base(r):
			u_ = eps/2.0*(1.0-np.tanh(kappa*(r-sigma)))
			du_ = -eps*kappa/2.0*(1.0/np.cosh(kappa*(r-sigma))**2)
			return u_, du_
	    
		u_shift, du_shift = u_base(r_cut)
		u, du = u_base(rvals)
	    
		return rvals, u - (rvals - r_cut)*du_shift - u_shift, -(du - du_shift)

	@staticmethod
	def force_shifted_yukawa(r_min, r_cut, eps, kappa, bins=1000):
		"""
		Compute a force-shifted Yukawa-like potential in LAMMPS RSQ layout.
		"""
		rvals = np.sqrt(np.linspace(np.min([r_min, r_cut])**2, r_cut**2, bins))
		u = -eps*np.exp(-kappa*r_cut)*(np.exp(-kappa*(rvals-r_cut)) - 1.0)
		f = -eps*kappa*np.exp(-kappa*rvals)
		
		return rvals, u, f

	@staticmethod
	def force_shifted_lennard_jones(r_min, r_cut, eps, sigma, alpha, bins=1000):
		"""
		Compute a discretized force-shifted Lennard-Jones potential in LAMMPS RSQ layout.
		"""
		rvals = np.sqrt(np.linspace(np.min([r_min, r_cut])**2, r_cut**2, bins))
		def lj(r):
			u = (4.0 * eps * ((sigma / r) ** (2 * alpha) - (sigma / r) ** (alpha)))
			f = (4.0 * eps * alpha / r * (2.0 * (sigma / r) ** (2 * alpha) - (sigma / r) ** (alpha)))
			return u, f

		shift, force_shift = lj(r_cut)
		u, f = lj(rvals)

		energy = u - shift + (rvals - r_cut) * force_shift
		force = f - force_shift

		return rvals, energy, force

	@staticmethod
	def rigid_body_min_comm(colloid, box, eps=0.1):
		"""
		Compute the minimum communication cutoff needed for LAMMPS.
	    
		Rigid particles get assigned to the particle closest to geometric
		center, so interaction "range" needs to increase as a result.

		https://docs.lammps.org/comm_modify.html
		https://docs.lammps.org/fix_rigid.html#fix-rigid-npt-command
		"""
		# Particle that LAMMPS considers the "core"
		unwrapped = Analysis.unwrap(colloid.coords, box)
		com = np.mean(unwrapped, axis=0)
		d = np.sqrt(np.sum((unwrapped-com)**2, axis=1))
		core_pk = np.argmin(d)

		# Max distance from the "core particle" to anyt particle in the rigid body
		max_d = np.max(np.sqrt(np.sum((unwrapped-unwrapped[core_pk])**2, axis=1)))

		# Additional epsilon/fudge factor
		return max_d + eps

	@staticmethod
	def read_lammps_potential(filename, i, j):
		"""Read LAMMPS potential from file."""
		f = open(filename, 'r')
		lines = f.readlines()
    
		for ctr,l in enumerate(lines):
			if '{}_{}'.format(i,j) == l.strip():
				start = ctr
				break
            
		N = int(lines[ctr+1].strip().split()[1])
		r_u_f = [lines[idx].strip().split()[1:4] for idx in np.arange(ctr+3, ctr+3+N)]
    
            
		return np.array(r_u_f, dtype=np.float64)

	@staticmethod
	def tabulate_potentials(colloid, filename='potentials.lammps', alpha=6, kappa_sigma=1.0, bins=1000, style='fslj'):
		"""
		Given a colloid, tabulate the pair potentials between all pairs.

		If style == 'fslj':
		Interactions are force-shifted LJ between all pairs of points, but
		are set to WCA cutoff when boundary points are (1) unlike (i!=j) or
		(2) involve a stop codon. No interaction (U=0) exists with motif
		points. LB mixing rules are applied.
		
		If style == 'fstanh':
		Interactions are force-shifted LJ cutoff at 2.5*sigma (WCA) when
		involving a motif point and stop codons.  Boundary points which are 
		unlike (i!=j) have U=0 (no interaction) while boundary points which
		have the same identity interact favorably.  LB mixing rules are 
		applied.

		Results are stored in the file under headings X_Y for interactions
		between types X and Y.

		Parameters
		----------
		colloid : Colloid
			Colloid object used to infer the particle types in the system.
		filename : str
			Name of file to write the results to.
		alpha : float
			fslj exponent to use. Ignored otherwise.
		kappa_sigma : float
			Dimensionless inverse decay length (kappa*sigma) for fstanh to use. Ignored otherwise.
		bins : int
			Number of bins to discretize energy and force into.
		style : str
			Type of potential to use: {'fslj': Force-Shifted Lennard-Jones, 'fstanh': Force-Shifted Tanh()}

		Returns
		-------
		max_rcut : float, rcut_dict : dict(str, (float, float))
			Maximum cutoff distance for all potentials; dictionary of (r_min, r_cut) for each particle type pairs.
		"""

		max_rcut = 0.0
		rcut_dict = {}
		with open(filename, "w") as fn:
			fn.write("# Tabulated pair potentials\n\n")
			for i in range(1, len(colloid.reverse.values())+1):
				for j in range(i, len(colloid.reverse.values())+1):
					fn.write("{}\n".format(str(i) + "_" + str(j)))
				
					if style == 'fslj':
						"""Cutoff at 2.5*sigma, rmin = 0.7*sigma, LB mixing, etc."""
						factor = 1.0
						sigma = (colloid.sigma(i) + colloid.sigma(j))/2.0
						eps = np.sqrt(colloid.eps(i)*colloid.eps(j))
						r_min = 0.7*sigma
						if (colloid.involves_motif(i, j)):
							# No interaction with motif (U = 0)
							factor = 0.0
							r_cut = 2.5*sigma
						elif (colloid.is_stop_codon(i) or colloid.is_stop_codon(j)):
							# WCA interaction between any boundary point any stop codon
							r_cut = 2.0**(1.0/alpha) * sigma
						elif (i != j):
							# Boundary points (not stop codons) that are different have WCA
							r_cut = 2.0**(1.0/alpha) * sigma
						else:
							# Boundary points with identical label interact favorably
							r_cut = 2.5*sigma

						rvals, energy, force = LAMMPS.force_shifted_lennard_jones(r_min, r_cut, eps=eps, sigma=sigma, alpha=alpha, bins=bins)
						energy *= factor
						force *= factor
					elif style == 'fstanh':
						sigma = (colloid.sigma(i) + colloid.sigma(j))/2.0
						eps = np.sqrt(colloid.eps(i)*colloid.eps(j))
						r_cut = 2.5*sigma
						r_min = 1.0e-6 # LAMMPS does not allow R=0 in tables
						if (colloid.involves_motif(i, j) or colloid.is_stop_codon(i) or colloid.is_stop_codon(j)):
							# WCA interaction with motif and with any stop codon
							r_min = 0.5*sigma
							r_cut = 2.0**(1.0/alpha) * sigma
							rvals, energy, force = LAMMPS.force_shifted_lennard_jones(r_min=r_min, r_cut=r_cut, eps=eps, sigma=sigma, alpha=alpha, bins=bins)
						elif (i != j):
							#  Boundary points (not stop codons) that are different have no interaction
							rvals, _, _ = LAMMPS.force_shifted_tanh(r_min=r_min, r_cut=r_cut, eps=eps, sigma=sigma, kappa=kappa_sigma/sigma, bins=bins)
							energy = np.zeros_like(rvals)
							force = -np.zeros_like(rvals)
						else:
							# Boundary points with identical label interact favorably - eps must be negative for this potential to be attractive
							rvals, energy, force = LAMMPS.force_shifted_tanh(r_min=r_min, r_cut=r_cut, eps=-1.0*eps, sigma=sigma, kappa=kappa_sigma/sigma, bins=bins)
					else:
						raise ValueError("unrecognized style {}".format(style))

					max_rcut = np.max([max_rcut, r_cut])
					rcut_dict['{}_{}'.format(i, j)] = (r_min, r_cut)
					fn.write("N {} RSQ {} {}\n\n".format(bins, r_min, r_cut))
					for idx, (r_, u_, f_) in enumerate(zip(rvals, energy, force)):
						fn.write("{} {} {} {}\n".format(idx + 1, '%.12f'%(r_**1), '%.12f'%u_, '%.12f'%f_))
					fn.write("\n")	
		return max_rcut, rcut_dict

	@staticmethod
	def minimum_bounding_box(coords, sigma):
		"""
		Minimize the bounding box for a set of coordinates.
		This rotates the coordinates as a rigid body to find the smallest area
		when the bounding box is aligned with the x and y axes.

		Parameters
		----------
		coords : ndarray
			Array of particle coordinates.
		sigma : float
			Padding to add to the bounding box; typically the diameter of
			boundary particle.
		Returns
		-------
		xyz : ndarray, bbox : ndarray
			Rotated, centered coordinates around which the bounding box is given.
		"""
		coords = np.array(coords)

		def bbox(points, sigma):
			"""
			Compute the bounding box.
			[xmin, xmax]
			[ymin, ymax]
			"""
			a = np.zeros((2, 2), dtype=np.float64)
			a[:, 0] = np.min(points, axis=0) - sigma / 2.0
			a[:, 1] = np.max(points, axis=0) + sigma / 2.0
			return a

		def center(points):
			com = np.mean(points, axis=0)
			return points - com

		def area(box):
			return np.prod(box[:, 1] - box[:, 0])

		def rotate(theta, points):  # Counterclockwise
			M = np.array(
				[
					[np.cos(theta), -np.sin(theta)],
					[np.sin(theta), np.cos(theta)],
				]
			)
			return np.matmul(M, points.T).T

		def func(x, points):
			return area(bbox(rotate(x[0], center(points[0])), sigma))

		result = minimize(
			fun=func,
			x0=[0],
			args=[
				coords.copy(),
				],
			method="L-BFGS-B",
			bounds=[(0, np.pi)],
		)  # not 2*pi

		if result.success:
			theta = result.x[0]
			xyz = rotate(theta, center(coords.copy()))
			return xyz, bbox(xyz, sigma)
		else:
			return result.message

	@staticmethod
	def tile(colloid, box, spacing, n=[0, 0]):
		"""
		Tile coordinates in a 2D simulation cell to make an initial configuration.

		Spreads out particles as much as possible and alternates chirality when
		placing particles.

		Parameters
		----------
		colloid : Colloid
			Colloid to use as a template.
		box : array-like
			Array of simulation box size; (L_x, L_y), for example.
		spacing : float
			Padding to add to the bounding box; typically the diameter of a 
			boundary particle.
		n : array-like
			Number of units of each enantiomorph to put in the box. This first
			will keep the orientation provided, the second will be "flipped".

		Returns
		-------
		positions : ndarray, identities : ndarray, bbox : ndarray
			Array of particle coordinates and identities in the simulation box,
			and the bounding box that goes around a single colloid.
		"""

		box = np.array(box)
		coords = colloid.coords
		ids = colloid.types
		template, bbox = LAMMPS.minimum_bounding_box(coords, sigma=spacing)
		
		def flip(template, bbox):
			"""Flips a template to change chirality."""
			enan = template.copy()
			enan[:, 0] = -enan[:, 0]  # Swap x coordinate
			return enan
		
		# Shift to have left hand at (0,0)
		template1 = template.copy() - bbox[:, 0].T
		template2 = flip(template1, bbox) + np.array([bbox[0][1]-bbox[0][0], 0.0])

		bbox = (bbox.T - bbox[:, 0]).T
		delta = bbox[:, 1] - bbox[:, 0]
		N = np.array(np.floor(box / delta), dtype=np.int32)

		def choose(total, n, total_each):
			# By default, alternate.
			choice = total % 2

			# If we've already got enough, choose the other
			if total_each[choice] >= n[choice]:
				choice = (total + 1) % 2

			# If the other is also complete, raise an error
			assert total_each[choice] < n[choice], "Already added all particles."
			return choice

		total = 0
		total_each = [0, 0]
		length = len(template)
		ntot = np.sum(n)
		positions = np.empty((ntot * length, 2), dtype=np.float64)
		identities = np.empty(ntot * length, dtype=np.int32)
		spacing = np.prod(N) // ntot
		counter = 0
		for j in range(N[1]):
			for i in range(N[0]):
				if counter % spacing == 0:  # Spread out
					shift = delta * np.array([i, j])
					enan = choose(total, n, total_each)
					if enan == 0:
						temp = template1
					else:
						temp = template2

					positions[total * length : total * length + length, :] = (
						temp + shift
					)
					identities[total * length : total * length + length] = np.array(
						ids, dtype=np.int32
					)
					total += 1
					total_each[enan] += 1

					if total == ntot:
						return positions, identities, bbox
				else:
					pass
				counter += 1

		raise Exception("Unable to place all particles.")

	@staticmethod
	def create_initial_configuration(colloid, box, spacing, n, filename, unit_mass=False):
		"""
		Create an initial configuration on a lattice.

		Parameters
		----------
		colloid : Colloid
			Colloid to use as a template.
		box : array-like
			Array of simulation box size; (L_x, L_y), for example.
		spacing : float
			Padding to add to the bounding box; typically the diameter of a 
			boundary particle.
		n : array-like
			Number of units of each enantiomorph to put in the box. This first
			will keep the orientation provided, the second will be "flipped".
		filename : str
			Name of file to write configuration to.
		unit_mass : bool
			If true, the mass of the entire colloid is set to 1.02 by evenly 
			dividing the mass over all constituent "atoms"; otherwise all
			atoms have a mass of 1.0.

		Returns
		-------
		pseudo_diameter : float
			"Diameter" of a colloid ((over)estimated as diagonal of bounding box).
		"""

		length = len(colloid.coords)
		positions, identities, bbox_ = LAMMPS.tile(colloid, box, spacing, n)

		with open(filename, "w") as f:
			f.write("# LAMMPS configuration file\n\n")
			f.write(
				"0 {} xlo xhi\n0 {} ylo yhi\n-0.5 0.5 zlo zhi\n".format(
					box[0], box[1]
				)
			)
			f.write("\n")
			f.write("{} atoms\n".format(len(positions)))
			f.write("\n")
			f.write("{} atom types\n".format(len(np.unique(identities))))
			f.write("\n")

			f.write("\nAtoms\n\n")  # "Molecular"=atom-ID molecule-ID atom-type x y z
			assert (len(positions) % length == 0), "Size of colloid disagrees with positions provided"
			for i in range(len(positions)):
				f.write(
					str(i + 1)
					+ "\t"
					+ str(i // length + 1)
					+ "\t"
					+ str(identities[i])
					+ "\t"
				)
				f.write(
					str(positions[i][0]) + "\t" + str(positions[i][1]) + "\t0\n"
				)

			f.write("\nMasses\n\n")  
			m = 1.0
			if unit_mass: # Set mass of whole colloid to unity
				m /= len(colloid.coords)
			for i in range(1, len(np.unique(identities)) + 1):
				f.write(str(i) + "\t" + str(m) + "\n")

			f.write("\nVelocities\n\n")  # Zero initial velocity, will be set later
			for i in range(1, len(positions) + 1):
				f.write(str(i) + "\t0.0\t0.0\t0.0\n")

		pseudo_diameter = np.sqrt(np.sum((bbox_[:,1] - bbox_[:,0])**2))
		return pseudo_diameter

