#! /usr/bin/env python

"""
Toy Metropolis-coupling example
@author: jembrown

Starting with a surface consisting of two uniform densities separated by a deep valley. The
width of the valley is bigger than the proposal window, but its density is > 0.

**Python 2.7.10**
"""

# Import dependencies
from scipy.stats import uniform
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import math
import random

# ----> ADJUST TARGET DISTRIBUTION HERE	<----

# Define the target distribution ("toy posterior")
# Needs to integrate to 1
peakOneBounds = (0,0.4)
peakOneHeight = 0.5
valleyHeight = 0.000001
peakTwoBounds = (0.9,1.0)
peakTwoHeight = (8.0-(valleyHeight*(peakTwoBounds[0]-peakOneBounds[1])))


# ----> FUNCTION AND CLASS DEFINITIONS <----

# Define probability distribution on parameters of interest
def posterior(theta):
	"""
	A function to return the 'posterior' density for a given value of theta.
	"""
	if theta >= peakOneBounds[0] and theta <= peakOneBounds[1]:		# Peak One
		return peakOneHeight
	elif theta >= peakTwoBounds[0] and theta <= peakTwoBounds[1]:	# Peak Two
		return peakTwoHeight
	elif theta >= peakOneBounds[1] and theta <= peakTwoBounds[0]:	# Valley
		return valleyHeight
	else:															# Outside range
		return 0.0


def drawTheta(thetaCurr):
	"""
	This function provides proposed values for theta, given a current value. It currently
	uses only a uniform distribution centered on the current value. The size of the window
	is specified below (propWindowSize).
	"""
	propWindowSize = 0.1
	newTheta = uniform.rvs(loc=thetaCurr-(propWindowSize/2.0),scale=propWindowSize) # min = loc, max = loc+scale
	return newTheta


class chain(object):
	"""
	A class to define instances of Metropolis-coupled Markov chains.
	"""

	def __init__(self,lam,num,name,theta):
		"""
		Initializes variables associated with a chain. The parameter value of interest is
		"theta". The strength of heating is adjusted with "lam".
		"""
		self.name = name
		self.theta = theta
		self.thetaPost = posterior(self.theta)
		self.thetaProp = drawTheta(self.theta)
		self.lam = lam
		self.num = num
		self.temp = 1.0 + (lam*(self.num-1))
		self.samples = []
		self.swapCount = 0
		self.reportState(debug=False)

		# Parameter for Wang-Landau
		# Estimate wight DISTRIBUTION
		self.f = 1
		self.f_min = 0.0001
		self.min = 0.0
		self.max = 1.0
		self.N = 100
		self.bin_length = (self.max-self.min)/self.N
		self.parameter_bins = np.arange(self.min, self.max, self.bin_length)
		self.weight_distribution = [0] * self.N
		chain.histogram = [0] * self.N
		self.flat_threshold = 0.95

	def reportState(self,debug=False):
		"""
		A function to report a chain's current status.
		"""
		print ""
		print "Chain Name: %s" % self.name
	 	print "Chain Number: %s" % self.num
		print "Chain Lambda: %s" % self.lam
		print "Chain Temperature: %s" % self.temp
		print "Chain Theta: %f" % self.theta
		print "Chain Probability Density: %f" % self.thetaPost
		print "Chain Swap Count: %d" % self.swapCount
		if (debug):		# Turn on debugging statements, if needed
			print "Chain %s Temp Test: %f" % (self.name,math.pow(0.9,(1.0/self.temp)))
		print ""

	def update(self):
		"""
		Updates a chain according to the usual Metropolis-Hastings rules.
		"""
		self.thetaPropPost = posterior(self.thetaProp)
		self.thetaPost = posterior(self.theta)
		if (self.thetaPropPost >= self.thetaPost): # If proposed value has posterior density > curr value
			self.theta = self.thetaProp
		elif (self.thetaPropPost < self.thetaPost): # If proposed value has posterior density < curr value
			self.r = math.pow((self.thetaPropPost/self.thetaPost),(1.0/self.temp))
			self.ranUnifDraw = uniform.rvs()
			if (self.ranUnifDraw <= self.r):
				self.theta = self.thetaProp
		else:
			print "Problem calculating proposal ratio for chain %s!" % self.name
			print (self.thetaProp,self.thetaPropPost)
			print (self.theta,self.thetaPost)
			print ""
		self.thetaProp = drawTheta(self.theta) # Finish by drawing new proposed value of theta
		self.thetaPropPost = posterior(self.thetaProp)	# Calculating densities for current and proposed thetas, in case
		self.thetaPost = posterior(self.theta)
		self.samples.append(self.theta)

	def WLupdate(self):
		"""
		Updates a chain according to the usual Metropolis-Hastings rules.
		"""
		self.thetaPropPost = posterior(self.thetaProp)
		self.thetaPost = posterior(self.theta)
		self.entropyProp = math.exp(self.weight_distribution[self.get_index(self.thetaProp)])
		self.entropy = math.exp(self.weight_distribution[self.get_index(self.theta)])
		acceptance_ratio = (self.thetaPropPost * self.entropy) / (self.thetaPost * self.entropyProp)
		self.ranUnifDraw = uniform.rvs()
		if (1 <= acceptance_ratio) or (self.ranUnifDraw <= acceptance_ratio): # If proposed value has posterior density > curr value
			self.theta = self.thetaProp
		elif (1 > acceptance_ratio): # If proposed value has posterior density < curr value
			pass
		else:
			print "Problem calculating proposal ratio for chain %s!" % self.name
			print (self.thetaProp,self.thetaPropPost)
			print (self.theta,self.thetaPost)
			print ""
		self.thetaProp = drawTheta(self.theta) # Finish by drawing new proposed value of theta
		self.thetaPropPost = posterior(self.thetaProp)	# Calculating densities for current and proposed thetas, in case
		self.thetaPost = posterior(self.theta)
		self.samples.append(self.theta)

	def get_index(self, theta):
		index = 0
		while (theta > self.parameter_bins[index + 1]):
			index = index + 1
			if index is (self.N - 1):
				break
		return index



	def WLsample(self):
		"""
		Updates a chain according to the usual Metropolis-Hastings rules.
		"""
		self.thetaPropPost = posterior(self.thetaProp)
		self.thetaPost = posterior(self.theta)
		self.entropyProp = math.exp(self.weight_distribution[self.get_index(self.thetaProp)])
		self.entropy = math.exp(self.weight_distribution[self.get_index(self.theta)])
		acceptance_ratio = (self.thetaPropPost * self.entropy) / (self.thetaPost * self.entropyProp)
		self.ranUnifDraw = uniform.rvs()
		if (1 <= acceptance_ratio) or (self.ranUnifDraw <= acceptance_ratio): # If proposed value has posterior density > curr value
			self.theta = self.thetaProp
		elif (1 > acceptance_ratio): # If proposed value has posterior density < curr value
			pass
		else:
			print "Problem calculating proposal ratio for chain %s!" % self.name
			print (self.thetaProp,self.thetaPropPost)
			print (self.theta,self.thetaPost)
			print ""
		index = self.get_index(self.theta)
		self.weight_distribution[index] = self.weight_distribution[index] + self.f
		self.histogram[index] = self.histogram[index] + 1
		self.thetaProp = drawTheta(self.theta) # Finish by drawing new proposed value of theta
		self.thetaPropPost = posterior(self.thetaProp)	# Calculating densities for current and proposed thetas, in case
		self.thetaPost = posterior(self.theta)
		self.samples.append(self.theta)

	def flat(self):
		mean = np.average(np.array(self.histogram))
		index = 0
		for h in self.histogram:
			if(h < mean * self.flat_threshold) or (h is 0):
				print("mean: " + str(mean))
				print("Index: " + str(index))
				print("parameter: " + str(self.parameter_bins[index]))
				print("weight: " + str(self.weight_distribution[index]))
				print("f: " + str(self.f))
				return False
			index = index + 1

		return True

# ----> DEFINE CHAIN CONDITIONS AND RUN <----

# Set chain length and run it
ngens = 2000000				# Total length of chain
sampleFreq = 100			# Change this to something >1 if you want to space out samples.
updateFreq = ngens*0.1		# Frequency of screen updates to make sure chain is running.
samples = []				# Vector to hold sampled values from cold chain

# Initiate chains
print ""
print "Starting chain information: "
chainLambda = 0
chain 		= chain(lam=chainLambda,num=1,name="One",theta=uniform.rvs(loc=0,scale=1))	# Cold chain

while chain.f > chain.f_min:
	chain.WLsample()
	if chain.flat() is True:
		chain.f = chain.f * 0.5
		chain.histogram = [0] * chain.N

# Run chains
for gen in range(ngens):
	#if (uniform.rvs() <= 0.5):		# Performs typical within-chain M-H update with 50% probability
	# Running Metropolis-Hastings sampling within each chain
	chain.WLupdate()

	if (gen % sampleFreq == 0):				# Appending latest sample from cold chain
		samples.append(chain.theta)
	if (gen % updateFreq == 0):				# Reporting update to screen
		print "Generation %s" % gen

# Report chain states at end of run
chain.reportState()

# Summarizing MCMC samples

# Numerical summaries
burnin = int((ngens/sampleFreq)*0.1)	# Using a 10% burn-in
postBurnSamples = samples[burnin+1:]

entropies = [math.exp(entropy) for entropy in chain.weight_distribution]
weight_mean = np.mean(chain.weight_distribution)
weight_sum = np.sum(entropies)

# sift to histogram
corrected_histogram = [0] * chain.N
for samp in postBurnSamples:
	corrected_histogram[chain.get_index(samp)] += 1
index = 0
figure_weight_distribution = [0] * chain.N
for bin in chain.parameter_bins:
	corrected_histogram[index] = corrected_histogram[index] * math.exp(chain.weight_distribution[index] - weight_mean) / weight_sum
	figure_weight_distribution[index] = 1 / (math.exp(chain.weight_distribution[index] - weight_mean) / weight_sum)
	index = index + 1
index = 0
sum_ch = np.sum(corrected_histogram)
for c in corrected_histogram:
	corrected_histogram[index] = corrected_histogram[index]/sum_ch
	index = index + 1

index = 0
peakOne = 0
peakTwo = 0
for bin in chain.parameter_bins:
	if(bin < peakOneBounds[1]):
		peakOne = peakOne + corrected_histogram[index]
	if(bin > peakTwoBounds[0]):
		peakTwo = peakTwo + corrected_histogram[index]
	index = index + 1

print("Posterior Mean: %f" % np.mean(postBurnSamples))
postBurnSamples.sort()  # post-burnin samples will be sorted after this is called
print("Posterior 95% credibility interval: "+"(%f,%f)" % (postBurnSamples[int(len(postBurnSamples)*0.025)],postBurnSamples[int(len(postBurnSamples)*0.975)]))
print("Estimated probability of peak one: %f" % (float(sum(1 for samp in postBurnSamples if samp <= peakOneBounds[1]))/float(len(postBurnSamples))))
print("Estimated probability of peak two: %f" % (float(sum(1 for samp in postBurnSamples if samp >= peakTwoBounds[0]))/float(len(postBurnSamples))))

print("Estimated probability of peak one: %f" % peakOne)
print("Estimated probability of peak two: %f" % peakTwo)

# Marginal histogram
pl.figure()
pl.xlim(0,1)
pl.hist(postBurnSamples, bins=50, color="gray")
pl.show()

plt.figure()

plt.plot(chain.parameter_bins, corrected_histogram)
plt.show()

plt.figure()
plt.plot(chain.parameter_bins, figure_weight_distribution)
plt.show()
# Trace plot
plt.figure()
plt.plot(range(ngens/sampleFreq),samples)
plt.ylim(0,1)
plt.ylabel("Theta")
plt.xlabel("Sample Number")
plt.show()
