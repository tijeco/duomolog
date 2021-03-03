class Duo:
	def __init__(self, leftName, left, rightName, right):
		self.left_union_rightName = leftName + "_union_" + rightName
		self.left_not_rightName = leftName + "_not_" + rightName
		self.left_intersect_rightName = leftName + "_intersect_" + rightName
		self.right_not_leftName = rightName + "_not_" + leftName

		self.subsets = {
			leftName: left,
			rightName: right,
			self.left_union_rightName: left | right,
			self.left_not_rightName: left - right,
			self.left_intersect_rightName: left & right,
			self.right_not_leftName: right - left
		}

	def dropEmpty(self):
		emptySubsets = []
		for subset in self.subsets:
			if self.subsets[subset] == set():
				print(subset, "is empty, removing it")
				emptySubsets.append(subset)
		for emptySubset in emptySubsets: 
			del self.subsets[emptySubset]

	def dropRedundant(self):
		nrSubSets = []
		redundantSubsets = []
		for subset in self.subsets:
			if self.subsets[subset] not in nrSubSets:
				nrSubSets.append(self.subsets[subset])
			else:
				print(subset, "is redundant, removing it")
				redundantSubsets.append(subset)
		for redundantSubset in redundantSubsets:
			del self.subsets[redundantSubset]
	def intersectOnly(self):
		self.subsets = {self.left_intersect_rightName:self.subsets[self.left_intersect_rightName]}
	
def test():
	s = {
		1:"MENDEL_01",
		2:"MENDEL_02",
		3:"MENDEL_03",
		4:"MENDEL_04",
		6:"MENDEL_06",
		8:"MENDEL_08"
		}
	lrDuo = Duo("blast",{1,2,3},"hmm",{4,6,8})
	print(lrDuo.subsets)
	# lrDuo.dropEmpty()
	# lrDuo.dropRedundant()
	lrDuo.intersectOnly()
	print(lrDuo.subsets)
	lrDuo.dropEmpty()
	print(lrDuo.subsets)

	
# test()

