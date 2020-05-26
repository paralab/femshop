from enum import Enum, unique

class State:
	def __init__(self, configDict):
		self.configDict = {}
		for param in RequiredParameters:
			if param in configDict:
				self.setParameter(param, configDict[param])
				del configDict[param]
			else:
				self.setParameter(param)

		#merge in options with string keys
		self.configDict = {**self.configDict, **configDict}

	def setParameter(self, id, value=None):
		key = id
		if isinstance(id, RequiredParameters):
			key = id.name
			if value is None:
				#use default
				value = id.value

		self.configDict[key] = value

	def getParameter(self, id):
		key = id
		if isinstance(id, RequiredParameters):
			key = id.name

		return self.configDict[key]

class Domain:
	#configure boundaries, boundary conditions, geometry, mesh
	def __init__(self, boundaries, geometry, ):
		self.boundaries = boundaries
		self.geometry = geometry
		#...and so on...

@unique
class SolutionMethods(Enum):
	CG = 1
	DG = 2
	HDG = 3

#enum value is the default if unset by user.
#can be simple value, enum, or full class.
class RequiredParameters(Enum):
	Dimensions = 2
	SolutionMethod = SolutionMethods.CG
	Domain = Domain([],0)
	#... and many more ...  