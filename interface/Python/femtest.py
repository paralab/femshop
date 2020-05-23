from state import *

#basic usage. Initialize with a dict of options. Everything left out will be initialized to default value
femstate = State({
    RequiredParameters.Dimensions : 2,
    RequiredParameters.Domain: Domain([], 0),
    #can also just use a magic string instead of enums
    'SolutionMethod': SolutionMethods.CG
    #RequiredParameters.SolutionMethod: SolutionMethods.CG,
})

test = 0
