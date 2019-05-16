from pyMarmoteMDP import *


def valueIteration(listAction,listCout,maxIter = 500,epsilon = 0.00001):
	critere="min"
	dim_SS = len(listAction[0])
	dim_SA = len(listAction)
	actionSpace = marmoteInterval(0,dim_SA-1)
	stateSpace = marmoteInterval(0,dim_SS-1)
	trans=sparseMatrixVector(dim_SA)
	reward = sparseMatrix(dim_SS,dim_SA)
	for i in range(len(listAction)) : 
		pi  = sparseMatrix(dim_SS)
		for j in range(len(listAction[i])):
			b=True
			for y in range(len(listAction[i][j])):
				reward.addToEntry(j,i,listCout[i])
				if(listAction[i][j][y] !=0):
					pi.addToEntry(j,y,listAction[i][j][y])
					b=False
				else:
					pi.addToEntry(j,y,0)
				if(b):
					pi.addToEntry(j,j,1)

		trans[i]=pi
	mdp1 = totalRewardMDP(critere, stateSpace, actionSpace, trans, reward)
	optimum = mdp1.valueIteration(epsilon,maxIter)
	optimum.writeSolution()

	return optimum

def policyIterationModified(listAction,listCout,maxIter = 500,epsilon = 0.00001,maxInIter = 1000,delta = 0.001 ):
	
	critere="min"
	dim_SS = len(listAction[0])
	dim_SA = len(listAction)
	actionSpace = marmoteInterval(0,dim_SA-1)
	stateSpace = marmoteInterval(0,dim_SS-1)
	trans=sparseMatrixVector(dim_SA)
	reward = sparseMatrix(dim_SS,dim_SA)
	for i in range(len(listAction)) : 
		pi  = sparseMatrix(dim_SS)
		for j in range(len(listAction[i])):
			for y in range(len(listAction[i][j])):
				reward.addToEntry(j,i,listCout[i])
				pi.addToEntry(j,y,listAction[i][j][y])			
		trans[i] = pi
	mdp1 = totalRewardMDP(critere, stateSpace, actionSpace, trans, reward)
	print("oui")
	optimum = mdp1.policyIterationModified(epsilon,maxIter,delta,maxInIter)
	optimum.writeSolution()


	return optimum
