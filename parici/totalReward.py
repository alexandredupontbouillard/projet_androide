from pyMarmoteMDP import *
import copy


def totalReward(listAction,listCout,maxIter = 500,epsilon = 0.00001):
	critere="min"
	dim_SS = len(listAction[0])+1
	dim_SA = len(listAction)+1
	actionSpace = marmoteInterval(0,dim_SA-2)
	stateSpace = marmoteInterval(0,dim_SS-2)
	trans=sparseMatrixVector(dim_SA)
	reward = sparseMatrix(dim_SS,dim_SA)
	#for i in range(len(listAction)):
		
	
	for i in range(len(listAction)) :
		for j in range(len(listAction[i])):
			if(j != len(listAction[i])-1):
				reward.addToEntry(j,i,listCout[i])
			else:
				reward.addToEntry(j,i,0.0)
	print("reward : OK")
	


	for i in range(len(listAction)) : 
		trans[i] = sparseMatrix(dim_SS)
		print("action : " ,i)
		for j in range(len(listAction[i])):
			print("##")
			b = True
			#verifie si tous les coefficients sont nuls si c'est le cas on dira que p(s|a,s)=1
			for y in	range(len(listAction[i][j])):
				if(listAction[i][j][y] != 0):
					b = False
			if(b):
				for y in range(len(listAction[i][j])):
					if(y!=j):
						print(0.0)
						trans[i].addToEntry(j,y,0.0)
					else : 
						print(1.0)
						trans[i].addToEntry(j,y,1.0)
			else:
				for y in range(len(listAction[i][j])):
					
					print(listAction[i][j][y])
					trans[i].addToEntry(j,y,listAction[i][j][y])

	print("actions : OK")
	mdp1 = totalRewardMDP(critere, stateSpace, actionSpace, trans, reward)
	print("totalReward : OK")
	mdp1.writeMDP()
	return mdp1

