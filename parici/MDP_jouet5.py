from pyMarmoteMDP import * 

critere = "min"
epsilon = 0.00001
maxIter = 500
dim_SS = 5
dim_SA = 4
actionSpace =marmoteInterval(0,2)
stateSpace = marmoteInterval(0,3)

trans=sparseMatrixVector(dim_SA)
reward = sparseMatrix(dim_SS,dim_SA)
P0 = sparseMatrix(dim_SS)
P0.addToEntry(0,0,0)
P0.addToEntry(0,1,0.875)
P0.addToEntry(0,2,0.0625)
P0.addToEntry(0,3,0.0625)
P0.addToEntry(1,0,0)
P0.addToEntry(1,1,0.75)
P0.addToEntry(1,2,0.125)
P0.addToEntry(1,3,0.125)
P0.addToEntry(2,0,0)
P0.addToEntry(2,1,0)
P0.addToEntry(2,2,0.5)
P0.addToEntry(2,3,0.5)
P0.addToEntry(3,0,0)
P0.addToEntry(3,1,0)
P0.addToEntry(3,2,0)
P0.addToEntry(3,3,1.0)
trans[0] = P0
P1 =sparseMatrix(dim_SS)
P1.addToEntry(0,0,0)
P1.addToEntry(0,1,0.875)
P1.addToEntry(0,2,0.0625)
P1.addToEntry(0,3,0.0625)
P1.addToEntry(1,0,0)
P1.addToEntry(1,1,0.75)
P1.addToEntry(1,2,0.125)
P1.addToEntry(1,3,0.125)
P1.addToEntry(2,0,0)
P1.addToEntry(2,1,1.0)
P1.addToEntry(2,2,0)
P1.addToEntry(2,3,0)
P1.addToEntry(3,0,0)
P1.addToEntry(3,1,0)
P1.addToEntry(3,2,0)
P1.addToEntry(3,3,1.0)
trans[1] = P1
P2 = sparseMatrix(dim_SS)
P2.addToEntry(0,0,0)
P2.addToEntry(0,1,0.875)
P2.addToEntry(0,2,0.0625)
P2.addToEntry(0,3,0.0625)
P2.addToEntry(1,0,0)
P2.addToEntry(1,1,0)
P2.addToEntry(1,2,0)
P2.addToEntry(1,3,0)
P2.addToEntry(2,0,0)
P2.addToEntry(2,1,0)
P2.addToEntry(2,2,0)
P2.addToEntry(2,3,0)
P2.addToEntry(3,0,0)
P2.addToEntry(3,1,0)
P2.addToEntry(3,2,0)
P2.addToEntry(3,3,0)
trans[2] = P2

Reward  = sparseMatrix(dim_SS, dim_SA);
Reward.addToEntry(0,0,0);
Reward.addToEntry(0,1,4000)
Reward.addToEntry(0,2,6000)
Reward.addToEntry(1,0,1000)
Reward.addToEntry(1,1,4000)
Reward.addToEntry(1,2,6000)
Reward.addToEntry(2,0,3000)
Reward.addToEntry(2,1,4000)
Reward.addToEntry(2,2,6000)
Reward.addToEntry(3,0,3000)
Reward.addToEntry(3,1,4000)
Reward.addToEntry(3,2,6000)


#stateSpace.enumerate()
#actionSpace.enumerate()

mdp1 =totalRewardMDP(critere, stateSpace, actionSpace, trans, Reward)
mdp1.writeMDP()

optimum = mdp1.valueIteration(epsilon, maxIter)
optimum.writeSolution()
