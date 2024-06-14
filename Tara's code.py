import numpy as np
import itertools
import math


import pandas as pd
import scipy.stats

def ent(data):
    """Calculates entropy of the passed `pd.Series`
    """
    p_data = data.value_counts()           # counts occurrence of each value
    entropy = scipy.stats.entropy(p_data)  # get entropy from counts
    return entropy


d=2

range_of_K=8

AvgNumTests2Solve=np.empty(range_of_K, float)
SDEV_NumTests2Solve=np.empty(range_of_K, float)
max_NumTests2Solve=np.empty(range_of_K, float)
min_NumTests2Solve=np.empty(range_of_K, float)
print(np.shape(AvgNumTests2Solve))

#K=16


for K in range(24, 26):
    print("hi")
    AllSolutions = np.transpose(np.array(list(itertools.product(np.arange(d), repeat=K))))
    (dim, num_All_soln) = np.shape(AllSolutions)

    print(num_All_soln)

    round_up_K = pow(2, math.ceil(np.log2(K)))
    # print(round_up_K)

    for j in range(K, round_up_K):
        AllSolutions = np.concatenate((AllSolutions, np.zeros((1, num_All_soln))), axis=0)

    (dim, num_All_soln) = np.shape(AllSolutions)

    from scipy.linalg import hadamard

    H = hadamard(round_up_K, dtype=int)
    meas_mat = 1 / 2 * (np.ones(np.shape(H)) + H)

    NumTests2Solve = np.empty(50, int)

    for trial in np.arange(50):

        asoln = AllSolutions[:, np.random.choice(num_All_soln, size=1)]
        # print(asoln)
        Observations = np.matmul(meas_mat, asoln)
        # print(Observations)
        remaining_solutions = AllSolutions
        (dim, num_remaining_soln) = np.shape(remaining_solutions)
        # print(num_remaining_soln)

        # maxnumtests= np.arange(K)
        for num_tests in np.arange(round_up_K):

            E = np.zeros((round_up_K, 1))
            PossibleObservations = np.array(np.matmul(meas_mat, remaining_solutions))
            # print(PossibleObservations[1,:])
            for i in range(round_up_K):
                z = PossibleObservations[i:i + 1, :]
                (n, m) = np.shape(z)
                y = np.reshape(z, (m,))
                data = pd.Series(y)
                E[i] = ent(data)
            # print(E)

            Selected_meas_index = np.argmax(E)
            # print(Observations[Selected_meas_index])
            # print(PossibleObservations[Selected_meas_index,:])
            x = np.where(PossibleObservations[Selected_meas_index, :] == Observations[Selected_meas_index])

            (n, num_rem_solns) = np.shape(x)
            x2 = np.reshape(x, (num_rem_solns,))

            all_rows = np.arange(dim)
            the_solutions = remaining_solutions[np.ix_(all_rows, x2)]
            if num_rem_solns == 1:
                NumTests2Solve[trial] = num_tests
                # print("found the solution after " + str(num_tests) + " tests") # and it is: " + str(the_solutions))

                break
            else:
                remaining_solutions = the_solutions
                # print(str(num_rem_solns+1) + " possible solutions left to go")

    AvgNumTests2Solve[K] = np.mean(NumTests2Solve)
    SDEV_NumTests2Solve[K] = np.std(NumTests2Solve)
    # max_NumTests2Solve[K] = np.max(NumTests2Solve)
    # min_NumTests2Solve[K] = np.min(NumTests2Solve)

    print("found the solution, when K=" + str(K) + ", after " + str(AvgNumTests2Solve[K]) +
          " tests on average; compared to 2*K*np.log2(d+1)/np.log2(K+1)=" + str(
        2 * K * np.log2(d + 1) / np.log2(K + 1)))

# plot(AvgNumTests2Solve)


X = np.arange(2, 26)
Y = AvgNumTests2Solve[2:26]
STD = SDEV_NumTests2Solve[2:26]

import matplotlib.pyplot as plt

#


z = 2 * X * np.log2(d + 1) / np.log2(X + 1)
# plt.plot(X, Y, color='r', label='Avg # of Tests Required')
plt.errorbar(X, Y, yerr=2 * STD, color='r', label='Avg # of Tests Required')
plt.plot(X, z, color='g', label='K*np.log2(d+1)/np.log2(K+1)')
plt.xlabel("K")
plt.ylabel("Avergage Number of Tests")
plt.title("Number of Tests to Uniquely Solve")

# Adding legend, which helps us recognize the curve according to it's color
plt.legend()

# To load the display window
plt.show()