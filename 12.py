import random
import numpy as np
import matplotlib.pyplot as plt
from functools import cache

def generate_random_binary_vector(n, k):
    if k > n:
        raise ValueError("Number of ones (k) cannot be greater than vector length (n)")

    # Create a list with k ones and (n-k) zeros
    vector = [1] * k + [0] * (n - k)

    # Shuffle the list to randomize the positions of ones and zeros
    random.shuffle(vector)

    return vector


def partition_and_find_max_ones(vector, l):
    if l < 1:
        raise ValueError("Number of partitions (l) must be at least 1")

    # Calculate the partition size
    partition_size = len(vector) // l

    # Initialize variables to track the maximum number of ones
    max_ones = 0
    num_of_zeros=0
    num_of_ones=0
    lis=[]
    # Iterate through the partitions
    for i in range(0, len(vector), partition_size):
        partition = vector[i:i + partition_size]
        ones_in_partition = sum(partition)
        lis.append(ones_in_partition)
        max_ones = max(max_ones, ones_in_partition)
    #     if ones_in_partition==0:
    #         num_of_zeros+=1
    #     if ones_in_partition==1:
    #         num_of_ones+=1
    # print("num of zeros:", num_of_zeros)
    # print("num of ones:", num_of_ones)
    # print("num of total zeros:", num_of_zeros*l)
    d = {x: lis.count(x) for x in lis}


    return d,max_ones

@cache
def num_of_tests_in_detecting_matrix(n,d):
    nu=2
    result=0
    while result<= n*np.log2(d+1):
        nu+=1
        result=(nu-2)*2**(nu-1)

    nu-=1
    return 2**nu

# Example usage:
n = 2**16# Length of the binary vector
alpha=.5
#k = int(n**alpha) # Number of ones
#l = int( k/np.log(k)) # Number of partitions
# l= int(k) # Number of partitions
k=2**12
# l=10
#coefficients=np.arange(1.7,min(3,n/k),0.2)
coefficients=np.arange(.8,min(3,n/k),0.2)
#coefficients=np.array([1/1.61, 1/1.65])

our_results=[]
our_result_errors=[]
ISIT=[]
iterations=10
for cof in coefficients:
    l=int(k*cof)
    print(cof)
    temp=[]
    for _ in range(iterations):
        binary_vector = generate_random_binary_vector(n, k)
        list_hist,max_ones_in_partition = partition_and_find_max_ones(binary_vector, l)
        summ=l

        for j in list_hist.keys():
            m_j=int(list_hist[j])
            # print(m_j)
            if j>0:
                if m_j>0:
                    ##using bound instead of real value
                    #summ+=np.ceil(j*np.log2(n/l))*min(m_j,2*m_j/np.log2(m_j))*np.ceil(np.log2(j))
                    #summ += (j * np.log2(np.ceil(n / l))) * min(m_j, (2 * m_j / np.log2(m_j)) * np.ceil(np.log2(j+1)))
                    ##real value for Bshouties construction
                    summ += (j * np.log2(np.ceil(n / l))) * min(m_j, num_of_tests_in_detecting_matrix(m_j,j))


        temp.append(summ)

    our_results.append(np.mean(temp))
    our_result_errors.append(np.std(temp))
    ISIT.append(k*np.log2(n/k))

print("k:",k)
plt.plot(coefficients,ISIT)
plt.errorbar(coefficients,our_results,yerr=our_result_errors)
plt.xlabel(("log(n)="+str(np.log2(n))+"log(k)="+str(np.log2(k))))
plt.show()
# print("Binary Vector:", binary_vector)
# print(f"Max Ones in Any Partition: {max_ones_in_partition}")
# print("2log(k):", 2*np.log(k))

# print("log2(n):", np.log2(n))
# print("sum:", sum)
# print("ISIT:", k*np.log2(n/k))
