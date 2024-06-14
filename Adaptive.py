import random
import numpy as np
import matplotlib.pyplot as plt

def generate_random_binary_vector(n, k):
    if k > n:
        raise ValueError("Number of ones (k) cannot be greater than vector length (n)")

    # Create a list with k ones and (n-k) zeros
    vector = [1] * k + [0] * (n - k)

    # Shuffle the list to randomize the positions of ones and zeros
    random.shuffle(vector)

    return vector


def partition_and_test(vector, l,k):
    if l < 1:
        raise ValueError("Number of partitions (l) must be at least 1")

    # Calculate the partition size
    partition_size = max(1,len(vector) // l)

    # Initialize lis to track the number of defective items in each group
    num_of_def_items=[]
    new_vector=np.array([])
    number_of_groups_with_one_def_item=0
    # Iterate through the partitions
    for i in range(0, len(vector), partition_size):
        partition = vector[i:i + partition_size]
        ones_in_partition = sum(partition)
        num_of_def_items.append(ones_in_partition)
        if ones_in_partition ==1:
            number_of_groups_with_one_def_item+=1
        if ones_in_partition>1:
            #new_vector.append(vector[i:i + partition_size])
            new_vector=np.concatenate([new_vector, vector[i:i + partition_size]],axis=0)


    number_of_tests=min(number_of_groups_with_one_def_item*np.log2(partition_size),((2*number_of_groups_with_one_def_item)/np.log2(number_of_groups_with_one_def_item))*np.log2(partition_size))+l
        # max_ones = max(max_ones, ones_in_partition)
    #     if ones_in_partition==0:
    #         num_of_zeros+=1
    #     if ones_in_partition==1:
    #         num_of_ones+=1
    # print("num of zeros:", num_of_zeros)
    # print("num of ones:", num_of_ones)
    # print("num of total zeros:", num_of_zeros*l)
    #d = {x: num_of_def_items.count(x) for x in num_of_def_items}
    new_k=k-number_of_groups_with_one_def_item


    return number_of_tests,new_vector,new_k


# Example usage:
n = 2**16# Length of the binary vector
alpha=.5
#k = int(n**alpha) # Number of ones
k=100
#l=10
#coefficients=np.arange(.5,min(3,n/k),0.2)
coefficients=np.arange(1,3,.1)
#coefficients=np.array([ 1, 1.1, 1.3,2])
#coefficients=np.arange(1,min(1.2,n/k),0.2)


our_results=[]
our_result_errors=[]
ISIT=[]
iterations=100
threshold_n=10*k
threshold_k=10
stages=[]
stages_err=[]
tests=[]
tests_err=[]

for cof in coefficients:
    print(cof)
    stages_itr = []
    num_of_test_itr = []
    for _ in range(iterations):

        binary_vector = generate_random_binary_vector(n, k)
        num_of_test_list=[]
        new_k=k
        new_vector=np.array(binary_vector)
        num_of_stages=0
        # while (k>threshold_k) & (len(new_vector)>0):
        while   (len(new_vector) > 0):
            l = int(new_k * cof)
            number_of_tests,new_vector,new_k= partition_and_test(new_vector, l,new_k)
            # print("length:",len(new_vector))
            # print("k:", new_k)
            num_of_test_list.append(number_of_tests)
            binary_vector=random.shuffle(new_vector)
            num_of_stages+=1


        if len(new_vector)>0:
            num_of_test_list.append(new_k*np.log2(len(new_vector)))
        num_of_test_itr.append(sum(num_of_test_list))
        stages_itr.append(num_of_stages)

    tests.append(np.mean(num_of_test_itr))
    tests_err.append(np.std(num_of_test_itr))
    stages.append(np.mean(stages_itr))
    stages_err.append(np.mean(stages_itr))
    ISIT.append(k*np.log2(n/k))

print("k:",k)
print("ISIT", ISIT)
print("Total number of tests:", sum(num_of_test_list))
plt.errorbar(coefficients,tests, yerr=tests_err)

plt.plot(coefficients,ISIT)

plt.show()
# print("Binary Vector:", binary_vector)
# print(f"Max Ones in Any Partition: {max_ones_in_partition}")
# print("2log(k):", 2*np.log(k))

# print("log2(n):", np.log2(n))
# print("sum:", sum)
# print("ISIT:", k*np.log2(n/k))
