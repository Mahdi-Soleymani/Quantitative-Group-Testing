import random
import numpy as np
import matplotlib.pyplot as plt
import scipy.special


def num_of_tests_in_detecting_matrix(d):
   #d is the list of integers represnting maxs for each coordinate

    d=np.array(d)+1
    d=sorted(d)
    n=len(d)
    nu_max=int(np.ceil(np.log2(n)))
    nu_f=nu_max

    for nu in range(2,nu_max+1):
        if (nu-2)*(2**(nu-1))>((2**nu)*np.log2(d[int(n-2**nu-1)])+np.sum((np.log2(d[0:int(n-2**nu-1)])))):
            nu_f=nu-1
            break

    return 2**nu_f


## ??????????????????
def num_of_tests_in_detecting_matrix_improved(d):
   #d is the list of integers represnting maxs for each coordinate

    d=np.array(d)+1
    d=sorted(d)
    n=len(d)
    k_max=n
    k_f = n
    for k in range(2,k_max):
        if  ((np.log2(k)-4)*k)/2 > k *np.log2(d[int(n-k-1)])+np.sum((np.log2(d[0:int(n-k-1)]))):
            k_f = k - 1
            break

    return k_f







def num_of_tests_required_to_resolution(vector,resolution):
    # resolution is the size of groups
    n=len(vector)
    num_of_rounds=int(np.ceil(np.log2(n/resolution)))
    num_of_tests=0
    for i in range(num_of_rounds):
        partition_size=max(1,n//2**(i+1))
        d=[]
        for j in range(0, n, 2*partition_size):
            partition = vector[j:j + partition_size]
            ones_in_partition = sum(partition)
            d.append(ones_in_partition)
        num_of_tests+=num_of_tests_in_detecting_matrix(d)
    print("last group size:", partition_size)
    return num_of_tests

def num_of_tests_required_from_agroupsize_with_certain_number_of_stages(vector,group_size1,num_of_stages):
    # resolution is the size of groups
    n = len(vector)
    num_of_tests = np.ceil(n/group_size1)-1

    num_of_rounds=num_of_stages-2
    if group_size1>=n:
        if num_of_stages >= np.log2(n):
            num_of_rounds=num_of_stages
        else:
            num_of_rounds=num_of_stages-1
    else:
        if num_of_stages-1>=np.log2(group_size1):
            num_of_rounds=num_of_stages-1

    final_size=group_size1
    for i in range(num_of_rounds):
        partition_size = max(1, group_size1//2**(i+1))
        d = []
        for j in range(0, n, 2 * partition_size):
            partition = vector[j:j + partition_size]
            ones_in_partition = sum(partition)
            d.append(ones_in_partition)
        num_of_tests += num_of_tests_in_detecting_matrix(d)
        final_size=partition_size
    print("final size:", final_size)
    print("number of stages:", num_of_stages)
    return num_of_tests,final_size







def generate_random_binary_vector(n, k):
    if k > n:
        raise ValueError("Number of ones (k) cannot be greater than vector length (n)")

    # Create a list with k ones and (n-k) zeros
    vector = [1] * k + [0] * (n - k)

    # Shuffle the list to randomize the positions of ones and zeros
    random.shuffle(vector)

    return vector


def partition_and_find_max_ones(vector, partition_size):
    if partition_size < 1:
        raise ValueError("Number of partitions (l) must be at least 1")



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
    d = {x: lis.count(x) for x in lis}


    return d,max_ones


def num_of_measurements_in_sparse_recovery(k,size):
    if size==1:
        return 0
    else:
        return min(size,np.ceil(k*np.log2(size)))

def adaptive_kronecker_with_stage(n,k,num_of_stages):
    # Example usage:
    # n = 2**16# Length of the binary vector
    # k=500
    #group_sizes=np.linspace(100,2000,num=50,dtype=int)
    power=np.arange(num_of_stages,np.log2(n)+1,dtype=int)
    group_sizes=2**power
    #group_sizes=np.linspace(4,50,num=20,dtype=int)
    print(group_sizes)
    #group_sizes=np.linspace((2**num_of_stages),n,num=10,dtype=int)
    print("group sizes:", group_sizes)
    #group_sizes=[n]
    #group_sizes=2**(np.arange(0,16))


    our_results=[]
    our_result_errors=[]
    ISIT=[]
    iterations=100
    for size in group_sizes:
        #print("group size:", size)
        temp=[]
        for _ in range(iterations):
            binary_vector = generate_random_binary_vector(n, k)
            summ,final_size=num_of_tests_required_from_agroupsize_with_certain_number_of_stages(binary_vector,size,num_of_stages)
            #print("number of measurments needed in the adaptive stage to achieve the resolution",summ)
            #final_size=int(np.ceil(size/2**(num_of_stages-2)))
            list_hist,max_ones_in_partition = partition_and_find_max_ones(binary_vector, final_size)
            #print("summ:",summ)
            for j in list_hist.keys():
                m_j=int(list_hist[j])
                # print(m_j)
                if j>0:
                    if m_j>0:
                        ##using bound instead of real value
                        #summ+=np.ceil(j*np.log2(n/l))*min(m_j,2*m_j/np.log2(m_j))*np.ceil(np.log2(j))
                        #summ += (j * np.log2(np.ceil(n / l))) * min(m_j, (2 * m_j / np.log2(m_j)) * np.ceil(np.log2(j+1)))
                        ##real value for Bshouties construction
                        summ += num_of_measurements_in_sparse_recovery(j,final_size)* num_of_tests_in_detecting_matrix([j]*m_j)
            temp.append(summ)

        our_results.append(np.mean(temp))
        our_result_errors.append(np.std(temp))
        ISIT.append(k*np.log2(n/k))

    if False:
        print("k:",k)
        plt.plot(group_sizes,ISIT)
        plt.errorbar(group_sizes,our_results,yerr=our_result_errors)
        plt.xlabel(("log(n)="+str(np.log2(n))+"log(k)="+str(np.log2(k))))
        print(our_results)
        print("best:",min(our_results), "ratio:",ISIT[0]/min(our_results) )
        print(ISIT)
        print("group sizes:", group_sizes)

        plt.show()
    best=min(our_results)
    min_index=np.argmin(best)
    error=our_result_errors[min_index]
    return best,error, ISIT[0]



def adaptive_kronecker():
    # Example usage:
    n = 2**16# Length of the binary vector
    k=100
    group_sizes=np.linspace(10000,30000,num=30,dtype=int)
    #group_sizes=2**(np.arange(0,9))

    our_results=[]
    our_result_errors=[]
    ISIT=[]
    iterations=10
    for size in group_sizes:
        print("group size:", size)
        temp=[]
        for _ in range(iterations):
            binary_vector = generate_random_binary_vector(n, k)
            summ=n//size
            #summ=min(l,((2*k)/np.log2(k))*np.log2(l))
            #summ=num_of_tests_required_to_resolution(binary_vector, size)
            print("number of measurments needed in the adaptive stage to achieve the resolution",summ)

            list_hist,max_ones_in_partition = partition_and_find_max_ones(binary_vector, size)

            for j in list_hist.keys():
                m_j=int(list_hist[j])
                print("m_j:",m_j,"j:",j)
                # print(m_j)
                if j>0:
                    if m_j>0:
                        ##using bound instead of real value
                        #summ+=np.ceil(j*np.log2(n/l))*min(m_j,2*m_j/np.log2(m_j))*np.ceil(np.log2(j))
                        #summ += (j * np.log2(np.ceil(n / l))) * min(m_j, (2 * m_j / np.log2(m_j)) * np.ceil(np.log2(j+1)))
                        ##real value for Bshouties construction
                        summ += num_of_measurements_in_sparse_recovery(j,size)*  num_of_tests_in_detecting_matrix([j]*m_j)
            temp.append(summ)

        our_results.append(np.mean(temp))
        our_result_errors.append(np.std(temp))
        ISIT.append(k*np.log2(n/k))
        print(list_hist)

    print("k:",k)
    plt.plot(group_sizes,ISIT)
    plt.errorbar(group_sizes,our_results,yerr=our_result_errors)
    plt.xlabel(("log(n)="+str(np.log2(n))+"log(k)="+str(np.log2(k))))
    print(our_results)
    print("best:",min(our_results), "ratio:",ISIT[0]/min(our_results) )
    print(ISIT)
    plt.show()


def adaptive_old(n,k):
        # Example usage:
        # n = 2 ** 16  # Length of the binary vector
        # k = 100

        group_sizes=np.array([n])

        our_results = []
        our_result_errors = []
        ISIT = []
        iterations = 100
        for size in group_sizes:

            temp = []
            for _ in range(iterations):
                binary_vector = generate_random_binary_vector(n, k)
                summ = num_of_tests_required_to_resolution(binary_vector,  1)
                print("summ", summ)


                temp.append(summ)

            our_results.append(np.mean(temp))
            our_result_errors.append(np.std(temp))
            ISIT.append(k * np.log2(n / k))

        print("k:", k)
        plt.plot(group_sizes, ISIT)
        plt.errorbar(group_sizes, our_results, yerr=our_result_errors)
        plt.xlabel(("log(n)=" + str(np.log2(n)) + "log(k)=" + str(np.log2(k))))
        print(our_results)
        print("best:", min(our_results), "ratio:", ISIT[0] / min(our_results))
        print(ISIT)
        info_bound=(np.log2(scipy.special.comb(n,k,exact=False)))/np.log2(k)
        print("info theory lowerbound:", info_bound)
        print("ratio to info bound:", )
        plt.show()


def zeros_and_max():
    #plots the number of groups with zeros items and mor than 1 item versus the group size
    #plots the number of measurments in a non-adaptive
    n=2**16
    k=150
    vector=generate_random_binary_vector(n, k)
    #group_sizes=(np.linspace(k/np.log(k), 500, num=500, dtype=int))
    group_sizes=[411]
    num_of_zeros=[]
    num_of_more_than_one=[]
    maxs=[]
    for size in group_sizes:
        zeros=0
        more_ones=0
        max_ones=0
        for j in range(0,n,size):
            num_of_ones=sum(vector[j:j+size])
            if num_of_ones==0:
                zeros+=1
            if num_of_ones > 1:
                more_ones+=1
            if num_of_ones>max_ones:
                max_ones=num_of_ones
        num_of_zeros.append(zeros)
        num_of_more_than_one.append(more_ones)
        maxs.append(max_ones)

    # plt.plot(group_sizes,num_of_zeros)
    # plt.plot(group_sizes,num_of_more_than_one)
    # plt.show()
    print("zeros:", num_of_zeros)
    print("more than one:", num_of_more_than_one)
    print("group sizes:", group_sizes)
    print("maxes:", maxs)
    print("bound:", 2*np.log(k))

    num_of_measurments=[]
    ISIT=[]
    for i in range(len(group_sizes)):
        tests=maxs[i]*np.log2(group_sizes[i]+1)*num_of_tests_in_detecting_matrix([maxs[i]]*int(np.ceil(n/group_sizes[i])))
        num_of_measurments.append(tests)
        ISIT.append(k * np.log2(n / k))

    plt.plot(group_sizes,num_of_measurments)
    plt.plot(group_sizes,ISIT)
    plt.show()
    print("best we can do",min(num_of_measurments))
    print("ISIT",ISIT[0])
    print(group_sizes[int(np.argmin(num_of_measurments))])








if __name__ == '__main__':
    #adaptive_kronecker()
    #adaptive_old(2**16, 100)
    m=16
    n=2**m
    k=500
    ISIT=[]
    ours=[]
    Bsh=[]
    error=[]
    x_axis=np.arange(2,m+1,1)
    for i in x_axis:
        us,err,isit=adaptive_kronecker_with_stage(n, k, i)
        ours.append(us)
        ISIT.append(isit)
        Bsh.append(2*isit/np.log2(k))
        error.append(err)

    plt.plot(x_axis,Bsh,linestyle='dashed',label='Bshouty, Fully adaptive')
    plt.legend(loc='right')
    plt.plot(x_axis,ISIT,linestyle='dashed',label='Non-adaptive, e.g., ISIT paper')
    plt.legend(loc='right')
    plt.errorbar(x_axis,ours,marker='*',yerr=error,label='Our new algorithm')
    plt.legend(loc='right')

    print(Bsh)
    print(ISIT)
    print(ours)
    print(error)


    plt.xlabel('Num of stages')
    plt.ylabel('Num of tests')
    plt.title('k='+str(k)+' n='+str(n))

    plt.show()



    num_of_measurements=[]
    num_of_measurements_improved=[]
    info_bound=[]
    for i in range(1000,10000,1000):
        a=[1]*i
        num_of_measurements.append(num_of_tests_in_detecting_matrix(a))
        num_of_measurements_improved.append(num_of_tests_in_detecting_matrix_improved(a))
        info_bound.append(2*i/np.log2(i))
        print(i)
    plt.plot(num_of_measurements)
    plt.plot(np.array(num_of_measurements_improved))
    plt.plot(info_bound)
    plt.show()


