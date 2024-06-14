import numpy as np
import numpy.polynomial.polynomial as npp
import pyldpc as ldpc
from scipy.sparse import rand
import scipy.sparse as sparse
import math
import scipy

class factor_graph:
    ##number of variable nodes
    n=32
    ## number of factor nodes
    m=1
    alpha=.5
    constant = 10
    d_v=1





    def __init__(self,n,alpha, constant):
        self.n=n
        self.alpha=alpha
        self.k = int(np.ceil(self.n ** self.alpha))
        #self.k=252
        ###sparse regime
        ######self.k=30
        ####sub-linear regime
        #self.k = int(np.ceil(self.n ** self.alpha))
        self.beliefs=(self.k/self.n)*np.ones([self.n,])
        self.constant=constant

        ######Constructing the measurment matrix######

        # constant gap to the information theoretic bound


##################   d_c should be an integer ########################
        ### d_c grows
        # d_c = int(2**2*self.n**(1-self.alpha)/constant)
        # # d_v is kinda fixed
        # d_v = int(2**2*np.ceil( (1 - self.alpha) / self.alpha))
        #
        # H, G = ldpc.make_ldpc(n, d_v, d_c, systematic=True, sparse=True)
        # self.m,q=np.shape(H)
        # self.d_v=d_v
        # self.d_c=d_c
        # self.H=H

        ##################   random parity-check/test matrix ########################
        ### d_c grows
        # d_c = int(np.floor(2 ** 3 * self.n ** (1 - self.alpha) / constant))
        # # d_v is kinda fixed
        # d_v = int(2 ** 3 * np.ceil((1 - self.alpha) / self.alpha))
        #
        # H, G = ldpc.make_ldpc(n, d_v, d_c, systematic=True, sparse=True)
        # self.m, q = np.shape(H)
        # self.d_v = d_v
        # self.d_c = d_c
        # self.H = H


        #### columnn regular matrix
        self.d_v=3
        self.m=int(constant*2*self.k*np.log2(self.n/self.k)/np.log2(self.k))
        # self.m=20
        H=np.zeros([self.m,n])
        for i in range(n):
            h_sparse = sparse.random(self.m, 1, density=self.d_v/self.m, data_rvs=np.ones)
            H[:,i]= h_sparse.toarray().reshape([self.m,])

        self.H = H



        #### row regular matrix
        # self.d_c=11
        # self.m=int(constant*2*self.k*np.log2(self.n/self.k)/np.log2(self.k))
        # # self.m=20
        # H=np.zeros([self.m,n])
        # for i in range(self.m):
        #     h_sparse = sparse.random(self.n, 1, density=self.d_c/self.n, data_rvs=np.ones)
        #     H[i, :]= h_sparse.toarray().reshape([self.n,])
        #
        # self.H = H



        #
        # self.m=int(constant*2*self.k*np.log2(self.n/self.k)/np.log2(self.k))
        # p=(np.log2(self.n)/self.n)*.85
        # H=np.random.choice([0, 1], size=(self.m,self.n), p=[1-p, p])
        #
        # self.H = H


















    # ## The belief on each variable node= prob of being equal to zero (not 1)
    # beliefs = np.random.random([n,])
    # belief_polynomials=[]
    # for i in range(n):
    #     belief_polynomials.append(npp.Polynomial([beliefs[i], 1-beliefs[i]]))





## generating the factor graph
fg=factor_graph(n=1000, alpha=.6, constant=4)
num_of_iterations=1000
### subset selection (multiplicative) gap
gap=1

## generate the random vector x
x_sparse= sparse.random(fg.n, 1, density=fg.k/fg.n, data_rvs=np.ones)
# Convert the sparse matrix to a full matrix
x=x_sparse.toarray()
y=np.matmul(fg.H,x)


v_node_adj=[]
# initialize the messages
messages_received_by_variable_nodes=[]
for j in range(fg.n):
    temp=[]
    message_temp_dict={}
    for i in range(fg.m):
        if fg.H[i,j]==1:
            temp.append(i)
            ### initializtion of messages
            #message_temp_dict[str(i)]=fg.k/fg.n
            # message_temp_dict[str(i)] = np.random.rand()
            message_temp_dict[str(i)] = .5

    v_node_adj.append(temp)
    messages_received_by_variable_nodes.append(message_temp_dict)


f_node_adj=[]
# initialize the messages
messages_received_by_factor_nodes=[]
for i in range(fg.m):
    message_temp_dict={}
    temp=[]
    for j in range(fg.n):
        if fg.H[i,j]==1:
            temp.append(j)
            ### initializtion of messages
            message_temp_dict[str(j)]=.5
    f_node_adj.append(temp)
    messages_received_by_factor_nodes.append(message_temp_dict)



######################  information #####

#
# print (np.count_nonzero(x.transpose()-prediction))
print('sparsity:')
print(np.count_nonzero(x))
print('test matrix size:')
print(fg.H.shape)
print('lower bound:')
print(2*fg.k*np.log2(fg.n/fg.k)/np.log2(fg.k))
#print(np.log2(scipy.special.binom(fg.n, fg.k))/np.log2(fg.k))
print('our ISIT paper:')
print(fg.k*np.log2(fg.n/fg.k)+fg.k)
# print(np.log2(scipy.special.binom(fg.n, fg.k)))
print('variable node degree:')
# print(fg.d_v)
# print('factor node degree:')
# print(fg.d_c)





belief_ave=0
t=0
diff_belief_ave=1
while t<num_of_iterations and np.abs(diff_belief_ave)>10**-10:

    print(t)
    t=t+1

    # print(messages_received_by_variable_nodes[0])
########## update messages emitting from variable nodes #########
    for j in range(fg.n):
        ### just compute the products
        probs1=np.array(list(messages_received_by_variable_nodes[j].values()))
        probs0=1-probs1
        product1=np.product(probs1)
        product0=np.product(probs0)
        result1=product1/(product0+product1)
        result0=1-result1


        ###update the messages emmiting from the variable nodes to the factor nodes (llr's)

        for i in v_node_adj[j]:
            prod0 = 1
            prod1 = 1
            for l in v_node_adj[j]:
                if l!=i:
                    prod1=prod1*messages_received_by_variable_nodes[j][str(l)]
                    prod0=prod0*(1-messages_received_by_variable_nodes[j][str(l)])

            messages_received_by_factor_nodes[i][str(j)]=prod1/(prod0+prod1)
            #messages_received_by_factor_nodes[i][str(j)]=llr_sum-messages_received_by_variable_nodes[j][str(i)]


    #update messages emitting from the factor nodes  ####
    for i in range(fg.m):
        # llr_vector=np.array(list(messages_received_by_factor_nodes[i].values())).flatten()
        # ####convert llr to p (prob of x =1)
        # p_vector=.5*(1+np.tanh(llr_vector))
        # prod_of_all_polynomials=1
        # ## constructiong the product of all polynomials
        # for l in range(len(p_vector)):
        #     prod_of_all_polynomials=prod_of_all_polynomials*npp.Polynomial([1-p_vector[l], p_vector[l]])

        for j in f_node_adj[i]:
            if y[i]==0:
                messages_received_by_variable_nodes[j][str(i)]=0 ##################
            elif y[i]==len(f_node_adj[i]):
                messages_received_by_variable_nodes[j][str(i)]=1  ##################
            else:

                result_of_multiplying_polynomials=1.0
                for w in f_node_adj[i]:
                    if w!= j:
                        prob= messages_received_by_factor_nodes[i][str(w)]
                        result_of_multiplying_polynomials=np.convolve(result_of_multiplying_polynomials,[1-prob, prob])
                # if len(result_of_multiplying_polynomials)==y[i]:
                #     messages_received_by_variable_nodes[j][str(i)] = 1 ################

                prob_be_1=result_of_multiplying_polynomials[int(y[i]-1)]+10**-12
                prob_be_0=result_of_multiplying_polynomials[int(y[i])]+10**-12
                # print("convolution:")
                # print(result_of_multiplying_polynomials)
                # print("")
                # print("y[i]")
                # print(y[i])
                # print("")

                messages_received_by_variable_nodes[j][str(i)]=prob_be_1/(prob_be_0+prob_be_1)


    belief = np.zeros([fg.n, 1])
    for j in range(fg.n):
        p1= np.prod(np.array(list(messages_received_by_variable_nodes[j].values())).flatten())+10**-6
        p0= np.prod(1-np.array(list(messages_received_by_variable_nodes[j].values())).flatten())+10**-6
        belief[j]=p1/(p0+p1)

    print("beleif_average:")
    print(np.average(belief))
    print
    diff_belief_ave=(np.average(belief) - belief_ave)
    belief_ave=np.average(belief)
    print("diff_llr_ave:")
    print(diff_belief_ave)


######  several types of hard decisions

    #### hard threshold
    #prediction = (np.sign(belief-.5)+1)/2

    # k-threshold
    sorted_ids=np.argsort(belief.reshape([fg.n,]))
    prediction=np.zeros([fg.n,1])

    for i in sorted_ids[fg.n-int(gap*fg.k):fg.n]:
        prediction[i]=1




    print("prediction")
    print(np.sum(np.abs(prediction)))

    print("prediction:")
    print(np.sum(np.abs(prediction)))
    print('mismatch:')
    print(np.count_nonzero(x - prediction))
    print("")

    # inconsistency=0
    # for j in range(fg.n):
    #     if prediction[j]==1:
    #         if x[j]!=1:
    #             inconsistency=inconsistency+1


    inconsistency=0
    for j in range(fg.n):
        if x[j] == 1:
            if prediction[j]!=1:
                inconsistency=inconsistency+1



    print("inconsistency:")
    print(inconsistency)
    print("")



    # prediction = []
    # for j in range(fg.n):
    #     prediction.append(((np.sign(sum(np.array(list(messages_received_by_variable_nodes[j].values())).flatten())) + 1) / 2))
    #
    #
    #
    # print("prediction llr's")
    # print(prediction)







# prediction=[]
# for j in range(fg.n):
#     prediction.append(((np.sign(sum(np.array(list(messages_received_by_variable_nodes[j].values())).flatten()))+1)/2))
#
#
#
# print('mismatch:')
#
# print (np.count_nonzero(x.transpose()-prediction))
print('sparsity:')
print(np.count_nonzero(x))
print('test matrix size:')
print(fg.H.shape)
print('lower bound:')
print(2*fg.k*np.log2(fg.n/fg.k)/np.log2(fg.k))
#print(np.log2(scipy.special.binom(fg.n, fg.k))/np.log2(fg.k))
print('our ISIT paper:')
print(fg.k*np.log2(fg.n/fg.k)+fg.k)
# print(np.log2(scipy.special.binom(fg.n, fg.k)))
print('variable node degree:')
print(fg.d_v)
print('beleifs:')
# print(belief.reshape([fg.n,]))
# print(fg.d_c)

# a=fg.belief_polynomials[0]*fg.belief_polynomials[1]
# #converts to vector of coefficient
# b=a.convert().coef
# #works fo multiplying polynomials so far!
# a=np.prod(fg.belief_polynomials)


#
# print(f_node_adj)
# print(messages_received_by_factor_nodes)
#
# print(v_node_adj)
# print(messages_received_by_variable_nodes)


