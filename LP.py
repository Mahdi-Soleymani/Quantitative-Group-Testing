import numpy as np
import scipy


# Import OR-Tools wrapper for linear programming
from ortools.linear_solver import pywraplp


## problem parameters
l=100
const=5
k=const*np.ones([l,])

for i in range(l):
    k[i]=np.random.randint(const)


k_sum=sum(k)
# Generate a random vector x
w=np.zeros([l,])
for i in range(l):
    w[i]=np.random.randint(k[i]+1)


print("unknown vector:")
print(w)

w_sum=sum(w)
print(w_sum)


# Create a solver using the GLOP backend
solver = pywraplp.Solver('Find (a) solurion(s)', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

# Create the variables we want to optimize
x=[]
for i in range(l):
    temp=solver.IntVar(int(0), int(k[i]), 'x')
    #temp=solver.NumVar(int(0), int(l), 'x')

    x.append(temp)


# Add constraints for each resource
constraint=0
for i in range(l):
    constraint=constraint+x[i]

solver.Add(constraint <= int(w_sum))
solver.Add(constraint >= int(w_sum))


# Maximize the objective function
cost_function=0
for i in range(l):
    cost_function=cost_function+x[i]

unique_solution=False
num_of_tests=0
beleifs= np.ones([l,])/l
previous_solution=np.zeros([l,])
while unique_solution==False:
    num_of_tests=num_of_tests+1
    ### Construct a new random test
    new_test=np.zeros([l,])
    for i in range(l):
        #new_test[i]=np.random.binomial(size=1, n=1, p= np.min([10*w_sum*beleifs[i],1]))
        #new_test[i]=np.random.binomial(size=1, n=1, p= w_sum*beleifs[i])
        new_test[i] = np.random.randint(2)
        #new_test[i]=np.random.binomial(size=1, n=1, p= .5)


    ##result of test
    new_result=np.inner(new_test, w)
    print("new test", new_test)
    print("new result", new_result)

    c=0
    for i in range(l):
        if new_test[i]==1:
            c = c + x[i]

    solver.Add(c <= int(new_result))
    solver.Add(c >= int(new_result))


    solver.Minimize(cost_function)


    # Solve problem
    status = solver.Solve()

    # If an optimal solution has been found, print results
    if status == pywraplp.Solver.OPTIMAL:
      print('================= Solution =================')
      print(f'Solved in {solver.wall_time():.2f} milliseconds in {solver.iterations()} iterations')

      print('Solution:')
      current_solution = np.zeros([l, ])
      for i in range(l):
          print('x'+str(i)+'value=',{x[i].solution_value()})
          current_solution[i]=(x[i].solution_value())

      if   sum(np.abs(previous_solution-current_solution))<=.01:
          unique_solution=True

      previous_solution=current_solution

    else:
      print('The solver could not find an optimal solution.')


    print("ground truth:")
    print(w)

    print("")
    print("Error:")
    print(sum(np.abs(current_solution-w)))


    print("number of tests:")
    print(num_of_tests)

    print('lower bound:')
    print((2*l/np.log2(l))*np.log2(1+(k_sum/l)))
