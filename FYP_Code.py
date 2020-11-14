#!/usr/bin/env python3
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys

# globals
re_run = 0
blacklist = []
removed_lines = []


# m denotes the model
def refactor(m):
    global re_run, blacklist, removed_lines

    # initialising some objects
    consts = m.getConstrs()
    A1 = []
    numVars = len(m.getVars())
    for o in range(numVars):
        A1.append(0)

    # start of algorithm
    for i in range(len(consts)):
        for j in range(len(consts)):
            # to stop indexing out of range
            num_rem_vars = len(removed_lines)
            i -= num_rem_vars
            j -= num_rem_vars
            # if comparing same constraints, skip
            if i == j:
                pass
            elif 'Removed' in str(consts[i]).split("*") or 'Removed' in str(consts[j]).split("*"):
                pass
            else:
                match = 0
                # create sim eq's array for computation
                A = []
                for n in range(2):
                    A.append(A1.copy())
                B = [0, 0]

                # fill the sim eq's arrays
                for k in range(numVars):
                    A[0][k] = m.getCoeff(consts[i], m.getVars()[k])
                    B[0] = m.getAttr("rhs")[i]
                    A[1][k] = m.getCoeff(consts[j], m.getVars()[k])
                    B[1] = m.getAttr("rhs")[j]

                # simultaneous equation solver
                A = np.array(A)
                B = np.array(B)
                C = np.linalg.lstsq(A, B, rcond=None)[0]
                # checking to see that solution is in positive quadrant
                if not (C[0] > 0 and C[1] > 0):
                    # indicator set to catch that set of constraints found
                    match = 1

                # set of constraints found or the function is being re-run
                if match == 1 or re_run == 1:
                    # checks that the inequalities are opposite
                    if consts[i].getAttr("sense") != consts[j].getAttr("sense"):
                        # pick a constraint to be removed
                        if consts[i].ConstrName in blacklist:
                            rm_index = j
                        elif consts[j].ConstrName in blacklist:
                            rm_index = i
                        elif (consts[i].ConstrName and consts[j].ConstrName) in blacklist:
                            rm_index = "miss"
                        else:
                            # if modelsense is -1 maximise, if 1 minimise
                            if m.getAttr("ModelSense") == -1:
                                if consts[i].getAttr("sense") == '<':
                                    rm_index = j
                                else:
                                    rm_index = i
                            elif m.getAttr("ModelSense") == 1:
                                if consts[i].getAttr("sense") == '>':
                                    rm_index = j
                                else:
                                    rm_index = i

                        # if both consts aren't in the blacklist
                        if rm_index != "miss":
                            # deconstruct consts, add to blacklist and remove from model
                            c = consts[rm_index]
                            lhs, sense, rhs, name = m.getRow(c), c.Sense, c.RHS, c.ConstrName
                            blacklist.append(name)
                            m.remove(c)

                            # reset the model so it can be re-run with changes
                            m.reset(0)
                            # try optimisation again
                            try:
                                m.optimize()
                            except AttributeError:
                                print("Encountered an attribute error")
                                sys.exit(1)
                            # if solution unbounded then add constraint
                            # that was removed and re-run the algorithm
                            if m.getAttr("status") == GRB.UNBOUNDED:
                                print("unbounded solution, adding const and re-running")
                                m.addConstr(lhs, sense, rhs, name)
                                # update it to recognise re-add const
                                m.update()
                            else:
                                removed_lines.append(name)
                                blacklist.remove(name)

                            # if solution now found output it and exit
                            if m.SolCount > 0:
                                for v in m.getVars():
                                    print('%s %g' % (v.varName, v.x))
                                print('Obj: %g' % m.objVal)
                                print("Constraints removed were:")
                                for i in removed_lines:
                                    print(i)
                                sys.exit(1)

        # if algorithm has reached end without a match, re-run with marker
        if i == len(consts)-1 and j == len(consts)-1 and m.SolCount == 0:
            re_run += 1
            print("Re-running algorithm as model still infeasible")
            refactor(m)


def setUp():
    try:
        # Create a new model
        m = gp.Model("mip1")
        # quiets output
        m.setParam('OutputFlag', False)
        # helps granularity of outputs
        m.setParam('DualReductions', 0)

        # Create variables
        no_vars = input("How many variables do you have? ")

        x = m.addVar(vtype=GRB.CONTINUOUS, name="x")
        y = m.addVar(vtype=GRB.CONTINUOUS, name="y")
        z = m.addVar(vtype=GRB.CONTINUOUS, name="z")

        # to read the lp formatting from a file
        # m.read('input.lp')

        # Set objective
        m.setObjective(x + 3*y + 2*z, GRB.MAXIMIZE)

        # Add constraints
        m.addConstr(x + y + z >= 6, "c0")
        m.addConstr(2*x + y + 4*z <= 5, "c1")
        m.addConstr(y <= 3, "c2")
        m.addConstr(x + z <= 5, "c3")

        # Optimize model
        m.optimize()

        # output for succesful model
        for v in m.getVars():
            print('%s %g' % (v.varName, v.x))

        print('Obj: %g' % m.objVal)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    # error produced when no feasible region
    except AttributeError:
        # if infeasible and no current solution
        if m.getAttr("status") == GRB.INFEASIBLE and m.SolCount == 0:
            # execute algorithm
            refactor(m)
        else:
            print('Encountered an attribute error')


setUp()
