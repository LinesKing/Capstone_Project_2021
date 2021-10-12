import ipopt
import numpy as np
import math
from scipy.sparse import csr_matrix

def OBCAOneObsTimeOptimalPathPlanning(xInner,yInner,xOuter,yOuter,warm,initPoint,object, dMin):
    #N = np.shape(xOuter)[0]-1
    N=len(xOuter)-1#horizon
    tMax = 0.05
    vMin = 0
    thetaMin = - math.inf
    aMin = -5
    omegaMin = -2*math.pi
    vMax = 2
    thetaMax = math.inf
    aMax = 5
    omegaMax = 2*math.pi

    #define weights
    weightDist = 0.1
    weightVel = 0.5
    weightTheta = 0.1
    weightTime = 1

    #define constants
    nSystemState = 4
    nControlState = 2
    nObcaVariables = 4 #[lamda1-4]
    nTimeVariable = 1
    nStates = nSystemState+nControlState+nObcaVariables+nTimeVariable #[x,y,v,theta,a,w,lamda1,...,lamda4]
    print(nStates)
    #Initialize decision variables
    x0 = []
    options = {}
    options['lb'] = []
    options['ub']= []
    for i in range(N+1):
        a = np.array([warm['xWarm'][i],warm['yWarm'][i]])
        b = np.zeros(nStates-2,)
        x0 = np.hstack([x0,a,b])
        if i==0:
            options['lb'] = np.append(options['lb'],[initPoint[0],initPoint[1],vMin,initPoint[2],aMin,0])
            options['ub'] = np.append(options['ub'],[initPoint[0],initPoint[1],vMax,initPoint[2],aMax,0])
        else:
            yLb = np.min([yOuter[i],yInner[i]])
            yUb = np.max([yOuter[i],yInner[i]])
            xLb = np.min([xOuter[i],xInner[i]])
            xUb = np.max([xOuter[i],xInner[i]])
            options['lb'] = np.hstack([options['lb'], np.array([xLb,yLb,vMin,thetaMin,aMin,omegaMin])])
            options['ub'] = np.hstack([options['ub'], np.array([xUb,yUb,vMax,thetaMax,aMax,omegaMax])])

        options['lb'] = np.hstack([options['lb'], np.zeros((nObcaVariables))])
        options['ub'] = np.hstack([options['ub'], np.full((nObcaVariables), np.inf)])
        options['lb'] = np.hstack([options['lb'], 0])
        options['ub'] = np.hstack([options['ub'], tMax])
    print(len(x0),len(options['lb']),len(options['ub']))

    options['cl']=[]
    options['cu']=[]

    for i in range(N):
        options['cl'] = np.append(options['cl'], np.zeros((4, 1)))
        options['cu'] = np.append(options['cu'], np.zeros((4, 1)))

        options['cl'] = np.append(options['cl'],0)
        options['cu'] = np.append(options['cu'],0)

        options['cl'] = np.append(options['cl'],[dMin,-math.inf])
        options['cu'] = np.append(options['cu'],[math.inf,1])

    options['auxdata'] = {}
    options['auxdata']['N']=N
    options['auxdata']['nState']=nStates
    options['auxdata']['weightDist']= weightDist
    options['auxdata']['weightVel'] = weightVel
    options['auxdata']['weightTheta'] = weightTheta
    options['auxdata']['weightTime'] = weightTime
    options['auxdata']['xOuter'] = xOuter
    options['auxdata']['yOuter'] = yOuter
    options['auxdata']['xInner'] = xInner
    options['auxdata']['yInner'] = yInner
    options['auxdata']['ObjectA']=object['A']
    options['auxdata']['ObjectB'] = object['b']

    #
    # options.auxdata.weightDist = weightDist
    # options.auxdata.weightVel = weightVel
    # options.auxdata.weightTheta = weightTheta
    # options.auxdata.weightTime = weightTime
    # options.auxdata.xOuter = xOuter
    # options.auxdata.yOuter = yOuter
    # options.auxdata.xInner = xInner
    # options.auxdata.yInner = yInner


    class hs071():
        def __init__(self):
            pass

        def objective(self,x):

            # [N, nStates, weightDist, weightVel, weightTheta, weightTime] = options['auxdata'][0: 5].values
            f = 0
            k = 1
            for i in range(N):
                f = f + weightDist * (x[k] - x[k + nStates]) **2 \
                    + weightDist * (x[k + 1] - x[k + 1 + nStates]) **2 \
                        + weightVel * (x[k + 2] - x[k + 2 + nStates]) **2 \
                            + weightTheta * (x[k + 3] - x[k + 3 + nStates]) **2 \
                                + weightTime * x[k + 10]

                k = k + nStates
            return f

        def gradient(self,x):
            # [N, nStates, weightDist, weightVel, weightTheta, weightTime] = options['auxdata'][0: 5]
            g = np.zeros(((N + 1) * nStates, 1))
            k = 0; #Column
            m = 0; # Row

            for i in range(N + 1):
                if i == 1:
                    g[m, 0] = 2 * weightDist * (x[k] - x[k + nStates])
                    g[m + 1, 0] = 2 * weightDist * (x[k + 1] - x[k + 1 + nStates])
                    g[m + 2, 0] = 2 * weightVel * (x[k + 2] - x[k + 2 + nStates]);
                    g[m + 3, 0] = 2 * weightTheta * (x[k + 3] - x[k + 3 + nStates]);
                    g[m + 10, 0] = weightTime;
                elif i == N:
                    g[m, 0] = -2 * weightDist * (x[k - nStates] - x[k]);
                    g[m + 1, 0] = -2 * weightDist * (x[k + 1 - nStates] - x[k + 1])
                    g[m + 2, 0] = -2 * weightVel * (x[k + 2 - nStates] - x[k + 2])
                    g[m + 3, 0] = -2 * weightTheta * (x[k + 3 - nStates] - x[k + 3])
                    g[m + 10, 0] = weightTime

                else:
                    g[m, 0] = -2 * weightDist * (x[k - nStates] - x[k]) + 2 * weightDist * (x[k] - x[k + nStates])
                    g[m + 1, 0] = -2 * weightDist * (x[k + 1 - nStates] - x[k + 1]) + 2 * weightDist * (x[k + 1] - x[k + 1 + nStates])
                    g[m + 2, 0] = -2 * weightVel * (x[k + 2 - nStates] - x[k + 2]) + 2 * weightVel * (x[k + 2] - x[k + 2 + nStates])
                    g[m + 3, 0] = -2 * weightTheta * (x[k + 3 - nStates] - x[k + 3]) + 2 * weightTheta * (
                                x[k + 3] - x[k + 3 + nStates])
                    g[m + 10, 0] = weightTime

                m = m + nStates
                k = k + nStates
            return g

        def constraints(self,x):
            # [N, nStates, weightDist,weightVel,weightTheta, weightTime, xOuter, yOuter, xInner, yInner, A, b] = options['auxdata'][:]

            nConstraintsSystem = 4
            nConstraintsTriplet = 1
            nConstraintsObject = 2
            nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject
            A = options['auxdata']['ObjectA']
            b = options['auxdata']['ObjectB']
            c = np.zeros((nConstraints * N, 1))
            k = 0 # Column
            m = 0 # Row
            for i in range(N):
                c[m, 0] = x[k + nStates] - x[k] - x[k + 10] * x[k + 2] * np.cos(x[k + 3])
                c[m + 1, 0] = x[k + 1 + nStates] - x[k + 1] - x[k + 10]* x[k + 2] * np.sin(x[k + 3])
                c[m + 2, 0] = x[k + 2 + nStates] - x[k + 2] - x[k + 10] * x[k + 4]
                c[m + 3, 0] = x[k + 3 + nStates] - x[k + 3] - x[k + 10] * x[k + 5]
                m = m + nConstraintsSystem

                if xInner[i] == xOuter[i] or i == 1:
                    c[m, 0] = 0
                else:
                    c[m, 0] = x[k + 1] - x[k] * (yOuter[i] - yInner[i]) / (xOuter[i] - xInner[i]) \
                              - (yInner[i] - (yOuter[i] - yInner[i]) / (xOuter[i] - xInner[i]) * xInner[i])
                m = m + nConstraintsTriplet

                c[m, 0] = x[k+6] * (A[1-1, 1-1] * x[k] + A[1-1, 2-1] * x[k+1] - b[1-1, 1-1])+ x[k+7] * (A[2-1,\
                    1-1] * x[k] + A[2-1, 2-1] * x[k+1] - b[2-1, 1-1]) + x[k+8] * (A[3-1, 1-1] * x[k] + A[3-1,\
                        2-1] * x[k+1] - b[3-1, 1-1]) + x[k+9] * (A[4-1, 1-1] * x[k] + A[4-1, 2-1] * x[k+1] - b[4-1, 1-1])

                c[m+1, 0]= (x[k+6] * A[1-1, 1-1] + x[k+7] * A[2-1, 1-1] + x[k+8] * A[3-1, 1-1] + x[k+9] * A[4-1,\
                    1-1])**2 + (x[k+6] * A[1-1, 2-1] + x[k+7] * A[2-1, 2-1] + x[k+8] * A[3-1, 2-1] + x[k+9] * A[4-1,\
                        2-1])**2

                m = m + nConstraintsObject
                k = k + nStates
            return c

        def jacobian(self,x):
            # [N, nStates, weightDist, weightVel, weightTheta, weightTime, xOuter, yOuter, xInner, yInner, A, b] = options['auxdata'][:]
            nConstraintsSystem = 4
            nConstraintsTriplet = 1
            nConstraintsObject = 2
            nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject
            A = options['auxdata']['ObjectA']
            ############
            from scipy.sparse import csc_matrix
            J = csc_matrix((nConstraints * N, nStates * (N + 1)))
            ##########

            k = 0 # Column
            m = 0 # Row

            for i in range(N):
                J[m, k] = -1
                J[m, k + 2] = -x[k + 10] * np.cos(x[k + 3])
                J[m, k + 3] = x[k + 10] * x[k + 2] * np.sin(x[k + 3])
                J[m, k + nStates] = 1
                J[m, k + 10] = -x[k + 2] * np.cos(x[k + 3])
                J[m + 1, k + 1] = -1
                J[m + 1, k + 2] = -x[k + 10] * np.sin(x[k + 3])
                J[m + 1, k + 3] = -x[k + 10] * x[k + 2] * np.cos(x[k + 3])
                J[m + 1, k + 1 + nStates] = 1
                J[m + 1, k + 10] = -x[k + 2] * np.sin(x[k + 3]);
                J[m + 2, k + 2] = -1
                J[m + 2, k + 4] = -x[k + 10]
                J[m + 2, k + 2 + nStates] = 1
                J[m + 2, k + 10] = -x[k + 4]
                J[m + 3, k + 3] = -1
                J[m + 3, k + 5] = -x[k + 10]
                J[m + 3, k + 3 + nStates] = 1
                J[m + 3, k + 10] = -x[k + 5]
                m = m + nConstraintsSystem
                if xInner[i] != xOuter[i]:
                    J[m, k] = -(yOuter[i] - yInner[i]) / (xOuter[i] - xInner[i])
                    J[m, k + 1] = 1

                m = m + nConstraintsTriplet
               # print(np.shape(b))

                J[m, k]= x[k + 6] * A[0, 0] + x[k + 7] * A[1, 0] + x[k + 8] * A[2, 0] + x[k + 9] * A[3, 0]
                J[m, k + 1] = x[k + 6] * A[0, 1] + x[k + 7] * A[1, 1] + x[k + 8] * A[2, 1] + x[k + 9] * A[3, 1]
                J[m, k + 6] = A[1-1, 1-1] * x[k] + A[1-1, 2-1] * x[k + 1] - b[1-1]
                J[m, k + 7] = A[2-1, 1-1] * x[k] + A[2-1, 2-1] * x[k + 1] - b[2-1]#, 1-1]
                J[m, k + 8] = A[3-1, 1-1] * x[k] + A[3-1, 2-1] * x[k + 1] - b[3-1]#, 1-1]
                J[m, k + 9] = A[4-1, 1-1] * x[k] + A[4-1, 2-1] * x[k + 1] - b[4-1]#, 1-1]

                J[m + 1, k + 6] = 2 * A[1-1, 1-1] * (x[k + 6] * A[1-1, 1-1] + x[k + 7] * A[2-1, 1-1] + x[k + 8] * A[3-1\
                , 1-1] + x[k + 9]* A[4-1, 1-1]) + 2 * A[1-1, 2-1] * (x[k + 6] * A[1-1, 2-1] + x[k + 7] * A[1, 1] +\
                        x[k + 8] * A[3-1, 2-1] + x[k + 9] * A[4-1, 2-1])
                J[m + 1, k + 7] = 2 * A[1, 0] * (x[k + 6] * A[0, 0] + x[k + 7] * A[1, 0] + x[k + 8] * A[2, 0] + x[k + 9]\
                                                 * A[3, 0]) + 2 * A[1, 1] * (x[k + 6] * A[0, 1]+ x[k + 7] * A[1, 1] \
                                                                             + x[k + 8] * A[2, 1] + x[k + 9] * A[3, 1])
                J[m + 1, k + 8] = 2 * A[2, 0] * (x[k + 6] * A[0, 0] + x[k + 7] * A[1, 0]+ x[k + 8] * A[2, 0]\
                                                 + x[k + 9] * A[3, 0]) + 2 * A[2, 1] * (x[k + 6] * A[0, 1]+ x[k + 7] \
                                                                * A[1, 1] + x[k + 8] * A[2, 1] + x[k + 9] * A[3, 1])
                J[m + 1, k + 9] = 2 * A[3, 0] * (x[k + 6] * A[0, 0] + x[k + 7] * A[1, 0] + x[k + 8] * A[2, 0] + \
                                                 x[k + 9] * A[3, 0]) + 2 * A[3, 1] * (x[k + 6] * A[0, 1] + x[k + 7] \
                                                                * A[1, 1] + x[k + 8] * A[2, 1] + x[k + 9] * A[3, 1]);

                m = m + nConstraintsObject
                k = k + nStates
            return J

        def J_struct(self,x):
            # [N, nStates] = options['auxdata'][1: 2]
            nConstraintsTriplet = 1
            nConstraintsSystem = 4
            nConstraintsObject = 2
            nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject

            J = np.zeros((nConstraints * N, nStates * (N + 1)));
            k = 0
            m = 0
            for i in range(N):
                J[m, k] = 1
                J[m, k + 2] = 1
                J[m, k + 3] = 1
                J[m, k + nStates] = 1
                J[m, k + 10] = 1
                J[m + 1, k + 1] = 1
                J[m + 1, k + 2] = 1
                J[m + 1, k + 3] = 1
                J[m + 1, k + 1 + nStates] = 1
                J[m + 1, k + 10] = 1
                J[m + 2, k + 2] = 1
                J[m + 2, k + 4] = 1
                J[m + 2, k + 2 + nStates] = 1
                J[m + 2, k + 10] = 1
                J[m + 3, k + 3] = 1
                J[m + 3, k + 5] = 1
                J[m + 3, k + 3 + nStates] = 1
                J[m + 3, k + 10] = 1
                m = m + nConstraintsSystem

                J[m, k] = 1
                J[m, k + 1] = 1
                m = m + nConstraintsTriplet

                J[m, k] = 1
                J[m, k + 1] = 1
                J[m, k + 6] = 1
                J[m, k + 7] = 1
                J[m, k + 8] = 1
                J[m, k + 9] = 1
                J[m + 1, k + 6] = 1
                J[m + 1, k + 7] = 1
                J[m + 1, k + 8] = 1
                J[m + 1, k + 9] = 1
                m = m + nConstraintsObject
                k = k + nStates
            J_stru = csr_matrix(J)
            return J_stru

    # print(len(x0),len(options['lb']),len(options['ub']))
    nlp = ipopt.problem(
        n=len(x0),
        m=len(options['cl']),
        problem_obj=hs071(),
        lb=options['lb'],
        ub=options['ub'],
        cl=options['cl'],
        cu=options['cu']
    )
    # nlp.add_option('mu_strategy', 'adaptive')
    # nlp.add_option('tol', 1e-7)

    sol, info = nlp.solve(x0)
    x = sol[0:(N + 1) * nStates].reshape((nStates, (N + 1)),order='F')
    return x






