#%% import
import osqp
from quadprog import solve_qp
# import quadprog
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import sparse
import scipy as sp
#import interpolate as ip
import socket, time
from scipy.linalg import block_diag
import scipy.sparse as sparse


def unicycleMatrice(xipred, n_upper, T):
    ad = np.zeros((n_upper, 4, 4))
    bd = np.zeros((n_upper, 4, 2))
    for i in range(n_upper):
        ad[i,0:4,0:4] = np.array([\
        [1,0,-T*xipred[3,i]*np.sin(xipred[2,i]), T*np.cos(xipred[2,i])],\
            [0, 1,T * xipred[3,i] * np.cos(xipred[2, i]),T * np.sin(xipred[2, i])],\
                [0,0,1,0],\
                    [0,0,0,1]])
        bd[i,0:4,0:2] = np.array([\
            [-(T**2*xipred[3,i]*np.sin(xipred[2,i]))/2, (T**2*np.cos(xipred[2,i]))/2], \
                [(T** 2 * xipred[3, i] * np.cos(xipred[2, i])) / 2, (T**2 * np.sin(xipred[2, i])) / 2],\
                    [T,0],\
                        [0,T]])
    return ad, bd


def LTVMPC2QP(n_upper,q, r, nstate, ninput, nposition, neta, upred, xipred, e, p, j_p,etamin,etamax,umin,umax,ad,bd):
    # print('called')
    pmat = np.zeros(((n_upper+1)*nstate+(n_upper+1)*ninput,(n_upper+1)*nstate+(n_upper+1)*ninput))

    for i in range(1,n_upper+1):
        pmat[i * nstate:(i + 1) * nstate, i * nstate:(i + 1) * nstate] = q

    I1blk = np.vstack([-np.eye(ninput),np.eye(ninput)])
    I1 = np.kron(np.eye(n_upper), I1blk)
    I2blk =np.vstack([np.eye(ninput),np.eye(ninput)])
    I2 = block_diag(np.eye(ninput), np.kron(np.eye(n_upper - 1), I2blk), np.eye(ninput))
    Rblkdiag = np.kron(np.eye(n_upper), r)

    pmat[(n_upper + 1) * nstate:pmat.shape[0], (n_upper + 1) * nstate:pmat.shape[1]] = \
        np.transpose(I2) @ I1 @ Rblkdiag @ np.transpose(I1) @ I2
    # np.matmul(np.matmul(np.matmul(np.matmul(np.transpose(I2),I1),Rblkdiag),np.transpose(I1))，I2)
    pmat=2*pmat


    qmat=np.zeros([(n_upper+1)*nstate+(n_upper+1)*ninput, 1])
# state section of q matrix
    for i in range(1,n_upper+1):
        error=np.vstack([e[0,i-1],e[1,i-1], 0, 0])
        qmat[i*nstate:(i+1)*nstate, 0]= (2*q@error).reshape(4)
# input section of q matrix
    for i in range(1, n_upper+2):
        # not sure
        if i==1:
            umiddle=np.array([upred[0:2, 1]-upred[0:2, 0]]).reshape(2,1)
            qmat[(n_upper+1)*nstate:(n_upper+1)*nstate+ninput, 0]= 2*np.matmul(-r,umiddle).reshape(2)
        elif i==n_upper+1:
            umiddle2 = np.array([upred[:, n_upper-1] - upred[:, n_upper]]).reshape(2,1)
            qmat[(n_upper+1) * nstate + (i-1)*ninput: (n_upper+1) * nstate + i*ninput, 0] \
                = 2*np.matmul(-r,umiddle2).reshape(2)
        elif i > 1 & i < n_upper+ 1:
            umiddle3=np.array([upred[:, i] - 2*upred[:, i-1]+upred[:, i-2]]).reshape(2,1)
            qmat[(n_upper+1) * nstate + (i-1)*ninput: (n_upper+1) * nstate + (i)*ninput, 0]\
                = 2 *np.matmul(-r,umiddle3).reshape(2)

    aeq= np.zeros((nstate+ninput+n_upper*nstate, (n_upper+1)*nstate+(n_upper+1)*ninput))

    aeq[0:nstate, 0:nstate]= np.eye(nstate)
    aeq[nstate+0:nstate+ninput, (n_upper+1)*nstate+0:(n_upper+1)*nstate+ninput]= np.eye(ninput)

    for i in range(n_upper):
        aeq[nstate+ninput+i*nstate:nstate+ninput+(i+1)*nstate, i*nstate:(i+1)*nstate]=ad[i]
        # print('ad',ad[i])
        aeq[nstate + ninput + i * nstate:nstate + ninput + (i + 1) * nstate, (n_upper+1)* nstate +\
            (i + 1) * ninput:(n_upper+1)*nstate+(i+2)*ninput] = bd[i]
        # print('bd',bd[i])
        aeq[nstate + ninput + i * nstate :nstate + ninput + (i + 1) * nstate, (i+1) * nstate :(i + 2) * nstate] =\
            -np.eye(nstate)

    c= np.zeros((nstate+ninput+n_upper*nstate, 1))
    aineq= np.zeros((n_upper*nstate+n_upper*ninput, (n_upper+1)*nstate+(n_upper+1)*ninput))

    for i in range(1,n_upper+1):
        aineq[(i-1)*nstate:i*nstate, i*nstate:(i+1)*nstate]= block_diag(j_p[0:2,0:2,i-1], np.eye(2))
        # print('middle',block_diag(j_p[0:2,0:2,i-1], np.eye(2)))
    aineq[n_upper*nstate+0:np.shape(aineq)[0], (n_upper+1)*nstate+ninput+0:np.shape(aineq)[1]]=\
        np.kron(np.eye(n_upper), np.eye(ninput))

    l=np.zeros((n_upper*nstate+n_upper*ninput,1))

    for i in range(1,n_upper+1):
        l[(i-1)*nstate:(i-1)*nstate+nposition, 0]=(np.array([[-0.18,-0.1]])-p[0:2,0, i-1])
        l[(i - 1) * nstate + nposition:(i - 1) * nstate + nposition+neta, 0]=\
            etamin-xipred[nposition:nposition+neta, i]

    for i in range(n_upper):
        l[n_upper*nstate+i*ninput:n_upper*nstate+i*ninput+ninput,0]=umin-upred[0:ninput,i+1]

    u=np.zeros((n_upper*nstate+n_upper*ninput,1))

    for i in range(1,n_upper+1):
        u[(i-1)*nstate:(i-1)*nstate+nposition,0] = (np.array([[0.18,0.1]])-p[0:2,0,i-1])
        u[(i-1)*nstate+nposition:(i-1)*nstate+nposition+neta,0]=etamax-xipred[nposition:nposition+neta,i]

    for i in range(n_upper):
        u[n_upper*nstate+i*ninput:n_upper*nstate+i*ninput+ninput,0]=umax-upred[0:ninput,i+1]
    # print('pmat',np.shape(pmat),'qmat',np.shape(qmat),'aineq',np.shape(aineq),'aeq',np.shape(aeq),'l',np.shape(l),'c',np.shape(c),'u',np.shape(u))
    # print('l',l,'u',u)
    for i in range(np.shape(l)[0]):
        if l[i]>u[i]:
            print('ind',i,'lower bound',l[i],'upper bound',u[i])
    # import numpy.linalg as LA
    # e_vals, e_vecs = LA.eig(pmat)
    # print('eigen values',e_vals)
    # print('eigen vectors',e_vecs)
    return pmat,qmat,aineq,l,u,aeq,c

def LTVPred(xi, u, T):

    dx = xi[3] * np.cos(xi[2])
    ddx = -xi[3] * np.sin(xi[2]) * u[0] + np.cos(xi[2]) * u[1]
    dy = xi[3] * np.sin(xi[2])
    ddy = xi[3] * np.cos(xi[2]) * u[0] + np.sin(xi[2]) * u[1]
    dtheta = u[0]
    dv = u[1]

    xiNext = np.zeros((4, 1))
    xiNext[0] = xi[0] + T * dx + T**2 / 2 * ddx
    xiNext[1] = xi[1] + T * dy + T**2 / 2 * ddy
    xiNext[2] = xi[2] + T * dtheta
    xiNext[3] = xi[3] + T * dv
    return xiNext

def QP2OSQP(pmat, qmat, aineq, l, u, aeq, c, n_upper, nstate, ninput):
    p = sparse.csc_matrix(pmat)
    # q=qmat.reshape(126)
    A_OSQP = sparse.vstack([aineq, aeq], format='csc')
    l_OSQP = np.vstack([l, c])
    u_OSQP = np.vstack([u, c])
    prob = osqp.OSQP()
    prob.setup(p, qmat, A_OSQP, l_OSQP, u_OSQP, verbose=True,
                                warm_start=True)
    res = prob.solve()
    # print('res',res.x[0:(n_upper+1)*nstate],'size',np.shape(res.x))

    if res.info.status != 'solved':
        raise ValueError('OSQP did not solve the problem!')

    delta_xiCol = res.x[0:(n_upper+1)*nstate]
    # print(delta_xiCol)
    # output = np.transpose(delta_xiCol)
    # delta_xi = output.reshape((nstate * (n_upper+1)), order='F')
    delta_xi = delta_xiCol.reshape((nstate, n_upper+1),order='F')
    # print(delta_xi)
    delta_uCol = res.x[(n_upper+1)*nstate:(n_upper+1)*nstate+(n_upper+1)*ninput]
    delta_u = delta_uCol.reshape((ninput, n_upper+1),order='F')
    # print('delta_xi',delta_xiCol,'delta_u',delta_uCol)
    return delta_xi, delta_u


def LTVMPC(q, r, n_upper, m_upper, T, xiref, centerpoints, xi, upred, constr):
    nstate = 4
    ninput = 2
    nposition = 2
    neta = nstate - nposition

    omega_guesscurr_stephis = np.vstack([upred[0,0:n_upper+1],np.zeros((m_upper+1, n_upper+1))])
    a_guesscurr_stephis = np.vstack([upred[1,0:n_upper+1],np.zeros((m_upper+1, n_upper+1))])

    xipred = np.hstack([xi,np.zeros((nstate, n_upper))])
    # xipred=np.zeros((xi.shape[0], n_upper + 1))
    # xipred[0:4, 0] = xi.reshape(4)
    # for i in range(n_upper):
    #     omega_guesscurr_stephis[0, i] = upred[0, i]
    #     a_guesscurr_stephis[0, i] = upred[1, i]

    # for j in range(m_upper - 1):
    #     if m_upper > 1:
    #         for i in range(n_upper):
    #             upred[0, i] = np.mean(omega_guesscurr_stephis[0:j, :], 0)
    #             upred[1, i] = np.mean(a_guesscurr_stephis[0:j, :], 0)



    j=0
    for i in range(1,n_upper+1):
        a = LTVPred(xipred[:, i-1], upred[:, i], T)
        xipred[:, i]=a.reshape(4)

    [ad, bd] = unicycleMatrice(xipred[:, 0:n_upper], n_upper, T)

    etamin = np.hstack([constr['thetaMin'],constr['vmin']])
    etamax = np.hstack([constr['thetaMax'], constr['vmax']])
    umin = np.hstack([constr['omegamin'], constr['amin']])
    umax = np.hstack([constr['omegamax'], constr['amax']])

    [e, p, j_p] = trackError(n_upper, xiref[:, 1:n_upper+1], centerpoints[:, 0:n_upper+1], xipred[:, 1:n_upper+1])
    # print('e',e,'p',p,'j_p',j_p)

    [p_qp, q_qp, aineq_qp, l_qp, u_qp, aeq_qp, c_qp] = \
    LTVMPC2QP(n_upper, q, r, nstate, nposition, neta, ninput, upred[:, 0:n_upper+1],
              xipred[:, 0:n_upper+1], e, p, j_p, etamin, etamax, umin, umax, ad, bd)

    [delta_xipred, delta_upred] = QP2OSQP(p_qp, q_qp, aineq_qp, \
                                              l_qp, u_qp, aeq_qp, c_qp, n_upper, nstate, ninput)

    upred[0:2, 1:n_upper + 1]= upred[0:2, 1:n_upper+1] + delta_upred[0:2, 1:n_upper+1]
    xipred[:, 1:n_upper+1] = xipred[:, 1:n_upper+1] + delta_xipred[:, 1:n_upper+1]

    a = saturation(upred[:, 1:n_upper+1], n_upper, constr)
    upred[:, 1:n_upper + 1]=a

    b = np.vstack([upred[0, 0:n_upper+1] - omega_guesscurr_stephis[j, 0:n_upper+1],upred[1,0:n_upper+1] - a_guesscurr_stephis[j,0:n_upper+1]])
    epscurr = np.linalg.norm(b,2)

    omega_guesscurr_stephis[j+1,:]=upred[0,0:n_upper+1]
    a_guesscurr_stephis[j+1,:]=upred[1,0:n_upper+1]

    return upred, xipred, epscurr

def saturation(xipred,n_upper,constr):
    for i in range(n_upper):
        if xipred[0, i]>constr['omegamax']:
            xipred[0, i]=constr['omegamax']
        elif xipred[0, i]<constr['omegamin']:
            xipred[0, i] = constr['omegamin']

        if xipred[1, i]>constr['amax']:
            xipred[1, i]=constr['amax']
        elif xipred[1, i]<constr['amin']:
            xipred[1, i] = constr['amin']
    return xipred

def trackError(n_upper, xiref, centerpoints, xipred):
    e = np.zeros((2, n_upper+1))
    p = np.zeros((2, 1, n_upper+1))
    j_p = np.zeros((2,2,n_upper+1))
    # print('centrepoints',centerpoints)
    for i in range(n_upper):
######################
       a=xipred[0:2, i]-xiref[0:2,i]
       e[0:2, i] = a.reshape(2)
       # print('e',e[0:2,i])
       # a=xipred[0,i]-xiref[0,i]
       # e[0,i]=a
       # e[1,i]=-a
       xDiff = np.zeros((2, 1))
       xDiff[0] = xipred[0, i]-centerpoints[0, i]
       xDiff[1] = xipred[1, i]-centerpoints[1, i]
       # print('DIFF',xDiff)
       angle = np.arctan((centerpoints[1, i+1]-centerpoints[1,i])/(centerpoints[0,i+1]-centerpoints[0,i]))
       if centerpoints[0,i+1]<centerpoints[0,i]:
           angle=angle-math.pi/2
       else:
           angle=angle+math.pi/2

       j_p[0,0,i]=np.cos(angle)
       j_p[0,1,i]=-np.sin(angle)
       j_p[1,0,i]=np.sin(angle)
       j_p[1,1,i]=np.cos(angle)
       # p[:,:,i]=j_p[:,:,i]@(xDiff.reshape(2,1))
       # print('jp',j_p[:,:,i])
       p[0,0,i]=np.cos(angle)*xDiff[0]-np.sin(angle)*xDiff[1]
       p[1,0,i] = np.sin(angle)*xDiff[0]+np.cos(angle)* xDiff[1]
       # print('p',p[:,:,i])
       # p[i,0,0]=np.array([np.cos(angle),-np.sin(angle)])*xDiff
       # p[i, 1, 0] = np.array([np.sin(angle), np.cos(angle)]) * xDiff
    return e, p, j_p

def unicycleODE(t, xi, v,w):
    dxidt = np.zeros((4,1))
    dxidt[0] = xi[3] * np.cos(xi[2])
    dxidt[1] = xi[3] * np.sin(xi[2])
    dxidt[2] = v
    dxidt[3] = w
    return dxidt.reshape(4,)




