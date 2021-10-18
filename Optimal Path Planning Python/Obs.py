import numpy as np

def obstacleMatrices(origin, theta, length, width):
    s = np.hstack([-np.cos(theta),-np.sin(theta)])
    d = np.hstack([np.cos(theta),np.sin(theta)])
    f = np.hstack([np.sin(theta),-np.cos(theta)])
    g = np.hstack([-np.sin(theta),np.cos(theta)])
    A = np.vstack([s,d,f,g])
    b = np.vstack([length / 2,length / 2,width / 2,width / 2])+np.matmul(A,origin)
    return A,b
