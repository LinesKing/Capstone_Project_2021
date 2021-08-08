function [th] = theta(x,s,bL,bU)
    th = sum(abs(res(x,s,bL,bU)));
