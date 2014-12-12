function [u, v, q, r] = pm4(w1, w2, w3, w4)

% [u, v, q, r] = pm4(w1, w2, w3, w4);
% u = (w1 + w2 - w3 - w4)/2;
% v = (w1 - w2 + w3 + w4)/2; 
% q = (w1 + w2 - w3 + w4)/2;
% r = (w1 + w2 + w3 - w4)/2;

u = (w1 - w2 - w3 - w4)/2;
v = (w1 - w2 + w3 + w4)/2; 
q = (w1 + w2 - w3 + w4)/2;
r = (w1 + w2 + w3 - w4)/2;
