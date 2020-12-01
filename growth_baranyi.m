function output = growth_baranyi(P,t)
 
A = t + (1/P(1)) * log(exp(-P(1)*t) + exp(-P(1)*P(2)) - exp(-P(1)*(t+P(2))));
output = P(4)  + P(1) * A - log(1 + (exp(P(1)*A) - 1)/exp(P(3)));
 
end