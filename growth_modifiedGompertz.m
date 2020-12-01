function output = growth_modifiedGompertz(P,t)

output = P(1) * exp( -exp(P(2)*exp(1)/P(1) .* (P(3)-t)+1) );

end
