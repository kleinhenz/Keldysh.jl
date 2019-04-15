function fermi(epsilon, beta)
  x = epsilon * beta
  x < 0.0 ? (1.0 / (1.0 + exp(x))) : (1.0 - 1.0 / (1.0 + exp(-x)))
end

