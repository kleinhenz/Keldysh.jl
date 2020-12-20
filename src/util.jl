fermi(x) = x <= 0 ? (1 / (1 + exp(x))) : 1 - fermi(-x)
dfermi(x) = - fermi(x) * fermi(-x)

fermi(ϵ, β) = fermi(ϵ*β)
dfermi(ϵ, β) = β * dfermi(ϵ*β)
