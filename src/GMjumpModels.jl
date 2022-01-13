# ==============================================================================
# Modele JuMP pour calculer la relaxation linéaire du 2SPA sur un objectif donne, avec eventuellement une ϵ-contrainte
function computeLinearRelax2SPA(nbvar::Int, nbcontraintes::Int, L::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, epsilon, obj::Int)

    model = Model(with_optimizer(GLPK.Optimizer))
#    model = Model(with_optimizer(Gurobi.Optimizer, GUROBI_ENV))
#    set_optimizer_attributes(model, "OutputFlag" => 0)

    @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
    @constraint(model, [i=1:nbcontraintes],(sum((x[j]*L[i,j]) for j in 1:nbvar)) == 1)
    if obj == 1
        @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= epsilon)
    else
        @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= epsilon)
    end
    optimize!(model)
    return objective_value(model), value.(x)
end