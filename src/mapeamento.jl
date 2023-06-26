
#
# projeção não linear
#
function x2proj(x::Vector{Float64}, β, η, ρ_min)
   
    # Projeção heaviside
    ρ_proj = Map_ρ_tanh(β,η,x)

    # Correção do rho mínimo
    return ρ_min .+ (1.0-ρ_min).*ρ_proj


end


#
# Correção das derivadas
#
function dproj2dx!(x::Vector{Float64},dproj::Vector{Float64},β, η, ρ_min)

    # Operador de correção da projeção (matriz diagonal)
    R = dρ_projdρ(β,η,x) 

    # Corrige as derivadas devido ao filtro
    dx = R*dproj

    # Corrige para o valor mínimo
    dx .= dx.*(1.0-ρ_min)

    dproj.= dx

end

function Map_ρ_tanh(β::Float64,η::Float64,ρ::Vector{Float64})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    ρ_proj = Array{Float64}(undef,length(ρ))

    # Cte
    cte = tanh(β*η)

    cte2 = (cte+tanh(β*(1-η)))

    cte2 !=0.0 || throw("Map_ρ_tanh:: check β and η")

    # Loop sobre os elementos, calculando a projeção
    @inbounds for ele in LinearIndices(ρ)

        ρ_proj[ele] = (tanh(β*(ρ[ele]-η))+cte)/cte2

    end

    return ρ_proj

end

function dρ_projdρ(β::Float64,η::Float64,ρ::Vector{Float64})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    operador = Array{Float64}(undef,length(ρ))

    # Cte
    cte =  (tanh(β*η)+tanh(β*(1-η)))

    # Assertion
    cte != 0.0 || throw("dρ_projdρ:: check  β and η")

    # Loop sobre os elementos, calculando a projeção
    @inbounds for ele in LinearIndices(ρ)

         operador[ele] = (β*sech(β*(ρ[ele]-η))^2)/cte

    end
    
    return Diagonal(operador)

end