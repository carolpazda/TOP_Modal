#
# Definição do operador <> que utilizamos no Lagrangiano Aumentado
#
function OP(a)
    max(0,a)
  end

#
# Função objetivo Lagrangiano Aumentado
#
function ObjetivoLa(f::Float64, gla::Vector{Float64},  c_atual::Float64, 
                    u_atual::Vector{Float64})

    # Realizamos um teste
    #length(gla)==length(u_atual)==m || throw("ObjetivoLa::dimensões das entradas estão erradas")

    # Somatório dos termos do LA
    La = f

    # Loop pelos termos das restrições
    for j=1:length(gla)

        La += (c_atual/2)*OP(u_atual[j]/c_atual + gla[j] )^2

    end #j

    # Retornamos a função La
    return La    

end

