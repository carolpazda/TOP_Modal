# Vamos implementar a norma de Torii, que é dada por
#
#    (∑_i (1/λ)^p)^(-1/p)
#
# Vamos construir a norma de Torii, temos de pegar os 'x' menores valores devolvidos pela
# nossa função auxiliar e construir a norma.
# pd -> Termo que iremos receber para a norma.
# x_autovalores -> Número de autovalores que vamos trabalhar.
function Miolo_NormTorii(pd, lamb, x_autovalores)

    # Vamos começar com zero
    soma = 0.0

    # Vamos fazer isso dentro do for
    for i=1:x_autovalores

        # Vamos calcular o termo do somatório da norma de Torii
        soma += lamb[i]^(-pd)

    end

    # Retornando vamos ter
    return soma
end

#
# Vamos implementar a função objetivo, a mesma é dada por
#
#     (∑_i (1/λ)^p)^(-1/p) 
#
#

# Perceba que eu já implementei a norma de Torii, então vai ser bem simples
function Restricao_norma_Torii(soma, pd)

    # Definindo que
    T = soma^(-1/pd) 


end