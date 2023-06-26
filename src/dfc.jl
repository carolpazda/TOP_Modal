#
# Calcula derivada de uma função por Diferenças Finitas Centrais
#
function DFC(f::Function, x0::Array, xmin::Array, xmax::Array, δ=1E-6)

    # Primeiro cuidado é identificar as dimensões do problema
 
    # Número de variáveis
    n = length(x0)
 
    # Primeira chamada da função
    f0 = f(x0)
 
    # Com isso, podemos dimensionar a saída da função
    D = zeros(n)
 
    # Loop entre as variáveis 
    for i=1:n
 
         # Perturba a variável i
         xb = x0[i]
 
         # Calculo para frente
         flag_frente = false
         if xb + δ <= xmax[i]
            flag_frente = true
            x0[i] = xb + δ
            ff = f(x0)
         end
 
         # Calculo para trás
         flag_tras = false
         if xb - δ >= xmin[i]
            flag_tras = true
            x0[i] = xb - δ
            fa = f(x0)
         end
 
         # Se for possível, calcula usando DFC
         if flag_frente && flag_tras
            
            D[i] = (ff-fa)/(2*δ)
    
         elseif flag_frente
            
            D[i] = (ff-f0)/(δ)
            
         elseif flag_tras
            
             D[i] = (f0-fa)/(δ)
             
         else
 
             println("Como isso é possível?!")
 
         end
 
         # Desfaz a perturbação nessa posição
         x0[i] = xb
 
    end #n
 
    # Retorna as derivadas
    return D
 
 end
 
 