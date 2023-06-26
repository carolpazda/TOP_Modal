# 
#   ESSA ROTINA SERVE PARA O CÁLCULO DE TODAS AS ITERAÇÕES DO LA
#   Utiliza o pacote WallE e o La/DerLa presentes no main
#
function Calculos_LA(x, u_atual, c_atual, c_limite,  La::Function, g::Function, 
                    DerLa::Function, b, numero_iteracoes_LA=50)
        
    # Restrições laterais do problema de otimização                
    ci = 1E-3*ones(length(x))
    cs = ones(length(x))
    
    # Define fora do loop para poder devolver depois
    rest_at = zeros(length(u_atual))

    # Restrições no começo do problema
    rest_at .= g(x)

    # Violação máxima no começo do problema
    max_viola = maximum(rest_at[rest_at.>0.0],init=0.0)

    # Loop externo do LA
    for k=1:numero_iteracoes_LA

        # Usando o WallE (Steepest com Bloqueio) -> Pacote carregado do GitHub
        options = WallE.Init()
        options["NITER"] = 1000
        output = WallE.Solve(La,DerLa,x,ci,cs,options)
        xe = output["RESULT"]
        flag_converged = output["CONVERGED"]
        opt_norm = output["NORM"]
        @show flag_converged, opt_norm

        # Atualizamos o valor de 'x' para a próxima iteração
        x .= xe

        # Calculamos as restrições no valor de 'x' atual
        # g_x(k+1)
        rest_at .= g(x) #[g[j](x) for j=1:m] #Rest_Teste(x)

        # Violação máxima atual
        viola_atual = maximum(rest_at[rest_at.>0.0],init=0.0)
        flag_viola = false
        if  viola_atual > 1.05*max_viola
            println("Violação aumentou de $(max_viola) para $(viola_atual)")
            max_viola = viola_atual
            flag_viola = true
        end

        # Print das restrições calculadas
        println("Restrições ",rest_at[rest_at.>=-1E-3])

        # Atualizamos o termo c_atual
        if flag_viola
            c_atual = min(c_atual * 1.1,c_limite)
        end
     

        # Atualizamos o u_atual para a próxima iteração
        u_atual .= OP.(u_atual .+ (c_atual * rest_at))
        
        
    end # Loop

    return x, u_atual, c_atual, rest_at
end #LA