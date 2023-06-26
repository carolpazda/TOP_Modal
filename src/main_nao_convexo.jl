#
# Otimização estrutural utilizando Lagrangiano Aumentado
# TCC: Verônica Caroline Herbst Pazda
#
#  PROBLEMA 1 DO TCC
#
#  Min V = L * A * x
#  x∈R^n
#  T.q
#     (K-w²M)φ = 0
#      -||W||min/w + 1.0 ≤ 0
#      xinf ≤ x ≤ xsup 
#     
#     De forma que,
#
#     ||W||min = (∑i=1^n (w_i^-P))^(-1/P) -> Norma
#
#
# Pacotes que estamos utilizando
using LinearAlgebra, SparseArrays, Plots, ProgressMeter, StaticArrays
using WallE, BMesh, TMeshes, LMesh, LFEM, LinearSolve, LFilter
using DelimitedFiles

# Rotinas carregadas para o main
include("Malha-tower.jl")
include("Malha-bracket.jl")
include("Malha-bracket3D.jl")
include("Malha-coluna2D.jl")
include("Malhas-TCC.jl")
include("Malha-bridge.jl")
include("objetivo_volume.jl")
include("mapeamento.jl")
include("Driver_global-Problema1.jl")
include("derivada_norma_restricao.jl")
include("restricaoNorma_Torii.jl")
include("La.jl")
include("dfc.jl")
include("Loop_LA.jl")
include("valida_dfc.jl")

#
# Programa principal do Problema 1
#
# Recebe um vetor de variáveis de projeto e um fator de 
# correção da função Heaviside aproximada
#
function main(x0::Vector{Float64},verifica_derivada=false)

    println("Entrando em main")


    # Parâmetros que vamos utilizar na projeção
    β = 1.0
    η = 0.5
    ρ_min = 1E-3

    # Define o expoente de penalização da nossa regra de misturas
    penal = 3.0

    # Expoente da norma P (Torii)
    pd = 6.0
    
    # Verifica consistência do expoente da norma
    pd >= 2.0 || throw("Expoente da norma deve ser maior ou igual a 2")
           
    # Número de autovalores que vamos utilizar para calcular a norma
    x_autovalores = 1

    # Número de iterações externas no LA
    NITER = 10

    # Carrega as parametrizações materiais
    #------------------------------------------------------------------------------------------------------------------------------------------#

    # Recebe a parametrização da rigidez e sua derivada
    kparam(xe::Float64,p=penal)=xe^p
    dkparam(xe::Float64,p=penal)=p*xe^(p-1)

    # Recebe a parametrização da massa e sua derivada
    function mparam(xe::Float64,corte=0.1)
        
        saida1 = xe

        rho_min = 1E-3
        C1 = -(6*rho_min - 6)/(corte^5)
        C2 = (5*rho_min - 5)/(corte^6)

        if xe < corte
            saida1 = C1*xe^6  + C2*xe^7
        end

        return saida1

    end
    function dmparam(xe::Float64,corte=0.1)
       
        saida1 = 1.0
        rho_min = 1E-3
        C1 = -(6*rho_min - 6)/(corte^5)
        C2 = (5*rho_min - 5)/(corte^6)

        if xe < corte
            saida1 = 6*C1*xe^5  + 7*C2*xe^6
        end

        return saida1
        
    end

    #------------------------------------------------------------------------------------------------------------------------------------------#
    # Aqui vamos carregar as malhas!
    #------------------------------------------------------------------------------------------------------------------------------------------#
    mesh = Malha3()
    #------------------------------------------------------------------------------------------------------------------------------------------#
    # Vamos preparar a nossa visualização do Gmsh!
    saida = "saida_modos.pos"
    Gmsh_init(saida,mesh)

    # Número de modos para calcular (e visualizar)
    ngls = Get_nn(mesh)*Get_dim(mesh) - (4+1)*2

    nmodos = min(ngls,15)

    # Número de modos para armazenar (para pós processamento)
    NA =  min(nmodos,10)

    massa_full = 7850*Volume(mesh,ones(Get_ne(mesh)))

    @show massa_full

    # Adiciona massas em cima da estrutura
    mesh.options[:Mass] = [77 1 3e3
                           77 2 3e3]
    # Alias para ne
    ne = Get_ne(mesh)

    # Se o vetor de variáveis de projeto for vazio, inicializamos com ums
    if isempty(x0) || length(x0)!=ne
        println("Redimensionando o vetor x0")
        resize!(x0,ne)


        fill!(x0,0.5)

        #x0 = max.(0.0,min.(1.0,ones(ne) .+ 0.01*randn(ne)))

        writedlm("variaveis_teste.txt",x0)

    end

    # Calcula o volume total da estrutura
    full_volume = Volume(mesh,ones(ne))
    @show full_volume

    # Com isso podemos calcular o valor limite de volume
    volume_limite = full_volume

    # Número de variáveis de projeto (primais)
    n = ne
   
    # Número de restrições -> Nesse caso temos apenas uma restrição: Norma!
    d = 1

    # Define os demais coeficientes do problema -> c e u
    c_atual = 10.0
    c_limite = 100.0
    u_atual = zeros(d)

    # Normalização - não estou utilizando no momento, mas vou deixar para teste
    VALS = [0.0]
    
    # Incializa w_esp para poder definir as funções. O valor certo
    # vai ser calculado depois
    w_esp = 1.0

    ###############################################################################################################################################
    # Construi um driver separado para isso -> LA -> DriverLA
    # Agora vamos particularizar o nosso driver para os cálculos que iremos fazer ao longo da solução
    La(x) = DriverGlobal(x, c_atual,  u_atual, mesh,  kparam, dkparam, mparam, dmparam,  penal,  pd,  x_autovalores, nmodos, volume_limite, w_esp, β, η, ρ_min, VALS,"fun_LA")
    modal(x) = DriverGlobal(x, c_atual,  u_atual, mesh,  kparam, dkparam, mparam, dmparam,  penal,  pd,  x_autovalores, nmodos, volume_limite, w_esp, β, η, ρ_min, VALS,"modal")
    g(x) = DriverGlobal(x, c_atual,  u_atual, mesh,  kparam, dkparam, mparam, dmparam,  penal,  pd,  x_autovalores, nmodos, volume_limite, w_esp, β, η, ρ_min, VALS,"restricoes")
    dLa(x) = DriverGlobal(x, c_atual,  u_atual, mesh,  kparam, dkparam, mparam, dmparam,  penal,  pd,  x_autovalores, nmodos, volume_limite, w_esp, β, η, ρ_min, VALS,"derivada")

    # Modos de referência para o problema
    omega,Phi = modal(x0)

    @show omega

    # Frequência especificada pelo usuário
    w_esp = 0.9*omega[1]

    @show w_esp
    @show omega[1]

    # Testa a chamada da função objetivo
    objetivo = La(x0)
    
    @show objetivo, omega[1], g(x0)

    # valida derivada na entrada
    #
    if verifica_derivada
        X0 = copy(x0)
        M0 = Valida_dfc(ne,X0,La,dLa)
        return M0
    end
     
    # Chama o cálculo pelo LA -> Rotina criada: CalculosLA.jl | Cálculo utilizando o WallE
    # Restrições laterais do problema de otimização                
    ci = zeros(length(x0))
    cs = ones(length(x0))
    
    # Armazena os modos
    modos = zeros(NITER+1,NA)
    norma = zeros(NITER+1)
    vol = zeros(NITER+1)

    # Frequências  iniciais da estrutura
    omega, PHI = modal(x0)
    modos[1,:] .= omega[1:NA]

    # Define fora do loop para poder devolver depois
    rest_at = zeros(length(u_atual))

    # Restrições no começo do problema
    rest_at .= g(x0)

    # Restrição
    ρ = x2proj(x0, β, η, ρ_min)
    norma[1] = rest_at[1]

    # Para vermos as variáveis de projeto
    Gmsh_element_scalar(mesh,ρ,saida,"Initial Design variables")


    # Para vermos os modos iniciais da estrutura
    for i=1:nmodos
        Gmsh_nodal_vector(mesh,Phi[:,i],saida,"Modo $i $(omega[i]/(2*pi))[Hz]")
    end


    # Calcula o volume inicial para normalização
    volum_inicia_normalizacao = Volume(mesh,ρ)

    vol[1] =  1.0

    # Violação máxima no começo do problema
    max_viola = maximum(rest_at[rest_at.>0.0],init=0.0)

    # Loop externo do LA
    for k=1:NITER

        # Usando o WallE (Steepest com Bloqueio) -> Pacote carregado do GitHub
        options = WallE.Init()
        options["NITER"] = 1000
        options["LS_ALPHA_INI"] = 10.0
        output = WallE.Solve(La,dLa,x0,ci,cs,options)
        x0 .= output["RESULT"]
        flag_converged = output["CONVERGED"]
        opt_norm = output["NORM"]
        @show flag_converged, opt_norm

        
        # Aqui temos que mapear, pois não estamos usando o driver
        ρ = x2proj(x0, β, η, ρ_min)
        println("Valor máximo das variáveis projetadas ", maximum(ρ))

        # Analise modal da solução atual
        omega, _ = modal(x0)
        modos[k+1,:] .= omega[1:NA]

        # Calculamos as restrições no valor de 'x' atual
        # g_x(k+1)
        rest_at .= g(x0) #[g[j](x) for j=1:m] #Rest_Teste(x)

        norma[k+1] = rest_at[1]
        vol[k+1] = Volume(mesh,ρ)/volum_inicia_normalizacao

        # Violação máxima atual
        viola_atual = maximum(rest_at[rest_at.>0.0],init=0.0)
        flag_viola = false
        if  viola_atual > max_viola
            println("Violação aumentou de $(max_viola) para $(viola_atual)")
            max_viola = viola_atual
            flag_viola = true
        end

        # Print das restrições calculadas
        writedlm("variaveis_projeto_($k).txt", ρ)

        # Atualizamos o termo c_atual
        if flag_viola
            c_atual = min(c_atual * 1.1,c_limite)
        end

        # Mostra quantas variáveis de projeto são intermediárias
        println("Variáveis de projeto intermediárias ",length([1.1E-3.>=ρ.<=0.9])," de ",ne)

        # Atualizamos o u_atual para a próxima iteração
        u_atual .= OP.(u_atual .+ (c_atual * rest_at))
        
        # Incrementa o β na maciota
        if k>5
            β += 1.0
        end

    end # Loop

    β -= 1.0

    # Aqui temos que mapear, pois não estamos usando o driver
    ρ = x2proj(x0, β, η, ρ_min)
 
    # Para vermos as variáveis de projeto
    Gmsh_element_scalar(mesh,ρ,saida,"Design variables")

    # Para vermos os modos da estrutura
    for i=1:nmodos
        Gmsh_nodal_vector(mesh,Phi[:,i],saida,"Modo $i $(omega[i]/(2*pi))[Hz]")
    end

    # Mostra quantas variáveis de projeto são intermediárias
    println("Variáveis de projeto intermediárias ",length(ρ[1.1E-3.>=ρ.<=0.9])," de ",ne)
    println("Valor máximo das variáveis de projeto ", maximum(ρ))

    @show rest_at[end]
        
    writedlm("variaveis_projeto_final.txt",ρ)
    writedlm("modos.txt",modos)

    # Retorna o vetor de variáveis de projeto de interesse
    if verifica_derivada
       return M0, G0, X0
    else
        @show w_esp
        return modos, norma, vol, x0
    end
        
end #main
