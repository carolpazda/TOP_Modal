#
# DRIVER APENAS DA RESTRIÇÃO GLOBAL - TESTE
#

function DriverGlobal(x::Vector, c_atual::Float64,  u_atual::Vector{Float64}, mesh::Mesh, 
                     kparam::Function, dkparam::Function, mparam::Function, dmparam::Function,
                     penal::Float64, pd::Float64,  x_autovalores::Int64, nmodos::Int64,
                     volume_limite::Float64, w_esp, β::Float64, η::Float64, ρ_min::Float64, 
                     VALS::Vector{Float64}, tipo::String)


    # Tipos válidos
    tipo in ["fun_LA", "restricoes", "modal", "derivada"] || throw("Driver_global:: tipo $tipo não é válido")
    
    # Mapeamento x -> rho
    ρ = x2proj(x, β, η, ρ_min)
   
    # Número de variáveis de projeto
    ne = length(ρ)

    # Dimensão
    dim = Get_dim(mesh)

    # Número de graus de liberdade total 
    ngls = dim*Get_nn(mesh)

    # Graus de liberdade livres do problema
    free_dofs = mesh.free_dofs[1]

    # Matriz de rigidez global
    K = Global_K(mesh,ρ,kparam)

    # Matriz de massa global
    M = Global_M(mesh,ρ,mparam)
    
    #
    # Soluciona o problema de autovalores e autovetores
    #
    @assert nmodos >= 2*x_autovalores "Calcular pelo menos o dobro dos autovalores que serão utilizados na norma"
   
    # Estou devolvendo um flag no Solve_Eigen, de modo que
    # Flag =  1 -> Deu tudo certo;
    # Flag = -1 -> Deu pau completamente, irrecuperável :c ;
    # Flag = -2 -> Não está satisfazendo a equação de equilíbrio, pânico no parquinho!;
    # Flag = -3 -> Os autovetores não são ortogonais - nada de mais.   
    flag, lamb, phi =  Solve_Eigen_(K[free_dofs,free_dofs],M[free_dofs,free_dofs],nmodos,positive=true)#,verbose=true)
    #@show flag

    # Número efetivo de autovalores calculados
    n_efetivo = length(lamb)
    n_efetivo >= x_autovalores || error("Driver_global::não temos um número suficiente de autovalores para calcular a norma ", n_efetivo) 

    # Aumenta o PHI
    PHI = zeros(ngls,size(phi,2))
    for i=1:size(phi,2)
        vtemp = zeros(ngls)
        Expand_vector!(vtemp,phi[:,i],free_dofs)
        PHI[:,i] .= vtemp
    end #i

    # Verifica se estamos satisfazendo o problema modal
    #for i=1:x_autovalores
    #    norm(Residue_Eigenpair(K[free_dofs,free_dofs], M[free_dofs,free_dofs], PHI[free_dofs,i], lamb[i]))<1E-0 || throw("Deu caca")
    #    @show norm(Residue_Eigenpair(K[free_dofs,free_dofs], M[free_dofs,free_dofs], PHI[free_dofs,i], lamb[i]))
    #end

    # Mais uma verificação, só para garantir
    #for i=1:x_autovalores
    #   norm(dot(PHI[:,i], Residue_Eigenpair(K, M, PHI[:,i], lamb[i])))<1E-1 || throw("Deu caca 2")
    #end

    # Converte os lambdas para frequência angular 
    omega = sqrt.(lamb)

    # Devolve as frequências naturais angulares (rad/s) e os modos
    if tipo=="modal"
        return omega, PHI
    end

    # Calcula o miolo da Norma Torii
    mioloNorma = Miolo_NormTorii(pd, omega, x_autovalores)

    # Calcula a função objetivo - volume
    obj_volume =  Volume(mesh,ρ)/volume_limite

    # Calculamos a primeira parte da restrição da norma
    norma_torii = Restricao_norma_Torii(mioloNorma, pd)
    
    # Normalização do objetivo
    #if sum(VALS)==0

    #   VALS[1] = abs(norma_torii) 

    #end

    # Normaliza
    #
    # w >= w_lim
    # w/w_lim  >= 1
    # - w/w_lim + 1 <=0 
    #
    rest_norma = -(norma_torii / w_esp) + 1.0

    # Concatena as restrições -> Nesse primeiro caso temos apenas uma - de frequência
    gla = [rest_norma]

    # Se o driver é para restrições, devolve os valores
    if tipo=="restricoes"
       return gla         
    end
    
    # Calcula a função Lagrangiano aumentado
    LA = ObjetivoLa(obj_volume, gla, c_atual, u_atual)
    
    if tipo=="fun_LA"
        return LA
    end

    #
    #                                       DERIVADAS
    #
    # A derivada da nossa restrição vai ser o seguinte: 
    #       g1 = -||w||/w + 1.0
    # dg11/dxm = -1/w d||w||/dxm
    #
    D_f = Derivada_norma_Torii( pd, omega, x_autovalores, mesh, ρ, PHI, dkparam, M, dmparam, free_dofs)
    #========================== DERIVADA DA RESTRIÇÃO  ===========================#
    gdf =  c_atual*OP(u_atual[1]/c_atual + rest_norma)*((-1/w_esp)*(D_f))

    #========================== DERIVADA DO VOLUME ===========================#

    # Derivada do volume -> Que é a nossa função objetivo
    gdV = dVolume(mesh,ρ) ./ volume_limite

    # Assim...o gradiente do LA em relação a x é (trocando o sinal do objetivo, pois é maximização)
    dLa = gdV .+ gdf 
  
    # Corrige a derivada 
    dproj2dx!(x,dLa,β, η, ρ_min)
    
    return dLa
end