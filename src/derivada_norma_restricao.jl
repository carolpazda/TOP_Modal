#
# Derivada da Norma
#
#
function Derivada_norma_Torii( pd::Float64, omega::Vector, x_autovalores::Int64,
                              mesh::Mesh, x::Vector, PHI::Matrix, dkparam::Function,
                              M::AbstractMatrix{Float64}, dmparam::Function,free_dofs)

    
    # Número de elementos
    ne = Get_ne(mesh)

    # Aloca o vetor de saída
    D = zeros(ne)

    # Vamos calcular o termo T_min 
    t_1 = Miolo_NormTorii(pd, omega, x_autovalores)
    T_min = t_1^((-1/pd) - 1)

    # Calcula os denominadores antes do loop 
    denominador = [2*omega[i]*dot(PHI[free_dofs,i],M[free_dofs,free_dofs],PHI[free_dofs,i]) for i=1:x_autovalores]

    # Loop externo
    for ele in mesh

        # Variáveis de projeto
        xe = x[ele]

        # Derivadas das parametrizações
        Dmparam = dmparam(xe)
        Dkparam = dkparam(xe)
        
        # Matriz local de massa do elemento
        m = Local_M(mesh,ele)

        # Matriz local de rigidez do elemento
        k = Local_K(mesh,ele)

        # Matriz de rotação
        R = T_matrix(mesh.bmesh,ele)

        # Derivada da matriz local de rigidez do elemento
        # o sistema global de referência
        dkl = Dkparam * transpose(R)*k*R #   To_global(k,mesh,ele) #dot(R,k,R)

        # Derivada da matriz local de massa do elemento
        # o sistema global de referência
        dml = Dmparam * transpose(R)*m*R #To_global(m,mesh,ele) #dot(R,m,R)

        # Gls do elemento
        gls = DOFs(mesh,ele)
        
        # Loop interno -> autovalores
        for i=1:x_autovalores

            # Autovetor i do elemento ele
            phi = PHI[gls,i]

            # Primeiro termo_a
            t1 = omega[i]^(-pd-1)            
 
            # Termo de cima da derivada do omega
            numerador = dot(phi,dkl.-dml*omega[i]^2,phi)

            # Derivada do omega_i em relação ao x_ele
            domega_i_ele = numerador/denominador[i]

            # Somatório do autovalor i na posição ele
            D[ele] += T_min*t1*domega_i_ele

        end # Loop interno (autovalor)

    end # Loop externo (elementos)

    # Retorna o valor de D
    return D

end # Function