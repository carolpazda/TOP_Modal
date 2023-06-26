#
# Rotina muito simples que retorna os esforços normais da estrutura
#
#     N = sxx * Ae
#
# Os esforços normais podem ser dados por -> N = sxx * A
#
function Esforcos_Normais(mesh::Mesh, x::Vector, Us::Vector,  sparam::Function)

    # Número de elementos
    ne = Get_ne(mesh)

    # Aloca o vetor com os esforços internos da malha
    N = zeros(ne)

    # Retorna os esforços normais pelo Stress_truss
    for ele in mesh

        # Retorna o sxx -> Levando em conta o etype: 3D ou 2D
        sxx0 = Stress(mesh, ele, Us)[1]

        # Retorna o A do elemento
        Ae = Get_geometry(mesh,ele).A

        # Retorna os esforços normais (parametrizados)
        N[ele] = sxx0 * Ae * sparam(x[ele])

    end

    # Retorna os esforços normais
    return N

end