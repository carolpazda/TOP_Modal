#
#
#
# Organize the eigenvalues and eigenvector . To be used as a post-processor
# for the Modal Analysis.
#
#
#
function Organize_Eigen(lambda::Vector,phi::Matrix,ngls::Int64,free_dofs::Vector)

    # Convert to real numbers
    lamb_before = real.(lambda)

    # Number of eigenvalues
    nev = length(lamb_before)

    # Make sure to get only the positive eigenvalues
    # in crescent order
    n_effective = 0
    lamb_a  = Float64[]
    pos_lamb = Int64[]
    for i=1:nev
        if lamb_before[i] > 0.0
            push!(lamb_a,lamb_before[i])
            push!(pos_lamb,i) 
            n_effective += 1
        end
    end

    # Avoid the situation of no positive eigenvalue
    n_effective >=1 || error("Organize_Eigen:: there is no valid positive solution")

    # sort
    ordem = sortperm(lamb_a)
    lamb = lamb_a[ordem]

    # Convert the eigenvectors to real numbers
    phi_real = real.(phi[:,pos_lamb[ordem]])

    # Alocate the matrix (full number of dofs)
    PHI = zeros(ngls,n_effective)

    # Expand the eigenvectores
    for i=1:n_effective

        # Expande esse modo para os gls globais
        Usf  = zeros(ngls)
        Expand_vector!(Usf,real.(phi_real[:,i]),free_dofs)
        PHI[:,i] .= Usf

    end

    # Return the positive eigenvalues and their eigenvectors
    return lamb, PHI

end

