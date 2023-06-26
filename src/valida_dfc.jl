function Valida_dfc(n, x0,  La::Function, dLa::Function, xin=0.0, delta=1E-6)

# Validação da derivada do LA por meio de diferenças finitas
    
        # Limites laterais para as variáveis de projeto
        xmin = zeros(n)
        xmax = ones(n)

        # Derivada analítica
        println("Analitico")
        referenciad_D = dLa(x0) 

        # Calculando por DFC
        println("DFC ", delta)
        dfD = DFC(La, x0, xmin, xmax, delta)
        println("FINAL DFC ")
    
        # Calculando por DFC
        println("DFC ", 1E-8)
        dfD8 = DFC(La, x0, xmin, xmax, 1E-8)
        println("FINAL DFC ")


        return [referenciad_D  dfD  dfD8]

end