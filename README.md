# TOP_Modal

![Capa](https://github.com/carolpazda/TOP_Modal/assets/107930972/b1e196a8-e4e5-4d47-882b-45551257e5bc)

_Trabalho de Conclusão de Curso - Bacharelado em Engenharia Mecânica_.

_Autoria: Verônica Caroline Herbst Pazda e Eduardo Lenz Cardoso_

_Email para contato_: veronicacpazda@outlook.com

# CÓDIGO COMPUTACIONAL

## O Problema
O objetivo do presente trabalho consiste no estudo da otimização topológica de treliças 2D e 3D considerando a resposta modal, trabalhando-se com um problema de minimização do volume com restrição de frequência (imposição de um valor mínimo para o menor autovalor). Para tal, o problema de otimização é definido por

```math
       \mathbf{P}_1 \left \lbrace \begin{array}{ccc}
       Min  &  V(\mathbf{{x}})&\\
        T.q. &  \left( \mathbf{K}-\omega^2 \mathbf{M} \right){\varphi} = \mathbf{0} & \\
        & -\frac{||{\omega}||_{min}}{\alpha \omega_{\text{min}}^{\text{full}}} + 1 \leq 0 & \\
        & \underline{\mathbf{x}} \leq \mathbf{x} \leq \bar{\mathbf{x}} & 
    \end{array} \right.
```
onde,
```math
  ||{\omega}||_{min} = \left( \sum_{i=1}^n \omega^{-P_d}_i  \right)^{-\frac{1}{P_d}}
```
é a norma que considera o menor valor entre as $n$ primeiras frequências naturais da estrutura (armazenadas em ordem crescente em um vetor $\omega$). O volume da estrutura, função objetivo do presente problema, é dado pela soma das contribuição de cada uma das barras

```math
       V\left(\mathbf{x}\right)=\sum \limits_{i=1}^{n} \rho_i(\mathbf{x}) A_i L_i.
```

## Entrada de dados - Código computacional
Os pacotes utilizados são:

using [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/), [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/), [Plots](https://docs.juliaplots.org/stable/), [ProgressMeter](https://docs.juliahub.com/ProgressMeter/3V8n6/1.3.1/), [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) <br/>
using [WallE](https://github.com/CodeLenz/WallE.jl), [BMesh](https://github.com/CodeLenz/BMesh.jl), [TMeshes](https://github.com/CodeLenz/TMeshes.jl), [LMesh](https://github.com/CodeLenz/LMesh.jl), [LFEM](https://github.com/CodeLenz/LFEM.jl), [LinearSolve](https://github.com/SciML/LinearSolve.jl), [LFilter](https://github.com/CodeLenz/LFilter.jl)<br/>
using [DelimitedFiles](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/)

### Utilização do _main_
Para entrar no _main_ do código você deve entrar com algum vetor _Float64_ inicial para as variáveis de projeto. O código posteriormente realiza a correção dessas varíaveis de entrada, então não se preocupe com esse valor :grin:!

```
main(x0::Vector{Float64},verifica_derivada=false)
```

Após isso, você pode arbitrar para o código os demais parâmetros. Mas cuidado aí em, você pode arbitrar, mas tenha coerência com seus valores :upside_down_face:.

```
# Parâmetros que vamos utilizar na projeção
β = 5.0
η = 0.5
ρ_min = 1E-3

# Define o expoente de penalização da nossa regra de misturas
penal = 3.0

# Expoente da norma P (Torii)
pd = 6.0
    
# Verifica consistência do expoente da norma
pd >= 2.0 || throw("Expoente da norma deve ser maior ou igual a 2")
           
# Número de autovalores que vamos utilizar para calcular a norma
x_autovalores = 6

# Número de iterações externas no LA
NITER = 10
```

Você pode inserir a malha de interesse no _main_ da seguinte forma

```
mesh = Tower3D()
```
Existem diversas malhas que podem ser utilizadas para teste, é só alterar essa entrada no _main_ com o nome da sua malha de interesse. Não esqueça que você também pode criar uma outra malha (caso for de seu interesse :blush:).<br/>

Para inserção de massas concentradas você deve inserir primeiro o nó, depois a direção (sendo 1 -> x, 2 -> y e 3 -> z) e após isso o valor da massa utilizada, conforme o exemplo

```
# Adiciona massas em cima da estrutura
mesh.options[:Mass] = [55 1 1e3
                       55 2 1e3
                       55 3 1e3
                       57 1 1e3
                       57 2 1e3
                       57 3 1e3
                       61 1 1e3
                       61 2 1e3
                       61 3 1e3
                       63 1 1e3
                       63 2 1e3
                       63 3 1e3]
```
Com relação à restrição, você pode modificar o $\alpha$ que é utilizado na restrição alterando a entrada
```
# Frequência especificada pelo usuário
w_esp = 0.9*omega[1]
```

