<img alt="GitHub Language Count" src="https://img.shields.io/github/languages/count/carolpazda/TOP_Modal" />

<img alt="GitHub Issues" src="https://img.shields.io/github/issues/carolpazda/TOP_Modal" />

<img alt="GitHub Contributors" src="https://img.shields.io/github/contributors/carolpazda/TOP_Modal" />

<img alt="GitHub Last Commit" src="https://img.shields.io/github/last-commit/carolpazda/TOP_Modal" />

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

### Saída de dados e visualização
O presente código computacional retorna os valores dos modos, norma, volume e das variáveis de projeto,
```
return modos, norma, vol, x0
```
Além disso, retorna um arquivo de saída para o [Gmsh](https://gmsh.info/), o arquivo é gravado como _saida_modos.pos_. Te juro, o gmsh é muito legal para realizar a visualização das malhas!!! :speak_no_evil:. Vamos ver algumas dicas de como você pode visualizar suas topologias no Gmsh?<br/>
Ao entrarmos no Gmsh, você pode abrir o arquivo clicando em _open_ e desabilitar todas as visualizações, para conseguirmos visualizar uma por vez. 

![gmsh-1](https://github.com/carolpazda/TOP_Modal/assets/107930972/e47a6659-0034-4b2f-bddf-b0348cf7d056)

Ao clicarmos em uma visualização, geralmente não vamos ter muitas informações da malha. Então, você pode clicar em _options_ para entrar nas configurações dessa visualização. Altere os parâmetros da visualização e puff... agora dá pra ver tudo, né? :eyes:

![gmsh-2](https://github.com/carolpazda/TOP_Modal/assets/107930972/0fd252e8-fc83-4076-8615-d148bb1abe7a)

Você ainda pode alterar a cor das barras com o _Map_ e visualizar os modos da estrutura com o _Transfo_. E olha só que legal, estamos exportando as topologias e os modos de todas as iterações, então podemos ver de que forma a topologia está mudando ao longo das iterações! :ok_hand:

![gmsh-3](https://github.com/carolpazda/TOP_Modal/assets/107930972/2c83600d-4f00-4462-8af6-eec4dfa92aa5)


Ficou alguma dúvida? Eu fiz um vídeo no YouTube explicando um pouco mais sobre o código e mostrando como você pode utilizá-lo, você pode dar uma conferida em :point_right: *Colocar assim que estiver pronto* :point_left: <br/>

Tu viu o vídeo e ficou com dúvidas ainda :eye_speech_bubble:? Então me contate pelo e-mail que deixei lá no início :point_up_2:! Grandes abraços! :wave: :wave:

