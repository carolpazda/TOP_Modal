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
