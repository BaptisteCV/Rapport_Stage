# Conditions de bord : formules de mise à jour des coefficients

**Configuration testée** : `input/testDFN/BoundaryConditions.txt`
```
LEFT:    MIXED_BC
RIGHT:   MIXED_BC  
BOTTOM:  MIXED_BC
TOP:     NEUMANN
```

---

## 1) Matrice poreuse (EPM)

Pour chaque cellule de bord `(i,j)`, la fonction `BoundCondCoeffs(...)` retourne `(coeff1, coeff3)` qui mettent à jour le système linéaire selon :

$$A_{cc} \leftarrow A_{cc} + \texttt{coeff1}, \quad b_c \leftarrow b_c + \texttt{coeff3}$$

avec $c = nb\_nodes\_fract + \text{index}(i,j)$.

### Facteur géométrique

$$F(border) = \begin{cases}
\dfrac{\Delta y}{\Delta x} & \text{if } border \in \{\text{LEFT, RIGHT}\} \\
\dfrac{\Delta x}{\Delta y} & \text{if } border \in \{\text{BOTTOM, TOP}\}
\end{cases}$$

### MIXED_BC (LEFT, RIGHT, BOTTOM)

$$\texttt{coeff1} = \begin{cases}
\dfrac{2\,\sigma\,\alpha_{bc}\,\Delta y}{2 + \alpha_{bc}\,\Delta x} & \text{LEFT/RIGHT} \\
\dfrac{2\,\sigma\,\alpha_{bc}\,\Delta x}{2 + \alpha_{bc}\,\Delta y} & \text{BOTTOM}
\end{cases}$$

$$\texttt{coeff3} = 0$$

où $\sigma$ est la conductivité et $\alpha_{bc}$ est calculé par `ReturnMixedBoundCondValue2(...)` en Fourier :

$$\alpha_{bc} = \omega \frac{K_1(\omega r)}{K_0(\omega r)} \times \begin{cases}
\dfrac{|x_b - x_s|}{r} & \text{LEFT/RIGHT} \\
\dfrac{|y_b - y_s|}{r} & \text{BOTTOM}
\end{cases}$$

avec $K_0, K_1$ fonctions de Bessel modifiées de seconde espèce, $x_s, y_s$ position source/électrode, $r$ distance source–point bord.

### NEUMANN (TOP)

$$\texttt{coeff1} = 0$$

$$\texttt{coeff3} = \sigma \, F(\text{TOP}) \, q_N = \sigma \frac{\Delta x}{\Delta y} \, q_N$$

---

## 2) Nœuds de fracture à la frontière

### MIXED_BC (LEFT, RIGHT, BOTTOM)

La fonction `BCDDPLinearSystemFourier(...)` ajoute à la diagonale du nœud $n$ :

$$A_{nn} \leftarrow A_{nn} + \frac{2\,e_f\,\sigma_f\,\alpha_n}{2 + \alpha_n\,\Delta x}$$

où :
- $e_f$ = ouverture de la fracture au nœud $n$
- $\sigma_f$ = conductivité de la fracture au nœud $n$  
- $\alpha_n$ calculé par `ReturnAlphaBCMixed(...)` :

$$\alpha_n = \omega \frac{K_1(\omega r)}{K_0(\omega r)} \times \begin{cases}
\dfrac{|x_n - x_s|}{r} & \text{LEFT/RIGHT} \\
\dfrac{|y_n - y_s|}{r} & \text{BOTTOM}
\end{cases}$$

### NEUMANN (TOP)

$$b_n \leftarrow b_n + q_N$$

---

## 3) Résumé des formules par région

| Type | Région | Terme diagonal ($A_{cc}$) | Terme RHS ($b_c$) |
|------|--------|---------------------------|-------------------|
| **MIXED_BC** | LEFT/RIGHT (poreuse) | $\dfrac{2\sigma\alpha\Delta y}{2+\alpha\Delta x}$ | $0$ |
| **MIXED_BC** | BOTTOM (poreuse) | $\dfrac{2\sigma\alpha\Delta x}{2+\alpha\Delta y}$ | $0$ |
| **MIXED_BC** | LEFT/RIGHT/BOTTOM (fracture) | $+\dfrac{2e_f\sigma_f\alpha_n}{2+\alpha_n\Delta x}$ | — |
| **NEUMANN** | TOP (poreuse) | $0$ | $\sigma\Delta x/\Delta y \cdot q_N$ |
| **NEUMANN** | TOP (fracture) | — | $+q_N$ |

---

