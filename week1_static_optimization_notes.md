# Week 1 Lecture Notes — Chapter 1: Static Optimization

**Course:** Optimal Control  
**Source:** Frank L. Lewis, Draguna L. Vrabie, and Vassilis L. Syrmos, *Optimal Control*, 3rd ed., Chapter 1 in `Lewis Ch1_ 2_ 6_ 11 with annotation.pdf`, PDF P10–P27.  
**Scope:** Sections 1.1–1.3 (PDF P10–P24), with selected practice based on PDF P25–P27. Page labels below refer to **PDF page numbers**, not the book's printed page numbers.

## Contents

1. [The problem and notation](#1-the-problem-and-notation)
2. [Unconstrained optimization](#2-unconstrained-optimization)
3. [Equality-constrained optimization](#3-equality-constrained-optimization)
4. [Numerical solution methods](#4-numerical-solution-methods)
5. [Connection to optimal control](#5-connection-to-optimal-control)
6. [One-page summary](#6-one-page-summary)
7. [Glossary](#7-glossary)
8. [Practice problems](#8-practice-problems)
9. [Answers to practice problems](#9-answers-to-practice-problems)

## 1. The problem and notation

**PDF P10.** Static optimization asks us to choose a **decision vector** $u\in\mathbb{R}^m$ that minimizes a scalar **performance index** $L(u)$, without time as a variable:

$$
\min_{u\in\mathbb{R}^m} L(u).
$$

For example, $u$ could be several converter settings, while $L$ measures a combination of voltage error and control effort. A smaller $L$ represents a better choice under the chosen engineering objective. **The objective must be specified first:** mathematics can minimize $L$, but it cannot decide whether $L$ captures what the engineer actually cares about.（先定義「好」的量化標準，再談最佳化。）

Chapter 1 builds the static tools used later when states and controls change over time. The notation used throughout these notes is:

| Symbol | Shape | Meaning |
|---|---:|---|
| $u$ | $m\times1$ | decision or control vector |
| $x$ | $n\times1$ | auxiliary/state vector in a constrained problem |
| $L$ | scalar | performance index |
| $L_u=\partial L/\partial u$ | $m\times1$ | gradient, defined as a **column** vector |
| $L_{uu}$ | $m\times m$ | Hessian or curvature matrix |
| $f(x,u)$ | $n\times1$ | equality constraint |
| $f_x$, $f_u$ | $n\times n$, $n\times m$ | constraint Jacobians |
| $\lambda$ | $n\times1$ | Lagrange multiplier |

An increment $du$ means a small **change** in $u$. The new point is $u+du$; these three expressions are different objects. Likewise, $x+dx$ is the new state, while $dx$ is only its change.（$dx$ 不是另一個狀態，而是狀態的微小變化。）

## 2. Unconstrained optimization

### 2.1 Gradient and Hessian

**PDF P10–P11, Section 1.1.** Near a chosen point, Taylor expansion describes how the cost changes:

$$
L(u+du)-L(u)
=dL
=L_u^Tdu+\frac12du^TL_{uu}du+O(\|du\|^3).
$$

The first term describes the local **slope**; the second describes **curvature**. For $u=(u_1,\ldots,u_m)^T$,

$$
L_u=\begin{bmatrix}\partial L/\partial u_1\\\vdots\\\partial L/\partial u_m\end{bmatrix},
\qquad
L_{uu}=\begin{bmatrix}
\partial^2L/\partial u_1^2&\cdots&\partial^2L/\partial u_1\partial u_m\\
\vdots&\ddots&\vdots\\
\partial^2L/\partial u_m\partial u_1&\cdots&\partial^2L/\partial u_m^2
\end{bmatrix}.
$$

**Shape check:** $L_u^Tdu$ is $(1\times m)(m\times1)$, hence a scalar; $du^TL_{uu}du$ is also a scalar. The gradient is a vector of partial derivatives, **not their sum**. The Hessian is a whole matrix; zero cross-partials only make its off-diagonal entries zero.（梯度是向量；Hessian 是矩陣，不能把元素直接加起來。）

**Worked example A — one decision variable.** Let $L(u)=(u-4)^2+3$. Then

$$
L_u=2(u-4),\qquad L_{uu}=2.
$$

Setting $L_u=0$ gives $u^*=4$. Since the curvature is positive, this is a minimum, and $L^*=L(4)=3$.

**Worked example B — two decision variables.** For $L(u_1,u_2)=u_1^2+2u_2^2$,

$$
L_u=\begin{bmatrix}2u_1\\4u_2\end{bmatrix},\qquad
L_{uu}=\begin{bmatrix}2&0\\0&4\end{bmatrix}.
$$

The only stationary point is $(0,0)$, and both curvature directions are positive, so it is a strict minimum with $L^*=0$.

### 2.2 What a stationary point tells us

**PDF P11–P12.** At an interior stationary point, every small direction has zero *first-order* cost change:

$$L_u(u^*)=0.$$

This is a **necessary condition** for a smooth interior local minimum. It is not, by itself, proof of a minimum. At $L_u=0$, the quadratic term supplies a useful second-order test:

| Hessian at the stationary point | Local conclusion |
|---|---|
| $L_{uu}\succ0$ (positive definite) | strict local minimum |
| $L_{uu}\prec0$ (negative definite) | strict local maximum |
| $L_{uu}$ indefinite | saddle point |
| $L_{uu}$ semidefinite or singular | second-order test is inconclusive; inspect higher-order terms or the function directly |

Here $L_{uu}\succ0$ means $v^TL_{uu}v>0$ for **every nonzero** direction $v$. A local result concerns nearby points; a global minimum additionally requires comparison with all feasible points.（$L_u=0$ 只表示「坡度暫時為零」；還要看各方向是向上彎、向下彎，或兩者都有。）

**Worked example C — saddle.** If $L(u_1,u_2)=u_1^2-u_2^2$, then

$$
L_u=\begin{bmatrix}2u_1\\-2u_2\end{bmatrix},\qquad
L_{uu}=\begin{bmatrix}2&0\\0&-2\end{bmatrix}.
$$

The origin is stationary. Moving along $u_1$ raises $L$, while moving along $u_2$ lowers it; therefore the origin is a saddle. Notice that the cross-partials are zero, **but the Hessian is not the zero matrix**.

### 2.3 Quadratic costs and contours

**PDF P11–P13, Examples 1.1-1 and 1.1-2.** Let $Q=Q^T\in\mathbb{R}^{m\times m}$ and $S\in\mathbb{R}^m$. For

$$
L(u)=\frac12u^TQu+S^Tu,
$$

matrix differentiation gives

$$
L_u=Qu+S,\qquad L_{uu}=Q.
$$

If $Q$ is invertible, the stationary point and its cost are

$$
u^*=-Q^{-1}S,\qquad L(u^*)=-\frac12S^TQ^{-1}S.
$$

These formulas locate and evaluate a stationary point. Calling it a *minimum* still requires $Q\succ0$. For a symmetric $2\times2$ matrix, $\det Q<0$ means one positive and one negative eigenvalue, hence a saddle. In the textbook, $|Q|$ means **$\det(Q)$**, the determinant, rather than an absolute value. A zero determinant makes $Q^{-1}$ unavailable and the second-order classification may be inconclusive.（先解 $Qu+S=0$，再用 $Q$ 判斷是不是最小值。）

**Worked example D — matrix calculation.** Choose

$$
Q=\begin{bmatrix}2&0\\0&4\end{bmatrix},\qquad
S=\begin{bmatrix}-4\\8\end{bmatrix}.
$$

Then $Qu+S=0$ gives $2u_1-4=0$ and $4u_2+8=0$, so

$$
u^*=\begin{bmatrix}2\\-2\end{bmatrix},\qquad
L^*=\frac12(2\cdot 2^2+4\cdot(-2)^2)-4(2)+8(-2)=-12.
$$

Since $Q\succ0$, this quadratic has a unique global minimum.

The source's two-dimensional example uses $Q=\left[\begin{smallmatrix}1&1\\1&2\end{smallmatrix}\right]$ and $S=(0,1)^T$. Its stationary point is $(1,-1)^T$, with $L^*=-1/2$ (PDF P12–P13). This gives a concrete interpretation of the contour plot: a **contour** is the set of points with one fixed value of $L$. Along a contour, $dL\approx L_u^Tdu=0$, so the gradient is perpendicular to its tangent direction. $L_u$ points toward increasing cost, and $-L_u$ gives the direction of steepest local decrease in the usual Euclidean scale.（等高線上數值不變；梯度垂直等高線，負梯度指向下降最快的方向。）

## 3. Equality-constrained optimization

### 3.1 Feasible changes couple $x$ and $u$

**PDF P13–P15, Section 1.2.** The new problem is

$$
\min_{x,u}L(x,u)\qquad\text{subject to}\qquad f(x,u)=0,
$$

where $x\in\mathbb{R}^n$, $u\in\mathbb{R}^m$, and $f\in\mathbb{R}^n$. A **feasible point** satisfies every component of $f=0$. An infeasible point is not an eligible answer even when its cost is lower.（限制式是入場資格；成本低也不能違反限制。）

The textbook assumes that $f_x$ is nonsingular near the point of interest, so the constraint locally determines $x$ for a given $u$. To stay on the constraint surface after small changes, we need

$$
df=f_xdx+f_udu=0,
\qquad
dx=-f_x^{-1}f_udu.
$$

Dimensions: $f_x^{-1}f_u$ is $(n\times n)(n\times m)=n\times m$, so multiplying by $du\in\mathbb{R}^m$ produces $dx\in\mathbb{R}^n$. Meanwhile,

$$dL=L_x^Tdx+L_u^Tdu.$$

Substituting the feasible $dx$ gives

$$
dL=\left(L_u^T-L_x^Tf_x^{-1}f_u\right)du.
$$

For this to vanish for every small control change $du$, a constrained stationary point must satisfy

$$
\boxed{L_u-f_u^Tf_x^{-T}L_x=0},
\qquad f_x^{-T}:=(f_x^{-1})^T.
$$

This is a **first-order necessary condition**, given the smoothness and nonsingular-$f_x$ assumptions. The direct gradient $L_u$ alone is insufficient because changing $u$ also changes $x$.（調整控制量時，狀態會被限制式牽動；要把兩者的成本變化一起算。）

### 3.2 Lagrange multiplier and Hamiltonian

**PDF P15–P17.** Introduce the **Lagrange multiplier** $\lambda\in\mathbb{R}^n$ and the chapter's **Hamiltonian**

$$
H(x,u,\lambda)=L(x,u)+\lambda^Tf(x,u).
$$

The stationary conditions are

$$
\boxed{H_\lambda=f(x,u)=0},
\qquad
\boxed{H_x=L_x+f_x^T\lambda=0},
\qquad
\boxed{H_u=L_u+f_u^T\lambda=0}.
$$

| Condition | What it does |
|---|---|
| $H_\lambda=0$ | enforces feasibility: the original constraint is satisfied |
| $H_x=0$ | balances the state gradient and determines the multiplier locally |
| $H_u=0$ | makes the cost stationary with respect to the decision/control |

Why does $H_\lambda=f$? $L$ does not depend on $\lambda$, while the gradient of $\lambda^Tf$ with respect to $\lambda$ is exactly $f$. **Keep the original sign of the constraint:** if $f=x+u-2$, then $H=L+\lambda(x+u-2)$ and $H_\lambda=x+u-2$.（對乘子微分，就是把原本的限制式拿回來。）

The multiplier method avoids explicitly substituting $dx=-f_x^{-1}f_udu$ each time. Indeed, $H_x=0$ gives $\lambda=-f_x^{-T}L_x$, and inserting this into $H_u=0$ reproduces the boxed condition in Section 3.1. These equations produce **candidate** optima; a curvature check is still needed.

**Worked example E — one equality constraint.** Minimize $L(x,u)=x^2+u^2$ subject to $x+u-2=0$. Construct

$$H=x^2+u^2+\lambda(x+u-2).$$

The three equations are

$$
H_x=2x+\lambda=0,\qquad
H_u=2u+\lambda=0,\qquad
H_\lambda=x+u-2=0.
$$

The first two give $x=u$; the third then gives $x^*=u^*=1$ and $\lambda^*=-2$. Thus $L^*=2$. Substitution along the feasible line, $x=2-u$, gives

$$L(2-u,u)=(2-u)^2+u^2=2(u-1)^2+2,$$

which confirms a unique global constrained minimum. The unconstrained point $(0,0)$ has a lower cost, $0$, but is **infeasible** because $0+0-2\ne0$.

### 3.3 Geometry, curvature, and sensitivity

**PDF P17–P23.** Write $z=(x^T,u^T)^T$. The stationary equations combine into

$$\nabla_zL+J_f^T\lambda=0.$$

Thus the objective gradient is balanced by a combination of constraint gradients. With one smooth equality constraint, the objective contour and constraint curve touch at the stationary point; the gradients are parallel. Feasible motion is tangent to the constraint, so the objective has zero first-order change in those tangent directions.（沿限制線走，最優點附近的一階成本變化為零。）

For a **sufficient local-minimum test**, check second-order curvature only along feasible directions. If $f_x$ is invertible, define

$$
T=\begin{bmatrix}-f_x^{-1}f_u\\I_m\end{bmatrix},
\qquad
\begin{bmatrix}dx\\du\end{bmatrix}=Tdu
\quad\text{to first order when }df=0.
$$

The relevant constrained curvature is

$$
L_{uu}^{f}=T^TH_{zz}T,
\qquad H_{zz}=\nabla_z^2H.
$$

If $L_{uu}^{f}\succ0$ at a stationary point, it is a strict local constrained minimum; if it is negative definite, a local constrained maximum; if indefinite, a constrained saddle. If the test is semidefinite, examine higher-order behavior. The key idea is that curvature in a **forbidden** direction cannot classify a constrained optimum.（只能沿「走得通」的方向檢查上彎或下彎。）

In worked example E, feasible changes satisfy $dx=-du$, so the tangent direction is $(-1,1)^T$. Because $H_{zz}=\operatorname{diag}(2,2)$, its curvature along that direction is $(-1,1)H_{zz}(-1,1)^T=4>0$, matching the direct substitution.

The multiplier also measures **sensitivity** to changing a constraint. With this note's convention $H=L+\lambda^Tf$, if the required constraint level is changed from $f=0$ to $f=\delta$, then, to first order at the optimum,

$$dL^*\approx-\lambda^{*T}d\delta.$$

For worked example E, requiring $x+u-2=\delta$ yields $L^*(\delta)=(2+\delta)^2/2$. Its derivative at $\delta=0$ is $2=-\lambda^*$. This is the **shadow price** idea: it tells how much the best achievable cost changes when the requirement shifts slightly. The sign depends on how the constraint is written; if a parameter is added *inside* $f(x,u,c)=0$, then the chain rule gives $dL^*/dc=\lambda^{*T}f_c$.（乘子的正負不能脫離限制式寫法單獨解讀。）

**A source example:** PDF P18–P20 fixes $x=3$ for $L=\tfrac12x^2+xu+u^2+u$. Solving its Hamiltonian equations gives $(x^*,u^*,\lambda^*)=(3,-2,-1)$, $L^*=1/2$, and positive feasible curvature $2$. This illustrates that the constrained optimum need not be the unconstrained optimum.

## 4. Numerical solution methods

**PDF P24, Section 1.3.** For nonlinear objectives and constraints, a closed-form stationary point may be unavailable. **Steepest descent** repeatedly moves opposite a gradient. For an unconstrained problem,

$$
u_{j+1}=u_j-\alpha_jL_u(u_j),
$$

where $j$ is the iteration index and $\alpha_j>0$ is the step size. Since $L_u^T(-\alpha_jL_u)=-\alpha_j\|L_u\|^2\le0$, the direction lowers the cost locally unless the gradient is zero. A large step can overshoot, so a practical method adjusts or checks the step size. Stop when the gradient is small and cost changes have become small, with feasibility checked separately.（負梯度給方向；步長決定一次走多遠。）

For the **constrained** problem in the source, the corresponding direction is $-H_u$, not generally $-L_u$:

1. Choose $u_j$ and solve $f(x_j,u_j)=0$ for a feasible $x_j$.
2. Obtain $\lambda_j=-f_x(x_j,u_j)^{-T}L_x(x_j,u_j)$.
3. Compute $H_u=L_u+f_u^T\lambda_j$.
4. Set $u_{j+1}=u_j-\alpha_jH_u$ and resolve the constraint; stop when the feasible cost change and $\|H_u\|$ are sufficiently small.

This is a conceptual algorithm, not a universal convergence guarantee. It assumes that the constraint can be solved locally and that an appropriate step size is used.

## 5. Connection to optimal control

**Chapter 1 versus later chapters.** Here, $f(x,u)=0$ is a **static** equality constraint. In a time-indexed system, the dynamics themselves become constraints, for example

$$x_{k+1}=f_k(x_k,u_k).$$

For a power converter, $x_k$ may contain inductor current and capacitor/output voltage, and $u_k$ may be a duty-ratio adjustment. A cost can penalize voltage/current tracking error and excessive changes in control effort. The same questions recur: Which choices are feasible under the system equations? What conditions identify a candidate optimum? Does the relevant curvature establish a minimum? The later chapters add time and the coupling between successive states.（第一章是「單一時刻」的數學地基；之後把每個時刻串成動態系統。）

## 6. One-page summary

**Static problem.** Choose $u$ to minimize scalar $L(u)$. The local change is $dL=L_u^Tdu+\tfrac12du^TL_{uu}du+O(\|du\|^3)$. $L_u$ is a column vector; $L_{uu}$ is a matrix. A smooth interior minimum needs $L_u=0$, but this only identifies a stationary point. Positive definite Hessian gives a strict local minimum; negative definite gives a maximum; indefinite gives a saddle; semidefinite requires more work.

**Quadratic problem.** For symmetric $Q$, $L=\tfrac12u^TQu+S^Tu$ has $L_u=Qu+S$ and $L_{uu}=Q$. When $Q$ is invertible, $u^*=-Q^{-1}S$. It is a minimum when $Q\succ0$. For the two-dimensional case, $|Q|$ in the source is $\det Q$.

**Equality constraint.** In $\min L(x,u)$ subject to $f(x,u)=0$, feasible changes obey $f_xdx+f_udu=0$. When $f_x$ is nonsingular, the first-order condition is $L_u-f_u^Tf_x^{-T}L_x=0$. Define $H=L+\lambda^Tf$; the equivalent candidate equations are $H_\lambda=f=0$, $H_x=L_x+f_x^T\lambda=0$, and $H_u=L_u+f_u^T\lambda=0$. Then check curvature along feasible directions. With the convention $f=\delta$, $dL^*/d\delta=-\lambda^*$ locally.

**Numerical idea.** Unconstrained descent uses $u_{j+1}=u_j-\alpha_jL_u$. Constrained descent first satisfies $f=0$, computes $\lambda$, and moves opposite $H_u$ while re-enforcing feasibility. The step size controls progress and stability.

**Four quick checks:** Is the proposed point feasible? Does it satisfy stationarity? Is the Hessian test applied in feasible directions? Are vector and matrix dimensions consistent?

## 7. Glossary

| English term | 繁體中文 | Short meaning |
|---|---|---|
| performance index | 性能指標／目標函數 | scalar cost to optimize |
| decision vector / control vector | 決策向量／控制向量 | adjustable variables $u$ |
| gradient | 梯度 | column vector of first partial derivatives |
| Hessian / curvature matrix | 海森矩陣／曲率矩陣 | matrix of second partial derivatives |
| stationary point | 駐點 | point where first-order change vanishes |
| positive definite | 正定 | positive curvature in every nonzero direction |
| negative definite | 負定 | negative curvature in every nonzero direction |
| indefinite | 不定 | both positive and negative curvature directions |
| equality constraint | 等式限制 | equation $f(x,u)=0$ to satisfy |
| feasible point | 可行點 | point satisfying all constraints |
| Lagrange multiplier | 拉格朗日乘子 | variable balancing constraint and cost gradients |
| Hamiltonian | 漢密頓函數 | here, $H=L+\lambda^Tf$ |
| necessary condition | 必要條件 | must hold at an optimum under its assumptions |
| sufficient condition | 充分條件 | guarantees the stated local conclusion |
| contour / level set | 等高線／等值集合 | points with equal $L$ |
| shadow price / sensitivity | 影子價格／敏感度 | change in optimal cost as a requirement changes |

## 8. Practice problems

These are paraphrased or adapted from the source exercises on **PDF P24–P27**, plus one short check of the worked-example method. Try them before opening Section 9.

1. **Gradient and curvature** (adapted from Problem 1.1-2, PDF P25). For $L(x_1,x_2)=x_1^2-x_1x_2+x_2^2+3x_1$, find the stationary point, the Hessian, its classification, and the minimum value.
2. **When the Hessian is singular** (adapted from Problem 1.1-3, PDF P25). For $g(x,y)=x^2+y^4$, show that the origin is stationary and that its Hessian is singular there. Decide from the original function whether it is a minimum.
3. **Constrained area** (adapted from Problem 1.2-5(a), PDF P26). A rectangle has positive side lengths $x,y$ and perimeter $p>0$. Maximize its area $xy$ subject to $2x+2y-p=0$. Find $x^*,y^*$, the maximum area, and the multiplier using $H=xy+\lambda(2x+2y-p)$.
4. **Feasibility and the Hamiltonian** (based on Section 1.2, PDF P13–P20). Minimize $x^2+u^2$ subject to $x+u-2=0$. Write all three stationary equations, solve them, and explain why $(0,0)$ cannot be accepted.
5. **Two equality constraints** (adapted from Problem 1.2-9, PDF P27). Find the point $z=(x,y,z_3)^T$ nearest the origin that satisfies $3x+2y+z_3=1$ and $x+2y-3z_3=4$. Hint: minimize $\tfrac12z^Tz$ subject to $Az=b$, where $A=\left[\begin{smallmatrix}3&2&1\\1&2&-3\end{smallmatrix}\right]$ and $b=(1,4)^T$.

## 9. Answers to practice problems

1. $L_{(x_1,x_2)}=(2x_1-x_2+3,\,-x_1+2x_2)^T$. Solving both entries equal to zero gives $(x_1^*,x_2^*)=(-2,-1)$. The Hessian is $\left[\begin{smallmatrix}2&-1\\-1&2\end{smallmatrix}\right]$, with eigenvalues $1$ and $3$, hence positive definite. Therefore it is the unique global minimum of this quadratic, with $L^*=-3$.
2. $\nabla g=(2x,4y^3)^T$, so the origin is stationary. $\nabla^2g(0,0)=\operatorname{diag}(2,0)$ is singular and the second-order test is inconclusive. Directly, $x^2+y^4\ge0$ with equality only at $(0,0)$, so the origin is a strict global minimum.
3. $H_x=y+2\lambda=0$, $H_y=x+2\lambda=0$, and $H_\lambda=2x+2y-p=0$. Thus $x=y=p/4$ and $\lambda=-p/8$. The maximum area is $p^2/16$. Along the feasible line $y=p/2-x$, area is $x(p/2-x)$, a concave quadratic, confirming the maximum.
4. $H=x^2+u^2+\lambda(x+u-2)$, so $H_x=2x+\lambda=0$, $H_u=2u+\lambda=0$, and $H_\lambda=x+u-2=0$. Hence $(x^*,u^*,\lambda^*)=(1,1,-2)$ and $L^*=2$. $(0,0)$ violates $x+u-2=0$, so its smaller unconstrained cost is irrelevant.
5. With $H=\tfrac12z^Tz+\lambda^T(Az-b)$, stationarity gives $z=-A^T\lambda$ and feasibility gives $-AA^T\lambda=b$. Since $AA^T=\left[\begin{smallmatrix}14&4\\4&14\end{smallmatrix}\right]$, $\lambda=\left(1/90,-13/45\right)^T$ and $z^*=\left(23/90,5/9,-79/90\right)^T$. Direct substitution verifies both equations. Its distance from the origin is $\|z^*\|=\sqrt{103/90}$. The squared norm is strictly convex and the constraints are linear, so this is the unique global minimum.
