## Discrete elastic Rod (DER)

Hangqi Cui - Nada Abdelwahab - Suneh Bhatia

---

### What is a Discrete Elastic Rod?

- Kirchhoff rod

  ![img](https://lh7-us.googleusercontent.com/2MaIrwJrbO8lpqlyAvPzHfpa1zQTUCvOBxygg2FZSjGa5_Otd5P9TYInMrpJ7wYkIMDfoTEJoONN65qHNycxDxk3dKA4BNWfQJk5NMw6NhB-qmFHnhO_IdB7mZvQ45gl7gAA_HEfquuDZKH2VyKX1-wdZw=s2048)
- Discrete polyline (primals and duals )

  ![img](https://lh7-us.googleusercontent.com/ba5RUxk5udFDIljpFOlmJCPK6Cg1X_Oz5y-9UsNPsXyWWMGDu9fG2SwcHAMvl_81Jo01CbQotjIeid2Q2LYCTFTTv74vfv1ULTIBpMFFF5b-y--LeopAXZtaiUtwakOoD4ZGnPt1lI9P2rxJB_aKwT3i5Q=s2048)

---

### Solve Dynamics!

- Cable/Rope
- Viscous thread
- Hair simulation

> Energy $\rightarrow$ Force $\rightarrow$ Accelaration $\rightarrow$ update position

---

### Bishop Frame

> Twist is independent from the curve

Tangnet is always certain in the discrete curve. 
We set up a reference frame (normal) at the beginning of the curve.
$$
N = (0,0,1),\quad B = T \times N
$$

---

### Parallel Transport
Construct the remaining normals iteratively. 

$$
N^{i+1} = P^i * N^i
$$

---

### Material Frame 
> Simply rotate u and v in twisting angle $\theta$

---

### Curvature and Darboux vector 

> Surprisingly Many ways to approximate curvature, and receive difference reuslt in the end.

$$
\kappa = 2 \tan \frac{\theta}{2}
$$

Darboux vector is curvature binormal. 
$$
\Omega = \kappa b
$$

---

### Dynamics from high school 
$$
F = dE \qquad E = \int F dx 
$$
$$
F = ma 
$$

---

### Energies and Forcies (continous case)
- Elastic energy : 0 
- Bending energy :

- Twisting energy :
  $$ {1\over2} \int \beta (\theta^ \prime)^2 ds $$

- $\alpha, \beta$: Modolus, Material properties

---

### Energies and Forcies (discrete case)
- Elastic energy : still 0 
- Bending energy : 
  $$ E_{\text {bend }}(\Gamma)=\frac{1}{2} \sum_{i=1}^n \alpha\left(\frac{\kappa \mathbf{b}_i}{\bar{l}_i / 2}\right)^2 \frac{\bar{l}_i}{2}=\sum_{i=1}^n \frac{\alpha\left(\kappa \mathbf{b}_i\right)^2}{\bar{l}_i} $$
- Twisting energy : 
  $$ E_{\text{twist}}(\Gamma)=\sum_{i=1}^n \beta \frac{\left(\theta^i-\theta^{i-1}\right)^2}{\bar{l}_i}=\sum_{i=1}^n \frac{\beta m_i^2}{\bar{l}_i} $$

- $l_i = |e^{i-1}| +|e^i|$ , 'weight' of the vertex. 

---

### Symplectic Euler

$$
a = {dv\over dt}, \qquad v = \int a dt
$$
$$
v = {dx\over dt}, \qquad x = \int v dt
$$
$$
F = ma \quad 
$$

---

### Constraints
- Manifold Projection 
  - Remeber our curve is inextenisble? 
- Twist holonomy 
  - update the twist angle on each vertex

---

### Experiments 
