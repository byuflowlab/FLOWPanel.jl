---
title: Lamb Vector in the Euler Pressure Equation
author: Ryan Anderson
date: 19 May 2026
geometry: margin=1in
output: pdf_document
---
# Derivation

Begin with the Euler momentum equation for inviscid, adiabatic flows:
$$\begin{align}
\frac{D \vec{u}}{Dt} &= -\nabla \left( \frac{p}{\rho} \right) + \vec{g}
\end{align}$$
where $\vec{u}$ is the velocity, $p$ is pressure, $\rho$ is density, and $\vec{g}$ acceleration due to additional field forces (like gravity). Assuming incompressible flow, expanding the total derivative, and ignoring field forces, we have

$$\begin{align}
\frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla) \vec{u} &= - \frac{1}{\rho} \nabla p
\end{align}$$
Then, taking the divergence of both sides provides us with a Poisson equation which may be solved for a pressure field $p$:

$$\begin{align}
\nabla \cdot \left( \frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla) \vec{u} \right) &= - \frac{1}{\rho} \nabla^2 p
\end{align}$$

## Enter the Lamb Vector

The last equation can be decomposed in terms of the Lamb vector $\vec{l} \equiv \vec{u} \times \vec{\omega}$ by deriving the identity 

$$\begin{align}
(\vec{u} \cdot \nabla) \vec{u} &= (\vec{u} \cdot \nabla) \vec{u} + \vec{u} \times (\nabla \times \vec{u}) - \vec{u} \times (\nabla \times \vec{u})\\
&= \frac12 \left( \vec{u} \cdot \vec{u} \right) - \vec{u} \times (\nabla \times \vec{u} )\\
&= \frac12 \left( \vec{u} \cdot \vec{u} \right) - \vec{u} \times \vec{\omega}
\end{align}$$
which means the Poisson equation above can be decomposed into 

$$\begin{align}
\nabla \cdot \left( \frac{\partial \vec{u}}{\partial t} + \frac{|\vec{u}|^2}{2} - \vec{u} \times \vec{\omega} \right) &= - \frac{1}{\rho} \nabla^2 p\\
\nabla \cdot \left( \rho \frac{\partial \vec{u}}{\partial t} + \frac12 \rho |\vec{u}|^2 - \rho \vec{u} \times \vec{\omega} \right) &= - \nabla^2 p
\end{align}$$

which for steady, irrotational flows, reduces to

$$\begin{align}
\nabla \cdot \left( \frac12 \rho |\vec{u}|^2 \right) &= - \nabla^2 p\\
\frac{1}{2}\rho |\vec{u}|^2 &= -\nabla p
\end{align}$$
which is dimensionally inconsistent with Bernoulli's equation

$$\begin{align}
\frac12 \rho |\vec{u}|^2 + \rho g h + p &= \mathrm{constant}
\end{align}$$

Why?