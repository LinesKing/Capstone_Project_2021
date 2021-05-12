# Capstone_Project_2021



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#Optimal-Path-Following">Optimal Path Following</a>
      <ul>
        <li><a href="#Linearised-NMPC">Linearised NMPC</a></li>
        <li><a href="#Nonlinear-MPC">Nonlinear MPC</a></li>
      </ul>
    </li>
  </ol>
</details>




<!-- OPTIMAL PATH FOLLOWING -->
## Optimal Path Following
We both use N = 20, 
$Q = \begin{bmatrix}
20 & 0 & 0 & 0\\
0 & 20 & 0 & 0\\
0 & 0 & 0 & 0\\
0 & 0 & 0 & 0
\end{bmatrix}$, 
$Q = \begin{bmatrix}
5 & 0\\
0 & 5
\end{bmatrix}$, 

$\begin{bmatrix} u\\e \end{bmatrix}  \sim 
N (\begin{bmatrix}0\\0\end{bmatrix},
\begin{bmatrix} G(\sigma_g) &0\\ 0 &R(\sigma_{\gamma})\end{bmatrix})$

### Linearised NMPC

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Linearised NMPC/history of position and track.jpg" alt="History of position and track in linearised NMPC" width="560" height="420">
  </a>
</p>

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Linearised NMPC/histroy of state.jpg" alt="History of state in linearised NMPC" width="560" height="420">
  </a>
</p>

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Linearised NMPC/history of input.jpg" alt="History of input in linearised NMPC" width="560" height="420">
  </a>
</p>

### Nonlinear MPC

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Nonlinear MPC/histroy of position and track.jpg" alt="History of position and track in nonlinear MPC" width="560" height="420">
  </a>
</p>

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Nonlinear MPC/history of state.jpg" alt="History of state in nonlinear MPC" width="560" height="420">
  </a>
</p>

<br />
<p align="center">
  <a href="https://github.com/LinesKing/Capstone_Project_2021">
    <img src="Optimal Path Following/Nonlinear MPC/history of input.jpg" alt="History of input in nonlinear MPC" width="560" height="420">
  </a>
</p>
