# Frequency-constrained optimization of large-scale cyclically symmetric domes using improved hybrid growth optimizer



## Abstract
This paper presents an efficient method for frequency-constrained optimization of large-scale cyclically symmetric domes. The approach integrates the improved hybrid growth optimizer (IHGO) algorithm with an eigenvalue decomposition method. IHGO incorporates the exploration mechanism of the improved arithmetic optimization algorithm (IAOA) into its learning phase, along with algorithm-specific modifications. While these modifications are general and problem-independent, their effectiveness in broader structural optimization tasks remains unexplored. To enhance computational efficiency, a decomposition-based method performs free vibration analysis. This method partitions the eigenvalue problem into smaller, decoupled sub-eigenproblems through block-diagonalization of structural matrices, significantly reducing CPU time and memory requirements compared to the standard method (which solves the full eigenvalue problem without decomposition). The performance of IHGO is demonstrated via optimization of two large-scale domes, comparing results against the original growth optimizer (GO) and literature-best solutions. These comparisons highlight the outstanding computational efficiency and accuracy of IHGO. The results confirm the robustness and computational advantages of IHGO, establishing it as a powerful tool for large-scale structural optimization under natural frequency constraints.



## Pseudocode of the Improved Hybrid Growth Optimizer (IHGO) Algorithm
![Figure 3](https://github.com/user-attachments/assets/6847fc70-4259-4691-8f88-5287d761701f)



## Flowchart of the Improved Hybrid Growth Optimizer (IHGO) Algorithm
![Figure 4](https://github.com/user-attachments/assets/9c0fc9be-902b-459a-a552-a621c4613b9a)



## Reference
Kaveh, A., Biabani Hamedani, K., and Hosseini, S.M. (2025). Frequency-constrained optimization of large-scale cyclically symmetric domes using improved hybrid growth optimizer. Scientia Iranica. https://doi.org/10.24200/sci.2025.66274.9952



## Acknowledgements
This work is based upon research funded by Iran National Science Foundation (INSF) under project No. 4024911. 
