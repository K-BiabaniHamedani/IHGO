# Improved Hybrid Growth Optimizer (IHGO)



## Abstract
This paper presents an efficient method for frequency-constrained optimization of large-scale cyclically symmetric domes. The approach integrates the improved hybrid growth optimizer (IHGO) algorithm with an eigenvalue decomposition method. IHGO incorporates the exploration mechanism of the improved arithmetic optimization algorithm (IAOA) into its learning phase, along with algorithm-specific modifications. While these modifications are general and problem-independent, their effectiveness in broader structural optimization tasks remains unexplored. To enhance computational efficiency, a decomposition-based method performs free vibration analysis. This method partitions the eigenvalue problem into smaller, decoupled sub-eigenproblems through block-diagonalization of structural matrices, significantly reducing CPU time and memory requirements compared to standard methods. The performance of IHGO is demonstrated via optimization of two large-scale domes, comparing results against the original growth optimizer (GO) and literature-best solutions. These comparisons highlight the outstanding computational cost and accuracy of IHGO. The results confirm the robustness and computational advantages of IHGO, establishing it as a powerful tool for large-scale structural optimization under natural frequency constraints.



## Reference
Kaveh, A., and Biabani Hamedani, K. (2024). A hybridization of growth optimizer and improved arithmetic optimization algorithm and its application to discrete structural optimization. Computers & Structures, 303, 107496. https://doi.org/10.1016/j.compstruc.2024.107496



## Acknowledgements
This work is based upon research funded by Iran National Science Foundation (INSF) under project No. 4024911. 
