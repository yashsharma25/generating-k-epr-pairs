1. To numerically examine all possible choices of all pairs for 10-qubit graph state and show that at least one choice of Pauli measurement generates target EPR state stabilizers on 2 pairs, using a python script. (10_qubit_verify.py)
2. Verify that no 9-qubit graph is 2-pairable. Using graph state database available at https://zenodo.org/record/3757948#.YxoKGS-B1QJ  (code in 9_qubit_verify.py)
3. Verifies that graph state circuit is indeed 2-pairable using Qiskit simulation  (code in verify.ipynb). This is done in 2 ways:\
    i) Save the statevectors in the simulator and verify if they are 2-pairable\
    ii) Reconstruct states using Quantum State Tomography for each measurement result (2^6 = 64 possible measurement), and verify if 2-pairable