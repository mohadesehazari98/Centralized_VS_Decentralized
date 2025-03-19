🚀 Average Entanglement Distribution Rate  
Centralized vs. Decentralized Architectures in Quantum Networks  

📡 This repository provides MATLAB code for analyzing entanglement distribution rates in quantum networks, comparing centralized and decentralized strategies for GHZ state distribution in 2D repeater architectures.  

---

🏗 Folder Structure  
📂 Average_Entanglement_Distribution_Rate/  
📄 Rate_Cent.m → Computes centralized GHZ entanglement rate.  
📄 Rate_Decent.m → Computes decentralized GHZ entanglement rate.  
📄 Rate_2D.m → Calculates rate for multi-hop/multi-level 2D architectures.  
📄 link_gen_prob.m → Computes link success probability (q_link) based on network geometry.  
📄 F_T_max.m → Recursive function for hierarchical GHZ fusion (CDF calculations).  
📄 README.txt → You’re reading this!  

---

🎯 Project Overview  

We analyze entanglement distribution in quantum repeater networks using two strategies:  

1️⃣ Centralized GHZ Distribution  
🔹 A central node (switch) collects photons and distributes GHZ states.  
🔹 Requires long-distance teleportation.  
🔹 Higher latency but potentially more efficient with high-fidelity links.  

2️⃣ Decentralized GHZ Distribution  
🔹 GHZ states are fused at intermediate repeater nodes.  
🔹 Uses local fusion operations rather than a central switch.  
🔹 Potentially faster and scalable, reducing bottlenecks.  

---

🏗 Mathematical Model  

The entanglement distribution rate is modeled as:  

    Rate = (GHZ teleportation success probability) / (Expected Total Time)

Rate functions incorporate:  
- Tail-sum formulas  
- Nested summations  
- Recursive CDF generation  

💡 Key Parameters:  
| Parameter | Description |
|-----------|-------------|
| N | Total nodes/qubits in the network |
| m | Number of hierarchy levels |
| q_link | Probability of successful entanglement link |
| q_BSM | Bell State Measurement success probability |
| q_Fuse | Fusion success probability |
| L₀_in | Distance between neighboring entangled nodes |

---

▶️ Running the Code  

Step 1: Clone the Repository  
    git clone https://github.com/mohadesehazari98/Centralized_VS_Decentralized.git  
    cd Centralized_VS_Decentralized/Average_Entanglement_Distribution_Rate  

Step 2: Run in MATLAB  
    Open MATLAB and execute:  
        example_plot_rates  

This will generate comparison plots showing the rate vs. distance for different architectures.  

---

📊 Sample Output  

The code produces visualizations comparing Centralized vs. Decentralized GHZ distribution:  
✅ Rate vs. Distance for 1D and 2D networks  
✅ Effect of hierarchy levels (m = 1, 2, 3)  
✅ Trade-offs between network designs  

---

🔧 Dependencies  
- MATLAB R2020+ (No additional toolboxes required)  

---

👩‍💻 Author  
Mohadeseh Azari  
📧 mohadesehazari1998@gmail.com  

---

📜 License  
Released under the MIT License – Feel free to use, modify, and cite!  

⭐️ If you find this useful, please give it a star!  

🚀 Quantum networks are the future—let’s optimize them together! 🚀  
