
📡 Average Entanglement Distribution Rate — Centralized vs Decentralized Architectures

This repository contains MATLAB code and analysis for evaluating the average entanglement distribution rate in quantum networks, focusing on Centralized vs Decentralized GHZ state distribution strategies across 1D and 2D repeater architectures.

📂 Folder Structure

Average_Entanglement_Distribution_Rate/
│
├── Rate_Cent.m             % 1D Centralized GHZ distribution rate calculation
├── Rate_Decent.m           % 1D Decentralized GHZ distribution rate calculation
├── Rate_2D_Cent.m          % 2D Centralized GHZ distribution rate (multi-level/multi-hop)
├── Rate_2D_Decent.m        % 2D Decentralized GHZ distribution rate (multi-level/multi-hop)
├── link_gen_prob.m         % Helper function to compute q_link based on geometry and type
├── F_T_max.m               % Recursive function to compute CDF F_T_max(k) for multilevel hierarchy
├── example_plot_rates.m    % Example script to generate plots comparing rate vs distance
└── README.txt              % Documentation file

🧠 Project Objective

To evaluate and compare the scaling behavior and performance of GHZ entanglement distribution under:

- Centralized GHZ generation (all photons routed to a central fusion node)
- Decentralized GHZ generation (photons fused hierarchically at distributed fusion nodes)

The study includes both 1D chains and 2D hierarchical repeater networks, with emphasis on the impact of parameters like:

- Number of nodes/qubits (N)
- Number of hierarchy levels (m)
- Link success probability q_link
- Bell State Measurement success q_BSM
- Fusion success probability q_Fuse
- Neighboring distance L_0

📈 Rate Model Description

Each rate function models the expected maximum completion time for successful GHZ distribution and computes the average rate as:

Rate = (GHZ success probability) / (Expected Total Time)

Rate functions include:
- Tail-sum formulas
- Nested summations
- Recursive CDF generation for hierarchical levels

▶️ How to Run

1. Clone the repository:
   git clone https://github.com/mohadesehazari98/Centralized_VS_Decentralized.git
   cd Centralized_VS_Decentralized/Average_Entanglement_Distribution_Rate

2. Open MATLAB and run:
   example_plot_rates

This will generate comparison plots of rate vs distance for different architecture types and hierarchy levels.

📊 Example Output

- Rate vs Distance curves for:
  - Centralized vs Decentralized (1D and 2D)
  - Different hierarchy levels (m = 1, 2, 3)

These help visualize the trade-offs between different network designs.

📎 Dependencies

- MATLAB R2020+ or later
- No external toolboxes required

✍️ Author

Mohadese Hazari  
📧 mohadesehazari98@gmail.com

📜 License

This project is released under the MIT License. Feel free to use, modify, and cite.

⭐️ If you find this work useful...
Please consider giving a star to the repository and sharing it with your peers!
