# mixed-onebit-doa-crb

This repository contains the source code and experimental implementations for the paper "Covariance Recovery and CRB Analysis for Mixed One-Bit Quantization".

## Experimental Instructions

### MUSIC Simulation Experiments

- **Entry point**: `main.m`
- **Usage**: Run this file using MATLAB and select the corresponding experiment type when prompted.

### Cramér-Rao Bound (CRB) Experiments

The following files contain CRB analysis experiments for different array configurations:

#### CRB vs SNR Experiments:
- `CRB_vs_SNR_uls.m` - CRB analysis vs SNR for Uniform Linear Array (ULA)
- `CRB_vs_SNR_nla.m` - CRB analysis vs SNR for Nested Linear Array (NLA)  
- `CRB_vs_SNR_cop.m` - CRB analysis vs SNR for Coprime Array (COP)

#### CRB vs Number of Snapshots Experiments:
- `CRB_vs_N_ula.m` - CRB analysis vs number of snapshots for Uniform Linear Array (ULA)
- `CRB_vs_N_nla.m` - CRB analysis vs number of snapshots for Nested Linear Array (NLA)
- `CRB_vs_N_cop.m` - CRB analysis vs number of snapshots for Coprime Array (COP)

## Getting Started

1. Ensure MATLAB is installed on your system
2. Clone this repository
3. For MUSIC simulations: Run `main.m` and follow the prompts
4. For CRB analysis: Execute the desired CRB experiment file directly

## Array Configurations

- **ULA**: Uniform Linear Array
- **NLA**: Nested Linear Array  
- **COP**: Coprime Array