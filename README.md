# Computational Repository: Modeling and Monitoring of Service Resilience and Reliability for Wireless Cellular Networks

This repository contains all implementation files and datasets supporting the research published in IISE Transactions under the title:  
**Modeling and Monitoring of Service Resilience and Reliability for Wireless Cellular Networks through a Probabilistic Geometrical Coverage Model**

## Repository Contents

### 📁 Data Files
- `'Dataset1.mat'`: Aggregated user demand data collected by base stations  
- `'Dataset2.mat'`: User demand data, base station geolocations, and reliability metrics  
- `'Dataset3.mat'`: User demand data, base station geolocations, reliability metrics, and alarm data  

### 🧠 Main Scripts
1. `'CodeFile1.m'`:  
   - Implements Equations 13-17  
   - Fits periodic functions to user demand patterns using `'Dataset1.mat'`  
   - Generates **Figure 4**  

2. `'CodeFile2.m'`:  
   - Implements Equations 18-19  
   - Models spatiotemporal demand distribution using `'Dataset2.mat'`  
   - Generates **Figure 5**  

3. `'CodeFile3.m'`:  
   - Implements Equations 21-24  
   - Computes holistic service reliability/resilience using `'Dataset2.mat'`  
   - Generates **Table 3** results  

4. `'CodeFile4.m'`:  
   - Implements Equations 21-24  
   - Computes reliability/resilience under dynamic alarms using `'Dataset3.mat'`  
   - Generates **Table 4** results  

### 🔧 Support Functions
#### Reliability Calculation (Equations 1-3)
- `'reliability.m'`, `'reliability_.m'`
#### Periodic Demand Fitting (Equations 13-17)
- `'f_num.m'`, `'f_num_.m'`, `'fit_f.m'`, `'fit_f_.m'`, `'fit_function.m'`, 
`'fit_function_.m'`, `'fit_plot.m'`, `'fit_plot_.m'`, `'regpoly0.m'`,
`'regpoly1.m'`, `'regpoly2.m'`,`'predictor.m'`
#### Spatiotemporal Correlation (Equations 18-19)
- `'corrcubic.m'`， `'correxp.m'`, `'correxpg.m'`, `'corrgauss.m'`, 
`'corrlin.m'`, `'corrspherical.m'`, `'corrspline.m'`, `'variogram.m'`
#### Visualization Tools
- `'plotcircle.m'`, `'plotkriging.m'`
#### Holistic Metrics (Equations 21-24)
- `'holistic_measure.m'`, `'holistic_plot.m'`
## 🖥 Experimental Environment (As Used in Section 7)
- Operating System: Windows 11 Pro

- MATLAB Version: 2022b

- Hardware Configuration:

  - CPU: Intel® Core™ i9-12900H (14 cores, 20 threads @ 2.50GHz)

  - GPU: NVIDIA GeForce RTX 3060 (6GB GDDR6)

  - RAM: 16GB DDR5

  - Storage: 1TB NVMe SSD
## Contact
For technical inquiries or dataset access issues, please open a GitHub Issue:
https://github.com/ZongqiXue/PGCM
