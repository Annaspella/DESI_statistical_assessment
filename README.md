# 📊 Unveiling the Digital Frontier: An In-Depth Analysis of the Digital Economy and Society Index (DESI)

**Author**: Anna Gotti  
**Degree**: MSc in Statistical Sciences  
**Institution**: University of Padova  
**Supervisors**: Dr. Manuela Scioni (Unipd), Dr. Paola Annoni (European Commission - DG CONNECT)  
📄 [Full thesis PDF](./Gotti_Anna.pdf)

---

## 🧠 Project Overview

This repository contains the full statistical and computational analysis supporting my master’s thesis:  
**“Unveiling the Digital Frontier: An In-Depth Analysis of the Digital Economy and Society Index”**.

The project stems from my research internship at the European Commission (DG CONNECT) and focuses on the **methodological validation and enhancement of the DESI 2022 framework**—a composite index used to measure digital performance across EU member states.

---

## 🌐 About DESI 2022

The **Digital Economy and Society Index (DESI)** monitors Europe’s overall digital performance and tracks the progress of EU countries in their digital competitiveness. The 2022 edition mainly reflects data collected in 2021 and includes recalculations for previous years due to structural changes in indicators, corrections in underlying data, and the removal of the UK from the EU average.

### DESI 2022 Dimensions

DESI is structured around **four key dimensions**, each representing a core area of digital progress:

| Acronym | Full Name                       | Description                                                                                 |
|---------|--------------------------------|---------------------------------------------------------------------------------------------|
| HC      | Human Capital                  | Measures digital skills and internet usage among individuals, including ICT specialists.   |
| CN      | Connectivity                  | Assesses the quality and coverage of digital infrastructure, such as broadband and 5G.      |
| IDT     | Integration of Digital Technology | Evaluates how businesses adopt digital technologies to improve operations and sales.       |
| DPS     | Digital Public Services       | Captures the availability and usage of online public services, including eGovernment tools.|

---

## 🎯 Research Goals

- Assess the **internal consistency** of DESI 2022 dimensions and sub-dimensions.
- Explore **latent structures** using PCA and **ROBPCA** to address outliers.
- Provide **robust score recalculations** for dimensions based on proposed adjustments.
- Conduct **sensitivity analysis** on indicator weights and new sub-dimensions (especially for Digital Public Services).
- Evaluate **ranking and score impacts** under various methodological scenarios.

---

## 🧪 Methods

- **Cronbach’s Alpha** – for internal reliability.
- **Principal Component Analysis (PCA)** – for latent dimension structure.
- **Robust PCA (ROBPCA)** – to detect multivariate outliers.
- **Sensitivity & Impact Analysis** – for weight testing and robustness checks.
- **Min-max scaling** – for score normalization using fixed min/max benchmarks.

---

## 🗂️ Repository Structure

- 📁 data/                        # Source data for DESI 2022 and revised indicators
- 📁 scripts/                     # MATLAB .m files for PCA, ROBPCA, scoring, and ranking
- 📁 results/                     # Outputs, scores, figures, and tables
- 📄 Gotti_Anna.pdf               # Full thesis document (link above)
- 📄 README.md                    # This file

---

## 🔍 Code & Data Description

### 📂 MATLAB Scripts

- **DESI_analysis_$.m**  
  Internal consistency analysis of `$` dimension (HC, CN, IDT, DPS).  
  DPS has two versions using 2021 and 2022 data separately.

- **DESI_Scores_$.m**  
  Score calculation scripts for `$` dimension and simulation of revised scoring scenarios.

- **get_rank.m**  
  Function to compute rankings from computed scores.

- **get_score.m**  
  Function to handle scoring in presence of missing (`NaN`) values.

---

### 📊 Excel Files

- **calculations_minmax_egov_2022data.xlsx / _2023data.xlsx**  
  Min-max thresholds used for normalization of new eGov indicators.

- **min_max.xlsx**  
  Aggregated table with min, max, and weights for all DESI and new indicators.

- **Weights_Subdim.xlsx**  
  Comparison of original and proposed sub-dimensional weights.

---

### 📈 CSV Files

- **DESI Y6.csv**  
  Raw dataset for DESI 2022 (no missing data).

- **DESI Y7_2.csv**  
  Updated dataset using new 2022 data and revised indicator 4a1 (`I_IGOVANYS`).

- **new_indicators_$.csv**  
  New indicator values for each `$` dimension (e.g. DPS, CN, etc.).

---

## 📌 Key Findings

- DESI dimensions vary significantly in **internal coherence**, raising concerns about conceptual validity.
- ROBPCA analysis detected **multivariate outliers** in several indicators, justifying need for robust statistics.
- Revised sub-dimensional weights and indicator inclusion strategies **altered country rankings**, especially in the DPS dimension.
- Sensitivity analysis highlighted **non-trivial trade-offs** between indicator structure and policy signaling.

---

## 📚 References

- European Commission’s [DESI 2022 Report](https://digital-strategy.ec.europa.eu/en/policies/desi)
- Hubert, M., Rousseeuw, P. J., & Vanden Branden, K. (2005). **ROBPCA: A new approach to robust principal component analysis**.
- OECD Handbook on Constructing Composite Indicators

---

## 📬 Contact

📧 anna.gotti16@gmail.com

---

## 📎 Citation

If you use this work or build upon it, please cite:

**Gotti, A. (2023).** *Unveiling the Digital Frontier: An In-Depth Analysis of the Digital Economy and Society Index (DESI)*. MSc Thesis, University of Padova.  