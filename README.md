# 📊 Unveiling the Digital Frontier: An In-Depth Analysis of the Digital Economy and Society Index (DESI)


**Author:** Anna Gotti  
**Affiliation:** VU Amsterdam | DG CONNECT, European Commission  
**Thesis Project | MSc in Statistical Sciences**  
**Academic Year:** 2024-2025

---

## 🧠 Overview

This repository contains the code and documentation for my Master's thesis project developed during a research internship at the **Directorate-General for Communications Networks, Content and Technology (DG CONNECT)** of the **European Commission**. The work critically analyses the **Digital Economy and Society Index (DESI) 2022**, a composite indicator that evaluates how digital transformation affects EU economies and societies.

My goal was to explore the **statistical soundness and internal consistency** of DESI’s methodological framework, focusing on latent variable modeling and robustness in the presence of outliers.

---

## 🎯 Objectives

- **Evaluate** the internal consistency of DESI and its dimensions.
- **Justify** or challenge the aggregation rules used in computing DESI scores.
- **Compare** classical PCA with **Robust PCA (ROBPCA)** for latent factor analysis.
- **Assess** how dimensional coherence and weighting schemes influence rankings.
- **Suggest** improvements for the structure and methodology of DESI to better support digital policy decision-making.

---

## 🔧 Methods

- **Cronbach’s Alpha** – to assess internal reliability of dimensions.
- **Principal Component Analysis (PCA)** – to explore latent structure.
- **Robust PCA (ROBPCA)** – to detect and handle outliers in high-dimensional space.
- **Pearson Correlation Coefficient** – for assessing indicator inter-relationships.
- **Sensitivity Analysis** – to evaluate the impact of structural and weight adjustments on country scores.

---

## 🗂️ Repository Structure

📁 data/ # Raw and cleaned DESI datasets
📁 scripts/ # Analysis scripts (R or Python)
📁 results/ # Output tables, plots, PCA loadings
📁 figures/ # Graphs and charts for the thesis
📄 thesis_summary.pdf # Summary of the thesis (Italian)
📄 README.md # Project overview and instructions

---

## 📈 Key Insights

- Some DESI dimensions exhibit **low internal coherence**, which may challenge their conceptual validity.
- PCA and ROBPCA reveal **different latent structures**, particularly in the presence of outliers.
- Adjustments to **subdimension structure and weights** significantly affect the **country rankings**, suggesting a need for a more data-driven and policy-aligned design.

---

## 📚 Background

Composite indicators like DESI are increasingly used in policy-making, yet their reliability depends on both **conceptual soundness** and **methodological rigor**. By grounding the analysis in **statistical validation techniques**, this project contributes to a more **transparent**, **robust**, and **policy-relevant** index design.

---

## 🔗 Related Work

- European Commission’s [DESI 2022 Report](https://digital-strategy.ec.europa.eu/en/policies/desi)
- Hubert, M., Rousseeuw, P. J., & Vanden Branden, K. (2005). **ROBPCA: A new approach to robust principal component analysis**.
- OECD Handbook on Constructing Composite Indicators

---

## 💡 Future Work

- Extend analysis to DESI 2023 and compare longitudinal effects.
- Incorporate machine learning for dimensionality reduction and classification.
- Develop a dashboard or R Shiny app for interactive policy simulations.

---

## 🤝 Acknowledgments

Thanks to my supervisor and the team at **DG CONNECT, European Commission**, for their guidance and data access. Special thanks to the VU Amsterdam Department of Economics for supporting this project.

---

## 📬 Contact

📧 ann.gotti16@gmail.com  

