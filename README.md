# resource_intxns

Code for the statistical models and analyses in

Arrabal A, San Román M, Diaz-Colunga J & Sanchez A (2026). Simple rules underlie complex interactions among bacterial resources. *biorXiv.*

Please refer questions to [juan.diaz\@ipla.csic.es](mailto:juan.diaz@ipla.csic.es) or [juan.diazcolunga\@gmail.com](mailto:juan.diazcolunga@gmail.com).

## Pipeline

The `./data` directory contains experimental data (measurements of bacterial fitness of 7 isolates in all 255 mixed-nutrient environments).

The `./scripts_lin_models` directory contains code to perform the analyses described in the main text:

-   `main.R` executes the main pipeline. In short, it loads the files in `./data` and performs the analyses described in the main text: splits each dataset into a training and a test set (with training fractions of varying sizes), uses the training set to train linear first-order linear models of the form

    $$
    \hat{f}(x) = f_0 + \sum_i f_i x_i
    $$

    or second-order models of the form

    $$
    \hat{f}(x) = f_0 + \sum_i f_i x_i + \sum_{j>i} f_{ij} x_i x_j
    $$

    using either Lasso or Ridge regularization to avoid overfitting. Model performance is then evaluated by comparing model predictions with empirical observations within the test set. These outputs are stored as .RData objects in the `./results` directory (files labeled as `eval_linmod_***.RData`).

-   `make_plots.R` loads these .RData outputs and produces the plots shown in the manuscript figures. Plots are saved in the `./plots/linmods` directory.

The `./scripts_lin_vs_fees` directory follows a similar structure. In this case, `main.R` does not only train linear models but also the global epistasis-based statistical method described in [1]. Again, outputs are saved in the `./results` directory (files labeled as `eval_data_***.RData`), and the `make_plots.R` script loads these files to produce the plots shown in the Supplementary Material of the manuscript. The `ecoFunctions.R` file is taken from [1] and contains auxiliary code to implement the the global epistasis-based statistical method.

[1] Diaz-Colunga J, Skwara A, Vila JCC, Bajic D & Sanchez A (2024). Global epistasis and the emergence of functionin microbial consortia. *Cell* **187**, 3108-3119.e30. [doi.org/10.1016/j.cell.2024.04.016](https://doi.org/10.1016/j.cell.2024.04.016)
