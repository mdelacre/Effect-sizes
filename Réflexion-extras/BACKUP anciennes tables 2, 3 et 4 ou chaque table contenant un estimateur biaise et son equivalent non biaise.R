---
title: "Backup"
output: html_document
---

# Table 2

\newpage
\blandscape

| Estimator|      expectancy    |          Bias             |               Variance                     |  
  |--------|:-----------------:  | :-----------------------:  | :------------------------------------------------------------:  |
  | \tiny$Glass's \; d_s$    | \tiny$\delta_{glass} \times \frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}$  |   \tiny$\delta_{glass} \times \left( \frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}-1 \right)$       |\tiny$\frac{df}{df-2} \times \left( \frac{1}{n_c} + \frac{\sigma^2_e}{n_e\sigma^2_c}\right) + \delta^2_{Glass} \left[ \frac{df}{df-2} - \left( \frac{\sqrt{\frac{df}{2}} \times \Gamma\left( \frac{df-1}{2}\right)}{\Gamma\left( \frac{df}{2}\right)} \right)^2 \right]$ |     |
|    |      |  |    |
| \tiny$Glass's \; g_s$    | \tiny$\delta_{glass}$  |   /     |\tiny$Var(Glass's \; d_s) \times \left[\frac{\Gamma(\frac{df}{2})}{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}\right]^2$ |     |
|    |      |  |    |
Table: Expentency, bias and variance of Glass's $d_s$ and Glass's $g_s$ under the assumptions that independent residuals are normally distributed.

*Note*. $df =n_c-1$. Equations in Table 2 require $df \ge 3$ (i.e. $n_c \ge 4$), and $n_e \ge 2$.

\elandscape

\newpage

# Table 3

\newpage
\blandscape

| Estimator|      expectancy    |          Bias         |                 Variance                       |  
|---------|:-----------------:  | :---------------------:  | :-------------------------------------------:  |
| \tiny$Cohen's \; d'_s$    | \tiny$\delta'_{Cohen} \times \frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}$  |   \tiny$\delta'_{Cohen} \times \left( \frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}-1 \right)$       |\tiny$\frac{df}{df-2} \times \frac{2\left( \frac{\sigma^2_1}{n_1} + \frac{\sigma^2_2}{n_2} \right)}{\sigma^2_1+\sigma^2_2} + (\delta'_{Cohen})^2 \left[ \frac{df}{df-2} - \left( \frac{\sqrt{\frac{df}{2}} \times \Gamma\left( \frac{df-1}{2}\right)}{\Gamma\left( \frac{df}{2}\right)} \right)^2 \right]$ |     |
  |    |      |  |    |
  |     | \tiny$\approx \delta'_{Cohen} \times \frac{4df-1}{4(df-1)}$  |   \tiny$\approx \delta'_{Cohen} \times \left( \frac{4df-1}{4(df-1)}-1 \right)$       |\tiny$\approx \frac{df}{df-2} \times \frac{2\left( \frac{\sigma^2_1}{n_1} + \frac{\sigma^2_2}{n_2} \right)}{\sigma^2_1+\sigma^2_2} + (\delta'_{Cohen})^2 \left[ \frac{df}{df-2} - \left( \frac{4 \;df-1}{4(df-1)}\right)^2 \right]$              ||
|    |      |  |    |
| \tiny$Cohen's \; g'_s$    | \tiny$\delta'_{Cohen}$  |   /     |\tiny$Var(Cohen's \; d'_s) \times \left[\frac{\Gamma(\frac{df}{2})}{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}\right]^2$|     |

\elandscape                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

# Table 4
  
  \newpage
\blandscape

| Estimator|      expectancy    |          Bias         |                 Variance                       |  
  |---------|:-----------------:  | :---------------------:  | :-------------------------------------------:  |
  |  \tiny$Shieh's \; d_s$  |\tiny$\delta_{Shieh} \times \frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}$|\tiny$\delta_{Shieh} \times \left(\frac{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}{\Gamma(\frac{df}{2})}-1 \right)$|\tiny$\frac{df}{(df-2)N}  + \delta^2_{Shieh} \left[ \frac{df}{df-2} - \left( \frac{\sqrt{\frac{df}{2}} \times \Gamma\left( \frac{df-1}{2}\right)}{\Gamma\left( \frac{df}{2}\right)} \right)^2 \right]$    |
|    |      |  |    |
|  \tiny$Shieh's \; g_s$  | \tiny$\delta_{Shieh}$  | / |\tiny$Var(Shieh's \; d_s) \times \left[\frac{\Gamma(\frac{df}{2})}{\sqrt{\frac{df}{2}} \times \Gamma(\frac{df-1}{2})}\right]^2$    |

Table: Expentency, bias and variance Shieh's $d_s$ and Shieh's $g_s$ under the assumptions that independent residuals are normally distributed. 

*Note*. $df \approx \frac{\left(\frac{\sigma^2_1}{n_1}+\frac{\sigma^2_2}{n_2} \right)^2}{\frac{(\sigma^2_1/n_1)^2}{n_1-1}+\frac{(\sigma^2_2/n_2)^2}{n_2-1}}$. Equations in Table 4 require $df \ge 3$ and at least 2 subjects per group.

\elandscape

\newpage
