### Stochastic Volatility Models - курс "Модели стохастической волатильности"
```lab1.ipynb```

## Описание 

Лабораторная работа состоит из двух тематических частей.

### Часть A — Call-Spread и хеджирование (Black–Scholes)

1. **Ценообразование call-spread.**
   Построены графики цены спрэда, как функции спота $$S\in[0.5,1.5]$$ при $$K_1=1,\ K_2=1.2,\ r=0$$.
   Для трёх волатильностей $\sigma\in{0.1,0.2,0.3}$.
2. **Качественная зависимость цены.**
3. **Репликация и ошибка хеджирования.**
   Реализована дискретная дельта-репликация (шаг $dt=T/n$, $n=123$) в модели БШ.
   Сымитированы $m=2^{10}-1$ траекторий; для каждой траектории посчитана ошибка репликации «портфель − выплата», построена гистограмма распределения ошибки.

4. **Модельный риск и ошибка хеджирования.**
   Пусть реальная динамика $dS_t=\alpha_tS_t,dt+\beta_tS_t,dW_t$.

   * **4.1.** Записан SDE стоимости хеджирующего портфеля $X_t$ при стратегии $H_t=\partial C/\partial S$.
   * **4.2.** По Ито получено SDE для $Y_t=C(t,S_t)$.
   * **4.3.** Показано, что ошибка $Z_t=X_t-Y_t$ удовлетворяет
     $$\frac{d}{dt}Z_t=rZ_t+\tfrac12 S_t^2\frac{\partial^2 C}{\partial S^2}\big(\hat\sigma^2-\beta_t^2\big)$$
   * **4.4.** Решение:
     $$Z_T=\int_0^T e^{r(T-s)}\tfrac12 S_s^2\frac{\partial^2 C}{\partial S^2}(s,S_s)\big(\hat\sigma^2-\beta_s^2\big)ds$$
   * **4.5.** Обсуждена возможность гарантированной прибыли при систематическом смещении $\hat\sigma^2-\beta_t^2$ c фиксированным знаком (и $\Gamma\ge 0$).

### Часть B — Имплайд- и локальная волатильности, Монте-Карло

1. **Построение IV-поверхности.**
   Из `data.txt` извлечены ($T_i,K_j,C_{ij}$). Для каждого узла вычислены $\sigma_{IV}(T_i,K_j)$ (устойчивый инвертор через граничное отсечение).
   Интер-/экстраполяция выполнена на поверхности полной дисперсии $w=\sigma^2 T$ в координатах $(T,y)$, $y=\log(K/F_T)$ (у нас $F_T=1$).
   Выведены:

   * **кривые IV по страйку** для каждого $T_i$;
   * **3D-поверхность IV** $\sigma_{IV}(T,K)$.

2. **Локальная волатильность по Дюпиру.**
   Реализована функция $\sigma_{\text{Dup}}(t,s)$ с зажимом аргументов в диапазоне данных и числовыми стабилизаторами.
   Построена **поверхность локальной волатильности** (\sigma_{\text{Dup}}(t,s)).

3. **Монте-Карло-калькулятор (T=0.25).**
   Под локальной волатильностью (лог-Эйлер, антидетические вариации) оценены цены $V(T,K)$ для $K={0.6,0.7,\dots,1.4}$ и соответствующие $\sigma_{IV}(T,K)$.
   Точность: абсолютная ошибка IV $\le 1%$ для «почти всех» страйков (при необходимости увеличены paths/steps).
   Построены:

   * **цены MC vs файл** на одном графике;
   * **IV MC vs IV из пункта (1)** на другом.

## Description

The lab contains two blocks.

### Part A — Call-Spread & Hedging (Black–Scholes)

1. **Call-spread pricing.**
   Plotted the price of call-spread.

2. **Qualitative dependence.**
   Discussed how the price depends on spot (S), volatility ($\sigma$) and time (T) (monotonicity, convexity, time value).

3. **Replication & hedging error.**
   Discrete delta-hedging with step $dt=T/123$; simulated $m=2^{10}-1$ paths; computed terminal replication error and plotted its histogram.

4. **Model risk: general Itô dynamics.**
   With $dS_t=\alpha_tS_tdt+\beta_tS_tdW_t$:

   * **4.1** SDE for portfolio value $X_t$ under $H_t=\partial C/\partial S$.
   * **4.2** Itô SDE for $Y_t=C(t,S_t)$.
   * **4.3** Derived
     $\tfrac{d}{dt}Z_t=rZ_t+\tfrac12 S_t^2 C_{SS}(\hat\sigma^2-\beta_t^2),
     \quad Z_t=X_t-Y_t$
   * **4.4** Solution
     $Z_T=\int_0^T e^{r(T-s)}\tfrac12 S_s^2 C_{SS}(s,S_s)(\hat\sigma^2-\beta_s^2)ds$
   * **4.5** Sufficient condition for P-a.s. profit: persistent sign of $\hat\sigma^2-\beta_t^2$ with non-negative gamma.

### Part B — Implied & Local Vol Surfaces, Monte Carlo

1. **Implied-vol surface.**
   From `data.txt` extracted $T_i,K_j,C_{ij}$; computed $\sigma_{IV}(T_i,K_j)$ with a robust BS inversion.
   Inter-/extrapolated **total variance**.
   Plotted **smiles by maturity** and the **3D IV surface**.

2. **Dupire local volatility (20 pts).**
   Implemented $\sigma_{\text{Dup}}(t,s)$ with clamping and numerical guards; plotted the **local-vol surface** on a $(t,s)$ grid.

3. **Monte-Carlo pricing at (T=0.25).**
   Plotted **MC vs file prices** and **MC vs question-(1) IVs**.


```lab2.ipynb```
# Heston Model & Calibration

This notebook contains a complete workflow for pricing and calibrating the **Heston stochastic volatility model** using the **COS (Fourier–cosine) method**, including truncation-range design, implied-vol conversion, calibration to synthetic market data, calibration with bid–ask spreads, and “smile + response” calibration across two dates.

---

## 1) Model setup (Heston)

Under the risk-neutral measure:
$$\frac{dS_t}{S_t} = r\,dt + \sqrt{V_t}\,dW_t^{(1)},\qquad
dV_t = \kappa(\theta - V_t)\,dt + \sigma\sqrt{V_t}\,dW_t^{(2)},\qquad
d\langle W^{(1)},W^{(2)}\rangle_t = \rho\,dt$$

Parameters:
- $v_0$: initial variance
- $\kappa$: mean-reversion speed
- $\theta$: long-run variance
- $\sigma$: vol-of-vol
- $\rho$: correlation

We implemented:
- `MarketState`, `HestonParams`, `StockOption` (call/put)
- Heston **log-characteristic function** $\log \phi_T(\omega)$

---

## 2) Pricing with COS method

For a European payoff $g(S_T)$:
$$C = e^{-rT}\,\mathbb{E}[g(S_T)]$$

COS approximation (Fang–Oosterlee):
$$C\approx e^{-rT} \bigl( \frac{V_0}{2}+\sum_{n=1}^{N}\Re \{ 
\phi_T\ ( \frac{n\pi}{b-a})\exp\ ( -in\pi\frac{x-a}{b-a} ) 
\} V_n \bigr)$$,
where $x=\log(S_0/K)$ and $V_n$ are payoff cosine coefficients (implemented analytically for vanilla put).
Call prices are obtained via put–call parity:
$$C_{\text{call}} = C_{\text{put}} + S_0 - Ke^{-rT}$$

---

## 3) Truncation range $[a,b]$ (Task 1)

COS requires truncating the integration domain to $[a,b]$. We tested cumulant-based truncation rules.

We used finite differences to approximate derivatives at $\omega=0$ and applied the Fang–Oosterlee rule:
$$[a,b]=\left[c_1 - L\sqrt{c_2+\sqrt{c_4}},\; c_1 + L\sqrt{c_2+\sqrt{c_4}}\right],\qquad L\approx 10$$

**Short-maturity test:** prices for calls and puts were plotted vs strike for $T\in\{7d,1m,2m,3m\}$, comparing:
- standard baseline (e.g. $c_2$-only rule),
- cumulant rule using $c_2$ and $c_4$.

---

## 4) From price to Black implied volatility (IV)

We compute forward and discount factors:
$$F = S_0 e^{rT},\qquad DF=e^{-rT},\qquad C^{\text{undisc}} = \frac{C}{DF}$$

For a call, no-arbitrage bounds in forward measure:
$$(F-K)^+ \le C^{\text{undisc}} \le F$$

Numerical COS pricing can slightly violate the lower bound, causing IV inversion failures.  
We stabilized IV inversion by clipping:
$$C^{\text{undisc}} \leftarrow \min\Big(\max(C^{\text{undisc}}, (F-K)^+ + \epsilon),\; F-\epsilon\Big),\qquad \epsilon\ll 1$$

---

## 5) Calibration to IV surface (Task 2)

Calibration objective (mean squared error in IV):
$$\mathrm{Obj}=
\frac{1}{|\mathcal{T}|}\sum_{T\in\mathcal{T}}
\frac{1}{|\mathcal{K}|}\sum_{K\in\mathcal{K}}
\big(IV_{\text{data}}(T,K)-IV_{\text{model}}(T,K)\big)^2$$

Implementation details:
- Model IV is produced by **COS pricing → IV inversion** on the same $(T,K)$ grid.
- Optimization: bounded `L-BFGS-B` with multiple starting points.
- Optional **Feller regularization**:
  $$2\kappa\theta \ge \sigma^2,\quad
  \text{penalty}=\alpha\max(0,\sigma^2-2\kappa\theta)^2$$

Plots:
- market IV points vs calibrated model IV curve
- for smooth visuals, model curves were also evaluated on a denser strike grid.

---

## 6) Calibration with bid–ask spread (Task 3)

Market provides $IV_{\text{bid}}(T,K)$ and $IV_{\text{ask}}(T,K)$.  
Objective:
$$\mathrm{Obj}_{BA}=
\frac{1}{|\mathcal{T}|}\sum_{T\in\mathcal{T}}
\frac{1}{|\mathcal{K}|}\sum_{K\in\mathcal{K}}
\frac{1}{2}\Big[
(IV_{\text{bid}}-IV_{\text{model}})^2 + (IV_{\text{ask}}-IV_{\text{model}})^2
\Big]$$

Result visualization:
- bid and ask IV points (different markers)
- calibrated model IV curve lying approximately between bid and ask.

---

## 7) Smile + response calibration (Task 4)

We observe one maturity $T$ but two dates:
- at $t$: spot $S_t$ and IV smile $IV_t(K)$
- at $t+1$: spot $S_{t+1}$ and smile $IV_{t+1}(K)$

We calibrate **one common parameter set** to satisfy:
1) fit the smile at $t$
2) match the **IV response** due to spot move $S_t\to S_{t+1}$

Define market and model responses:
$$\Delta IV_{\text{mkt}}(K)=IV_{t+1}(K)-IV_t(K),
\qquad
\Delta IV_{\text{model}}(K)=IV_{\text{model}}(K;S_{t+1})-IV_{\text{model}}(K;S_t)$$

Outputs:
- plot of $IV_t(K)$ and $IV_{t+1}(K)$ with corresponding model curves
- plot of $\Delta IV(K)$ (market vs model)

Interpretation:
- Heston typically yields a **smooth** $\Delta IV(K)$ response.
- Market $\Delta IV(K)$ can be more strike-dependent; $\lambda$ controls the trade-off between smile fit and response fit.

---

## Notes / Troubleshooting

- **IV inversion errors (BelowIntrinsic):** fixed by enforcing call bounds $(F-K)^+\le C^{undisc}\le F$.
- **COS divide-by-zero / NaNs during calibration:** prevented by enforcing a minimum truncation width.

---


