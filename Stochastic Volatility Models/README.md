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
