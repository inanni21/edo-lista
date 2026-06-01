#include <stdio.h>
#include <math.h>

/* ============================================================
   FUNÇÕES DO LADO DIREITO DAS EDOs  f(t, y) = y'
   ============================================================ */

/* Questão 1: ty' + 2y = t^2  =>  y' = (t^2 - 2y) / t */
double f1(double t, double y) {
    return (t * t - 2.0 * y) / t;
}

/* Solução exata Q1: y(t) = (t^4 + 1) / (4t^2) */
double exata1(double t) {
    return (pow(t, 4) + 1.0) / (4.0 * t * t);
}

/* Questão 2: y' + y^2*sen(x) = 0  =>  y' = -y^2 * sen(x) */
double f2(double x, double y) {
    return -y * y * sin(x);
}

/* Solução exata Q2: y(x) = 1 / (2 - cos(x)) */
double exata2(double x) {
    return 1.0 / (2.0 - cos(x));
}

/* Questão 3: y' + y^2 = t  =>  y' = t - y^2 */
double f3(double t, double y) {
    return t - y * y;
}

/* Questão 4: dy/dt = k*y  =>  y' = k*y */
double f4(double t, double y) {
    double k = 0.438203;
    return k * y;
}

/* Solução exata Q4: y(t) = y0 * e^(kt) */
double exata4(double t) {
    double k = 0.438203;
    double y0 = 3.0;
    return y0 * exp(k * t);
}

/* ============================================================
   MÉTODOS NUMÉRICOS
   ============================================================ */

/* Euler explícito
   t0, y0: condição inicial
   tf    : tempo final
   n     : número de passos
   f     : função f(t, y)
   Imprime a tabela de resultados                              */
void euler(double t0, double y0, double tf, int n,
           double (*f)(double, double),
           double (*exata)(double)) {

    double h = (tf - t0) / n;
    double t = t0, y = y0;

    printf("  i |     t      |   y (Euler)  |  y (exata)   |   erro abs\n");
    printf("----+------------+--------------+--------------+------------\n");

    for (int i = 0; i <= n; i++) {
        if (exata != NULL)
            printf(" %2d | %10.6f | %12.8f | %12.8f | %12.2e\n",
                   i, t, y, exata(t), fabs(y - exata(t)));
        else
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);

        if (i < n) {
            y = y + h * f(t, y);
            t = t + h;
        }
    }
}

/* Runge-Kutta de ordem 2 (método do ponto médio) */
void rk2(double t0, double y0, double tf, int n,
         double (*f)(double, double),
         double (*exata)(double)) {

    double h = (tf - t0) / n;
    double t = t0, y = y0;

    printf("  i |     t      |    y (RK2)   |  y (exata)   |   erro abs\n");
    printf("----+------------+--------------+--------------+------------\n");

    for (int i = 0; i <= n; i++) {
        if (exata != NULL)
            printf(" %2d | %10.6f | %12.8f | %12.8f | %12.2e\n",
                   i, t, y, exata(t), fabs(y - exata(t)));
        else
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);

        if (i < n) {
            double k1 = f(t, y);
            double k2 = f(t + h / 2.0, y + (h / 2.0) * k1);
            y = y + h * k2;
            t = t + h;
        }
    }
}

/* Runge-Kutta de ordem 4 (clássico) */
void rk4(double t0, double y0, double tf, int n,
         double (*f)(double, double),
         double (*exata)(double)) {

    double h = (tf - t0) / n;
    double t = t0, y = y0;

    printf("  i |     t      |    y (RK4)   |  y (exata)   |   erro abs\n");
    printf("----+------------+--------------+--------------+------------\n");

    for (int i = 0; i <= n; i++) {
        if (exata != NULL)
            printf(" %2d | %10.6f | %12.8f | %12.8f | %12.2e\n",
                   i, t, y, exata(t), fabs(y - exata(t)));
        else
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);

        if (i < n) {
            double k1 = f(t,           y);
            double k2 = f(t + h / 2.0, y + (h / 2.0) * k1);
            double k3 = f(t + h / 2.0, y + (h / 2.0) * k2);
            double k4 = f(t + h,       y + h * k3);
            y = y + (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            t = t + h;
        }
    }
}

/* ============================================================
   MAIN
   ============================================================ */
int main(void) {

    /* ----------------------------------------------------------
       QUESTÃO 1
       ty' + 2y = t^2,  y(1) = 0.5,  t em [1, 10],  n = 10
       Solução exata: y(t) = (t^4 + 1) / (4t^2)
       ---------------------------------------------------------- */
    printf("==========================================================\n");
    printf("  QUESTÃO 1: ty' + 2y = t^2,  y(1) = 0.5,  t in [1,10]\n");
    printf("==========================================================\n\n");

    printf("--- (a) Método de Euler, n = 10 ---\n");
    euler(1.0, 0.5, 10.0, 10, f1, exata1);

    printf("\n--- (b) Runge-Kutta Ordem 2, n = 10 ---\n");
    rk2(1.0, 0.5, 10.0, 10, f1, exata1);

    printf("\n--- (c) Runge-Kutta Ordem 4, n = 10 ---\n");
    rk4(1.0, 0.5, 10.0, 10, f1, exata1);

    /* ----------------------------------------------------------
       QUESTÃO 2
       y' + y^2*sen(x) = 0,  y(0) = 1,  x em [0, 3],  n = 10
       Solução exata: y(x) = 1 / (2 - cos(x))
       ---------------------------------------------------------- */
    printf("\n==========================================================\n");
    printf("  QUESTÃO 2: y' + y^2*sen(x) = 0,  y(0)=1,  x in [0,3]\n");
    printf("==========================================================\n\n");

    printf("--- (a) Método de Euler, n = 10 ---\n");
    euler(0.0, 1.0, 3.0, 10, f2, exata2);

    printf("\n--- (b) Runge-Kutta Ordem 2, n = 10 ---\n");
    rk2(0.0, 1.0, 3.0, 10, f2, exata2);

    printf("\n--- (c) Runge-Kutta Ordem 4, n = 10 ---\n");
    rk4(0.0, 1.0, 3.0, 10, f2, exata2);

    /* ----------------------------------------------------------
       QUESTÃO 3
       y' + y^2 = t,  y(0) = 0,  t em [0, 8],  h = 0.1  =>  n = 80
       Sem solução exata fechada simples (não é fornecida)
       ---------------------------------------------------------- */
    printf("\n==========================================================\n");
    printf("  QUESTÃO 3: y' + y^2 = t,  y(0) = 0,  t in [0,8], h=0.1\n");
    printf("==========================================================\n\n");

    int n3 = 80;   /* h = 0.1 => n = (8 - 0) / 0.1 = 80 */

    printf("--- (a) Método de Euler, h = 0.1 (n = 80) ---\n");
    printf("  i |     t      |   y (Euler)\n");
    printf("----+------------+------------\n");
    {
        double h = 0.1, t = 0.0, y = 0.0;
        for (int i = 0; i <= n3; i++) {
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);
            if (i < n3) { y = y + h * f3(t, y); t += h; }
        }
    }

    printf("\n--- (b) Runge-Kutta Ordem 2, h = 0.1 (n = 80) ---\n");
    printf("  i |     t      |    y (RK2)\n");
    printf("----+------------+------------\n");
    {
        double h = 0.1, t = 0.0, y = 0.0;
        for (int i = 0; i <= n3; i++) {
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);
            if (i < n3) {
                double k1 = f3(t, y);
                double k2 = f3(t + h/2.0, y + (h/2.0)*k1);
                y = y + h * k2; t += h;
            }
        }
    }

    printf("\n--- (c) Runge-Kutta Ordem 4, h = 0.1 (n = 80) ---\n");
    printf("  i |     t      |    y (RK4)\n");
    printf("----+------------+------------\n");
    {
        double h = 0.1, t = 0.0, y = 0.0;
        for (int i = 0; i <= n3; i++) {
            printf(" %2d | %10.6f | %12.8f\n", i, t, y);
            if (i < n3) {
                double k1 = f4(t, y);   /* atenção: f3 para Q3 */
                double k1_ = f3(t, y);
                double k2_ = f3(t + h/2.0, y + (h/2.0)*k1_);
                double k3_ = f3(t + h/2.0, y + (h/2.0)*k2_);
                double k4_ = f3(t + h,     y + h*k3_);
                y = y + (h/6.0)*(k1_ + 2.0*k2_ + 2.0*k3_ + k4_);
                t += h;
                (void)k1; /* silencia warning variável não usada */
            }
        }
    }

    /* ----------------------------------------------------------
       QUESTÃO 4
       dy/dt = k*y,  y(0) = 3,  k = 0.438203,  t em [0, 30],  n = 20
       Solução exata: y(t) = 3 * e^(k*t)
       ---------------------------------------------------------- */
    printf("\n==========================================================\n");
    printf("  QUESTÃO 4: dy/dt = k*y,  y(0)=3,  k=0.438203, t in [0,30]\n");
    printf("==========================================================\n\n");

    printf("--- Euler, n = 20 ---\n");
    euler(0.0, 3.0, 30.0, 20, f4, exata4);

    printf("\n--- Runge-Kutta Ordem 2, n = 20 ---\n");
    rk2(0.0, 3.0, 30.0, 20, f4, exata4);

    printf("\n--- Runge-Kutta Ordem 4, n = 20 ---\n");
    rk4(0.0, 3.0, 30.0, 20, f4, exata4);

    return 0;
}