#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100

/* Geometría del bosque */
#define TAM_CELDA_M 5.0 // puede ser 10.0

/* Temperaturas en grados Celsius */
#define T_AMBIENTE_C 25.0
#define T_MAXIMA_C 390.0
#define T_FINAL_C 40.0
#define DELTA_T_MINIMO_C 0.1

/* Duración de la curva de un árbol */
#define TIEMPO_MAXIMO_S 250.0
#define TIEMPO_TOTAL_S 500.0

/* Temperatura de ignición en Kelvin y Celsius */
#define T_IGNICION_K 560.9278
#define T_IGNICION_C (T_IGNICION_K - 273.15)

#define PI 3.14159265358979323846
#define SIGMA 5.670400e-8
#define RADIO_ARBOL_M 2.0

/*
 * PARÁMETROS AJUSTADOS PARA PROPAGACIÓN:
 * - MASA_COPA_KG: Reducida drásticamente a 2 kg (hojas secas).
 * - ALPHA: Reducido a 0.05 para que la radiación llegue más lejos.
 * - DT: Aumentado a 20 segundos para acumular energía.
 */
#define MASA_COPA_KG 2.0
#define CALOR_ESPECIFICO_J_KG_K 1700.0
#define ALPHA 0.05
#define DT 20.0

typedef enum {
    VACIO,
    FRIO,
    CALENTANDO,
    ARDIENDO
} Estado;

typedef struct {
    int hay_arbol;
    Estado estado;
    double temperatura_c;
    int tiempo_ignicion_s;
    double temperatura_siguiente_c;
} Arbol;

double temperatura_arbol(double tiempo_s) {
    if (tiempo_s < 0.0) tiempo_s = 0.0;
    if (tiempo_s <= TIEMPO_MAXIMO_S) {
        double argumento = PI * tiempo_s / (2.0 * TIEMPO_MAXIMO_S);
        return T_IGNICION_C + (T_MAXIMA_C - T_IGNICION_C) * sin(argumento);
    }
    double proporcion_final = (T_FINAL_C - T_AMBIENTE_C) / (T_MAXIMA_C - T_AMBIENTE_C);
    double k = -log(proporcion_final) / TIEMPO_MAXIMO_S;
    double tiempo_enfriamiento = tiempo_s - TIEMPO_MAXIMO_S;
    return T_AMBIENTE_C + (T_MAXIMA_C - T_AMBIENTE_C) * exp(-k * tiempo_enfriamiento);
}

double radiacion_emitida(double temperatura_c) {
    double temperatura_k = temperatura_c + 273.15;
    return SIGMA * pow(temperatura_k, 4.0);
}

double radiacion_recibida(double temperatura_a_c, double distancia_m, double alpha) {
    if (distancia_m <= 0.0 || alpha < 0.0) return 0.0;
    double temperatura_k = temperatura_a_c + 273.15;
    double flujo_emitido = SIGMA * pow(temperatura_k, 4.0);
    double factor_geometrico = (RADIO_ARBOL_M * RADIO_ARBOL_M) / (distancia_m * distancia_m);
    double factor_atenuacion = exp(-alpha * distancia_m);
    return flujo_emitido * factor_geometrico * factor_atenuacion;
}

double potencia_interceptada(double temperatura_a_c, double distancia_m, double alpha) {
    double flujo = radiacion_recibida(temperatura_a_c, distancia_m, alpha);
    double area_frontal = PI * RADIO_ARBOL_M * RADIO_ARBOL_M;
    return flujo * area_frontal;
}

double energia_absorbida(double temperatura_a_c, double distancia_m, double alpha, double intervalo_s) {
    if (intervalo_s < 0.0) return 0.0;
    double potencia = potencia_interceptada(temperatura_a_c, distancia_m, alpha);
    return potencia * intervalo_s;
}

double distancia_arboles(int fila_a, int columna_a, int fila_b, int columna_b) {
    double df = (double)(fila_b - fila_a);
    double dc = (double)(columna_b - columna_a);
    return TAM_CELDA_M * sqrt(df * df + dc * dc);
}

double actualizar_temperatura(double temperatura_anterior_c, double energia_j, double masa_kg, double calor_especifico_j_kg_k) {
    if (energia_j < 0.0 || masa_kg <= 0.0 || calor_especifico_j_kg_k <= 0.0) return temperatura_anterior_c;
    double aumento_c = energia_j / (masa_kg * calor_especifico_j_kg_k);
    return temperatura_anterior_c + aumento_c;
}

void inicializar_bosque(Arbol bosque[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bosque[i][j].hay_arbol = 0;
            bosque[i][j].estado = VACIO;
            bosque[i][j].temperatura_c = T_AMBIENTE_C;
            bosque[i][j].tiempo_ignicion_s = -1;
            bosque[i][j].temperatura_siguiente_c = T_AMBIENTE_C;
        }
    }
}

void generar_vegetacion(Arbol bosque[N][N]) {
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            double aleatorio = (double)rand() / ((double)RAND_MAX + 1.0);
            if (aleatorio > 0.40) {
                bosque[i][j].hay_arbol = 1;
                bosque[i][j].estado = FRIO;
                bosque[i][j].temperatura_c = T_AMBIENTE_C;
                bosque[i][j].temperatura_siguiente_c = T_AMBIENTE_C;
            }
        }
    }
}

int contar_arboles(Arbol bosque[N][N]) {
    int total = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (bosque[i][j].hay_arbol) total++;
    return total;
}

void iniciar_incendio(Arbol bosque[N][N], int fila, int columna) {
    if (!bosque[fila][columna].hay_arbol) bosque[fila][columna].hay_arbol = 1;
    bosque[fila][columna].estado = ARDIENDO;
    bosque[fila][columna].temperatura_c = T_IGNICION_C;
    bosque[fila][columna].temperatura_siguiente_c = T_IGNICION_C;
    bosque[fila][columna].tiempo_ignicion_s = 0;
}

void calcular_energia_recibida(Arbol bosque[N][N], double energia[N][N], int tiempo_s) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            energia[i][j] = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (bosque[i][j].estado != ARDIENDO || bosque[i][j].tiempo_ignicion_s < 0) continue;

            double edad = (double)(tiempo_s * DT) - (double)(bosque[i][j].tiempo_ignicion_s * DT);
            double temperatura_emisor = temperatura_arbol(edad);

            // Aumentar el radio de influencia di,dj ∈ [-2,2] → radio 2 celdas = 20 
            // Amplía a [-4,4] para que el frente de fuego tenga más alcance.
            for (int di = -4; di <= 4; di++) {          // Cambiar en este renglón los 2 por 4 
                for (int dj = -4; dj <= 4; dj++) {      // Cambiar en este renglón los 2 por 4
                    if (di == 0 && dj == 0) continue;
                    int ni = i + di;
                    int nj = j + dj;
                    if (ni < 0 || ni >= N || nj < 0 || nj >= N) continue;
                    if (!bosque[ni][nj].hay_arbol || bosque[ni][nj].estado == ARDIENDO) continue;

                    double distancia = distancia_arboles(i, j, ni, nj);
                    double energia_j = energia_absorbida(temperatura_emisor, distancia, ALPHA, DT);
                    energia[ni][nj] += energia_j;
                }
            }
        }
    }
}

void calcular_temperaturas_siguientes(Arbol bosque[N][N], double energia[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Arbol *arbol = &bosque[i][j];
            arbol->temperatura_siguiente_c = arbol->temperatura_c;
            if (!arbol->hay_arbol || arbol->estado == ARDIENDO) continue;

            arbol->temperatura_siguiente_c = actualizar_temperatura(
                arbol->temperatura_c, energia[i][j], MASA_COPA_KG, CALOR_ESPECIFICO_J_KG_K);
        }
    }
}

int actualizar_estados_ignicion(Arbol bosque[N][N], int tiempo_s) {
    int nuevas_igniciones = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Arbol *arbol = &bosque[i][j];
            if (!arbol->hay_arbol || arbol->estado == ARDIENDO) continue;

            arbol->temperatura_c = arbol->temperatura_siguiente_c;

            if (arbol->temperatura_c >= T_IGNICION_C) {
                arbol->estado = ARDIENDO;
                arbol->tiempo_ignicion_s = tiempo_s;
                nuevas_igniciones++;
            } else if (arbol->temperatura_c > T_AMBIENTE_C + DELTA_T_MINIMO_C) {
                arbol->estado = CALENTANDO;
            } else {
                arbol->estado = FRIO;
            }
        }
    }
    return nuevas_igniciones;
}

int simular_paso(Arbol bosque[N][N], int tiempo_s, double energia[N][N]) {
    calcular_energia_recibida(bosque, energia, tiempo_s);
    calcular_temperaturas_siguientes(bosque, energia);
    return actualizar_estados_ignicion(bosque, tiempo_s);
}

void contar_estados(Arbol bosque[N][N], int *vacios, int *frios, int *calentando, int *ardiendo) {
    *vacios = 0; *frios = 0; *calentando = 0; *ardiendo = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (!bosque[i][j].hay_arbol) { (*vacios)++; continue; }
            switch (bosque[i][j].estado) {
                case FRIO: (*frios)++; break;
                case CALENTANDO: (*calentando)++; break;
                case ARDIENDO: (*ardiendo)++; break;
                default: (*vacios)++; break;
            }
        }
    }
}

void guardar_resultados(FILE *archivo, Arbol bosque[N][N], int tiempo_s) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(archivo, "%d,%d,%d,%d,%.6f,%d\n",
                tiempo_s, i, j, bosque[i][j].estado,
                bosque[i][j].temperatura_c, bosque[i][j].tiempo_ignicion_s);
        }
    }
}

void guardar_resumen(FILE *archivo, Arbol bosque[N][N], int tiempo_s) {
    int vacios, frios, calentando, ardiendo;
    contar_estados(bosque, &vacios, &frios, &calentando, &ardiendo);
    fprintf(archivo, "%.1f,%d,%d,%d,%d\n", tiempo_s * DT, vacios, frios, calentando, ardiendo);
}

int main(void) {
    static Arbol bosque[N][N];
    static double energia[N][N];

    srand(42);
    inicializar_bosque(bosque);
    generar_vegetacion(bosque);

    const int fila_inicial = 43;
    const int columna_inicial = 51;
    iniciar_incendio(bosque, fila_inicial, columna_inicial);

    printf("Bosque inicializado. Arboles: %d\n", contar_arboles(bosque));
    printf("Temperatura de ignicion: %.2f C\n", T_IGNICION_C);
    printf("Masa de la copa: %.2f kg\n", MASA_COPA_KG);
    printf("Paso de tiempo (DT): %.1f s\n", DT);

    FILE *archivo = fopen("resultados_incendio.csv", "w");
    FILE *resumen = fopen("resumen_incendio.csv", "w");
    if (!archivo || !resumen) { perror("Error al crear CSV"); return 1; }

    fprintf(archivo, "tiempo_s,fila,columna,estado,temperatura_c,tiempo_ignicion_s\n");
    fprintf(resumen, "tiempo_s,vacios,frios,calentando,ardiendo\n");

    guardar_resultados(archivo, bosque, 0);
    guardar_resumen(resumen, bosque, 0);

    int total_igniciones = 0;
    int iteraciones = (int)(TIEMPO_TOTAL_S / DT);

    printf("\nEjecutando simulacion...\n");

    for (int paso = 1; paso <= iteraciones; paso++) {
        int nuevas_igniciones = simular_paso(bosque, paso, energia);
        total_igniciones += nuevas_igniciones;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (bosque[i][j].estado == ARDIENDO) {
                    int edad = (paso - bosque[i][j].tiempo_ignicion_s) * DT;
                    bosque[i][j].temperatura_c = temperatura_arbol((double)edad);
                }
            }
        }

        guardar_resultados(archivo, bosque, paso);
        guardar_resumen(resumen, bosque, paso);

        if (paso % 5 == 0) {
            printf("Paso %d (Tiempo real: %.0f s). Nuevas igniciones: %d\n", paso, paso * DT, nuevas_igniciones);
        }
    }

    fclose(archivo);
    fclose(resumen);
    printf("\nSimulacion terminada. Total de igniciones: %d\n", total_igniciones);
    return 0;
}