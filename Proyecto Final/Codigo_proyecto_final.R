#Codigó Realizado por Sergio Diaz Vera y Julián Espinosa


# Teoria ------------------------------------------------------------------

# Parámetros
N <- 10
p <- 18/37
q <- 19/37

# Construcción de la matriz P
P <- matrix(0, nrow = N+1, ncol = N+1)

P[1,1] <- 1     # estado 0 absorbente
P[N+1,N+1] <- 1 # estado N absorbente

for(i in 2:N){
  P[i, i-1] <- q
  P[i, i+1] <- p
}

# Distribución inicial: empezar con 5 dólares
lambda <- rep(0, N+1)
lambda[6] <- 1   

# Calcular lambda P^n para n <= 100
num_rounds <- 100
distributions <- matrix(0, nrow = num_rounds, ncol = N+1)

current <- lambda

for(n in 1:num_rounds){
  current <- current %*% P
  distributions[n,] <- current
}



# Graficos ----------------------------------------------------------------
library(ggplot2)
library(reshape2)


df <- as.data.frame(distributions)
df$round <- 1:num_rounds

df.m <- melt(df, id.vars = "round")
df.m$variable <- factor(df.m$variable,
                        labels = paste("Estado", 0:10))
ggplot(df.m, aes(x = round, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  labs(
    title = "Evolución de la distribución del capital (Empezar con 5$)",
    x = "Número de rondas jugadas",
    y = "Probabilidad",
    fill = "Estado"
  ) +
  theme_minimal()



# Simulación Monte Carlo --------------------------------------------------

#librerias
library(dplyr)
library(tidyr)
library(gridExtra)

# Ruleta Europea (37 números: 0-36)
prob_europea <- list(
  color =18/37,          
  docena =12/37,
  seisena =6/37,
  cuadro =4/37,
  semipleno =2/37,
  pleno = 1/37)


# Ruleta Americana (38 números: 0, 00, 1-36)
prob_americana <- list(
  color =18/38,          
  docena =12/38,
  seisena =6/38,
  cuadro =4/38,
  semipleno =2/38,
  pleno = 1/38)


# Pagos según tipo de apuesta
pagos <- list(
  color = 1,           
  docena =2,
  seisena =5,
  cuadro =8,
  semipleno =17,
  pleno =35
)



simular_ruleta <- function(capital_inicial, objetivo, tipo_apuesta, 
                           tipo_ruleta = "europea", n_sim = 100000) {
  
  # Seleccionar probabilidades según tipo de ruleta
  probs <- if(tipo_ruleta == "europea") prob_europea else prob_americana
  p_ganar <- probs[[tipo_apuesta]]
  pago <- pagos[[tipo_apuesta]]
  
  # Vectores para almacenar resultados
  exito <- numeric(n_sim)
  tiempo_absorcion <- numeric(n_sim)
  
  for(sim in 1:n_sim) {
    capital <- capital_inicial
    pasos <- 0
    
    while(capital > 0 && capital < objetivo) {
      # Realizar apuesta de 1 unidad
      if(runif(1) < p_ganar) {
        capital <- capital + pago  # Gana
      } else {
        capital <- capital - 1      # Pierde
      }
      pasos <- pasos + 1
      
      # Límite de seguridad (evitar loops infinitos)
      if(pasos > 100000) break
    }
    
    exito[sim] <- ifelse(capital >= objetivo, 1, 0)
    tiempo_absorcion[sim] <- pasos
  }
  
  return(list(
    prob_exito = mean(exito),
    tiempo_medio = mean(tiempo_absorcion),
    tiempo_sd = sd(tiempo_absorcion),
    todos_tiempos = tiempo_absorcion,
    capital_inicial = capital_inicial,
    objetivo = objetivo,
    tipo_apuesta = tipo_apuesta,
    tipo_ruleta = tipo_ruleta
  ))
}

##Función teoria

calcular_teorico <- function(capital_inicial, objetivo, p_ganar) {
  q <- 1 - p_ganar
  N <- objetivo
  i <- capital_inicial
  
  if(p_ganar == 0.5) {
    # Caso juego justo
    h_i <- (N - i) / N  # Prob. ruina
    k_i <- i * N - i^2  # Tiempo esperado
  } else {
    # Caso juego no justo
    ratio <- q / p_ganar
    h_i <- (ratio^i - ratio^N) / (1 - ratio^N)
    k_i <- i/(q - p_ganar) - N/(q - p_ganar) * (1 - ratio^i)/(1 - ratio^N)
  }
  
  return(list(
    prob_exito = 1 - h_i,
    tiempo_esperado = k_i
  ))
}

## Función matriz transicion

construir_matriz_transicion <- function(capital_inicial, objetivo, p_ganar, pago) {
  N <- objetivo
  n_estados <- N + 1  # Estados de 0 a N
  P <- matrix(0, nrow = n_estados, ncol = n_estados)
  
  # Estado 0 (arruinado)
  P[1, 1] <- 1
  
  # Estado N (objetivo alcanzado)
  P[n_estados, n_estados] <- 1
  
  # Estados intermedios
  for(i in 1:(n_estados - 2)) {
    capital_actual <- i
    # Si gana
    nuevo_capital_ganar <- min(capital_actual + pago, N)
    P[i + 1, nuevo_capital_ganar + 1] <- p_ganar
    
    # Si pierde
    nuevo_capital_perder <- max(capital_actual - 1, 0)
    P[i + 1, nuevo_capital_perder + 1] <- 1 - p_ganar
  }
  
  return(P)
}



analizar_matriz <- function(P, capital_inicial, n_pasos = 100) {
  n_estados <- nrow(P)
  
  # Distribución inicial
  lambda <- rep(0, n_estados)
  lambda[capital_inicial + 1] <- 1  # +1 por indexación en R
  
  # Calcular distribuciones para n pasos
  distribuciones <- matrix(0, nrow = n_pasos + 1, ncol = n_estados)
  distribuciones[1, ] <- lambda
  
  P_n <- P
  for(n in 1:n_pasos) {
    lambda <- lambda %*% P_n
    distribuciones[n + 1, ] <- lambda
    P_n <- P_n %*% P
  }
  
  return(list(
    matriz = P,
    distribuciones = distribuciones,
    prob_absorcion_0 = distribuciones[n_pasos + 1, 1],
    prob_absorcion_N = distribuciones[n_pasos + 1, n_estados]
  ))
}

# ----------------------------------------------------------------------------
# FUNCIÓN: VISUALIZAR MATRIZ DE TRANSICIÓN
# ----------------------------------------------------------------------------

visualizar_matriz <- function(P, titulo = "Matriz de Transición") {
  n <- nrow(P)
  
  # Convertir a formato largo para ggplot
  df_matriz <- expand.grid(
    Estado_Origen = 0:(n-1),
    Estado_Destino = 0:(n-1)
  )
  df_matriz$Probabilidad <- as.vector(t(P))
  
  # Filtrar solo probabilidades > 0 para mejor visualización
  df_matriz <- df_matriz %>% filter(Probabilidad > 0.001)
  
  p <- ggplot(df_matriz, aes(x = Estado_Origen, y = Estado_Destino, 
                             fill = Probabilidad)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "#FFF5E1", high = "#C41E3A", 
                        limits = c(0, 1)) +
    labs(
      title = titulo,
      x = "Estado Actual (Capital)",
      y = "Estado Siguiente (Capital)",
      fill = "Probabilidad"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}



# Simulaciones ------------------------------------------------------------


# Configuraciones a probar
configuraciones <- expand.grid(
  capital_inicial = c(5, 50),
  tipo_apuesta = c("color", "docena","seisena","cuadro","semipleno", "pleno"),
  tipo_ruleta = c("europea", "americana"),
  stringsAsFactors = FALSE
)

# Ejecutar simulaciones
cat("Ejecutando simulaciones...\n")
resultados <- list()

set.seed(11)
for(i in 1:nrow(configuraciones)) {
  config <- configuraciones[i, ]
  objetivo <- config$capital_inicial * 2  # Duplicar capital
  
  cat(sprintf("Simulación %d/%d: %s en ruleta %s (Capital: $%d)\n", 
              i, nrow(configuraciones), 
              config$tipo_apuesta, config$tipo_ruleta, config$capital_inicial))
  
  sim <- simular_ruleta(
    capital_inicial = config$capital_inicial,
    objetivo = objetivo,
    tipo_apuesta = config$tipo_apuesta,
    tipo_ruleta = config$tipo_ruleta,
    n_sim = 100000
  )
  
  resultados[[i]] <- sim
}



# Resultados --------------------------------------------------------------


# Crear dataframe con resultados
df_resultados <- data.frame(
  Capital_Inicial = sapply(resultados, function(x) x$capital_inicial),
  Objetivo = sapply(resultados, function(x) x$objetivo),
  Tipo_Apuesta = sapply(resultados, function(x) x$tipo_apuesta),
  Tipo_Ruleta = sapply(resultados, function(x) x$tipo_ruleta),
  Prob_Exito = sapply(resultados, function(x) x$prob_exito),
  Tiempo_Medio = sapply(resultados, function(x) x$tiempo_medio),
  Tiempo_SD = sapply(resultados, function(x) x$tiempo_sd)
)

# Añadir etiqueta descriptiva
df_resultados$Escenario <- paste0(
  "$", df_resultados$Capital_Inicial, "→$", df_resultados$Objetivo
)

# Imprimir tabla de resultados


print(df_resultados %>% 
        arrange(desc(Prob_Exito)) %>%
        mutate(
          Prob_Exito = sprintf("%.2f%%", Prob_Exito * 100),
          Tiempo_Medio = sprintf("%.1f ± %.1f", Tiempo_Medio, Tiempo_SD)
        ) %>%
        select(Escenario, Tipo_Apuesta, Tipo_Ruleta, Prob_Exito, Tiempo_Medio))

##Grafico 1

p1 <- ggplot(df_resultados, aes(x = Tipo_Apuesta, y = Prob_Exito * 100, 
                                fill = Tipo_Ruleta)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  facet_wrap(~ Escenario) +
  labs(
    title = "Probabilidad de Duplicar el Capital",
    subtitle = "Comparación por tipo de apuesta y ruleta",
    x = "Tipo de Apuesta",
    y = "Probabilidad de Éxito (%)",
    fill = "Tipo de Ruleta"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("europea" = "#2E86AB", "americana" = "#A23B72"))

##Grafico 2

p2 <- ggplot(df_resultados, aes(x = Tipo_Apuesta, y = Tiempo_Medio, 
                                fill = Tipo_Ruleta)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = Tiempo_Medio - Tiempo_SD, ymax = Tiempo_Medio + Tiempo_SD),
    position = position_dodge(0.9),
    width = 0.2
  ) +
  facet_wrap(~ Escenario, scales = "free_y") +
  labs(
    title = "Tiempo Medio hasta Absorción",
    subtitle = "Número promedio de apuestas hasta duplicar capital o arruinarse",
    x = "Tipo de Apuesta",
    y = "Número de Apuestas (pasos)",
    fill = "Tipo de Ruleta"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("europea" = "#2E86AB", "americana" = "#A23B72"))

# Mostrar gráficos
print(p1)
print(p2)


# Simulado vs Teorico -----------------------------------------------------



# Comparar resultados para apuestas al color
comparacion <- data.frame()

for(i in 1:nrow(df_resultados)) {
  if(df_resultados$Tipo_Apuesta[i] == "color") {
    config <- df_resultados[i, ]
    
    probs <- if(config$Tipo_Ruleta == "europea") prob_europea else prob_americana
    p_ganar <- probs[["color"]]
    
    teorico <- calcular_teorico(config$Capital_Inicial, config$Objetivo, p_ganar)
    
    comparacion <- rbind(comparacion, data.frame(
      Escenario = config$Escenario,
      Ruleta = config$Tipo_Ruleta,
      Prob_Simulacion = config$Prob_Exito,
      Prob_Teorica = teorico$prob_exito,
      Tiempo_Simulacion = config$Tiempo_Medio,
      Tiempo_Teorico = teorico$tiempo_esperado,
      Error_Prob = abs(config$Prob_Exito - teorico$prob_exito),
      Error_Tiempo = abs(config$Tiempo_Medio - teorico$tiempo_esperado)
    ))
  }
}

cat("Apuestas al COLOR (Comparación Simulación vs Teoría):\n\n")
print(comparacion %>%
        mutate(
          Prob_Simulacion = sprintf("%.4f", Prob_Simulacion),
          Prob_Teorica = sprintf("%.4f", Prob_Teorica),
          Error_Prob = sprintf("%.4f", Error_Prob),
          Tiempo_Simulacion = sprintf("%.2f", Tiempo_Simulacion),
          Tiempo_Teorico = sprintf("%.2f", Tiempo_Teorico),
          Error_Tiempo = sprintf("%.2f", Error_Tiempo)
        ) %>%
        select(Escenario, Ruleta, Prob_Simulacion, Prob_Teorica, 
               Error_Prob, Tiempo_Simulacion, Tiempo_Teorico, Error_Tiempo))

cat(sprintf("\nError promedio en probabilidad: %.6f (%.4f%%)\n", 
            mean(comparacion$Error_Prob), 
            mean(comparacion$Error_Prob) * 100))
cat(sprintf("Error promedio en tiempo: %.4f pasos\n", 
            mean(comparacion$Error_Tiempo)))




# Conclusiones ------------------------------------------------------------



# Mejor estrategia por probabilidad de éxito
mejor_prob <- df_resultados %>% 
  group_by(Escenario) %>%
  slice_max(Prob_Exito, n = 1)

cat("🎯 MEJOR ESTRATEGIA PARA DUPLICAR CAPITAL:\n\n")
for(i in 1:nrow(mejor_prob)) {
  cat(sprintf("  • %s: Apostar a %s en ruleta %s (%.2f%% éxito, %.0f apuestas)\n",
              mejor_prob$Escenario[i],
              toupper(mejor_prob$Tipo_Apuesta[i]),
              mejor_prob$Tipo_Ruleta[i],
              mejor_prob$Prob_Exito[i] * 100,
              mejor_prob$Tiempo_Medio[i]))
}

# Mayor duración (más apuestas)
mayor_duracion <- df_resultados %>% 
  group_by(Escenario) %>%
  slice_max(Tiempo_Medio, n = 1)

cat("\n⏱️  ESTRATEGIA CON MAYOR DURACIÓN (más apuestas):\n\n")
for(i in 1:nrow(mayor_duracion)) {
  cat(sprintf("  • %s: Apostar a %s en ruleta %s (%.0f apuestas, %.2f%% éxito)\n",
              mayor_duracion$Escenario[i],
              toupper(mayor_duracion$Tipo_Apuesta[i]),
              mayor_duracion$Tipo_Ruleta[i],
              mayor_duracion$Tiempo_Medio[i],
              mayor_duracion$Prob_Exito[i] * 100))
}

# Comparación europea vs americana
cat("\n🌍 COMPARACIÓN RULETA EUROPEA VS AMERICANA:\n\n")
comp <- df_resultados %>%
  group_by(Escenario, Tipo_Apuesta) %>%
  summarise(
    Dif_Prob = diff(Prob_Exito)[1] * 100,
    .groups = "drop"
  )
cat(sprintf("  • La ruleta Americana ofrece en promedio %.2f%% menos probabilidad de éxito\n", 
            mean(comp$Dif_Prob)))


