# Cargar las librerías necesarias
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpmisc)

setwd("/vast/palmer/pi/powell/DB_upload/BCF_MASC/LATB/VCFtools/LD_LARE//")
rm(list=ls())

# Definir los valores de prefix
prefix_values <- c("CHM","DGT","MAG","YEM")
prefix_values <- c("L00","TB0")


# Iterar sobre cada prefijo
for (prefix in prefix_values) {
  # Definir el nombre del archivo de entrada
  input_file <- paste0(prefix, ".geno.ld.clean.gz")
  
  # Verificar si el archivo existe antes de procesarlo
  if (file.exists(input_file)) {
    
    # Leer los datos usando data.table
    data <- fread(input_file, header = T, sep = "\t", dec = ".",  col.names = c("CHR", "POS1", "POS2", "N_INDV", "R^2"))
    
    # Filtrar por cromosomas específicos
    data <- data %>%
      filter(CHR %in% c("scaffold_1", "scaffold_2", "scaffold_3"))
    
    # Calcular la distancia entre POS2 y POS1
    data$Dist <- data$POS2 - data$POS1
    
    # Crear una variable de distancia categorizada en intervalos de 1000 unidades
    data$Dist_cat <- cut(data$Dist, breaks = seq(from = min(data$Dist) - 1, to = max(data$Dist) + 1, by = 500))
    
    # Calcular las medias y medianas de R^2 para cada intervalo de distancia
    data_summarized <- data %>%
      group_by(Dist_cat) %>%
      summarise(mean_R2 = mean(`R^2`, na.rm = TRUE), 
                median_R2 = median(`R^2`, na.rm = TRUE)) %>%
      mutate(start = as.numeric(str_extract(str_replace_all(Dist_cat, "[\\(\\)\\[\\]]", ""), "^[0-9-e+.]+")),
             end = as.numeric(str_extract(str_replace_all(Dist_cat, "[\\(\\)\\[\\]]", ""), "[0-9-e+.]+$")),
             mid = start + ((end - start) / 2))
    
    # Definir el nombre del archivo de salida
    output_file <- paste0(prefix, "_resumen.csv")
    
    # Guardar los resultados en un archivo CSV
    fwrite(data_summarized, file = output_file, sep = ",", dec = ".", quote = FALSE)
    
    # Mensaje indicando que el archivo ha sido procesado
    cat("Archivo procesado:", input_file, "->", output_file, "\n")
  } else {
    cat("Archivo no encontrado:", input_file, "\n")
  }
}

#############################################################################################################################################
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

rm(list=ls())
# Definir los valores de prefix y cromosoma
prefix_values <- c("L00","TB0") #,"JHM","MAD","MEC","NAJ","SAH")
prefix_values <- c("CHM","DGT","MAG", "YEM")

# Inicializar un dataframe vacío para almacenar todos los datos concatenados
final_data <- data.frame()

# Iterar sobre cada combinación de prefix y cromosoma
for (prefix in prefix_values) {
  
  # Crear el nombre del archivo basado en el prefix y el cromosoma
  input_file <- paste0(prefix, "_resumen.csv")
  
  # Leer el archivo procesado si existe
  if (file.exists(input_file)) {
    data <- read.csv(input_file, header = TRUE, sep = ",")
    data$Source <- prefix
    # Concatenar al dataframe final
    final_data <- rbind(final_data, data)
  }
}

# Eliminar todas las filas con cualquier valor NA
final_data_clean <- na.omit(final_data)

n_map <- c(
  "CHM" = 5,
  "DGT" = 5,
  "MAG" = 5,
  "YEM" = 3,
  "L00" = 7,
  "TB0" = 15
)

final_data_clean <- final_data_clean %>%
  mutate(
    n = n_map[Source],
    mean_R2_adj = mean_R2 - ((1 - mean_R2) / (n - 2))
  )

ggplot(data = final_data_clean, aes(x = mid, y = mean_R2_adj, colour = Source)) +
  geom_point(size = 3, shape = 21, fill = NA) +  # Relleno vacío y contorno de color
  labs(x = "Distance (Kb)", y = expression(LD ~ (r^{2})), colour = "Population") +  # Cambiar el título de la leyenda
  scale_x_continuous(limits = c(0, 50000), breaks = c(0, 1*10^4, 2*10^4, 3*10^4, 4*10^4, 5*10^4),
                     labels = c("0", "10", "20", "30", "40", "50")) +
  facet_wrap(. ~ Source, ncol=4, nrow = 4) +  # Usar facet_grid para dividir por cromosoma
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.8),
        axis.line = element_line(size = 0.8, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",  # Colocar las etiquetas de las facetas fuera del panel
        strip.text = element_text(size = 8)) + # Ajustar el tamaño de las etiquetas
  geom_smooth(method = "nls", 
              formula = y ~ SSasymp(x, Asym, R0, lrc), color = "black",
              se = FALSE, aes(group = Source)) +
  guides(colour = guide_legend(override.aes = list(size = 4)))  # Agrandar los puntos en la leyenda

ggsave("LD_distance_LATB_50a5_R2adj.pdf", width = 8, height = 4)

levels_order <- c("DGT","CHM","MAG","YEM")
levels_order <- c("L00", "TB0")

library(dplyr)
library(tidyr)
library(purrr)

fit_ssasymp <- possibly(
  function(df) {
    nls(mean_R2_adj ~ SSasymp(mid, Asym, R0, lrc), data = df)
  },
  otherwise = NULL
)

half_decay <- final_data_clean %>%
  group_by(Source) %>%
  nest() %>%
  mutate(
    model = map(data, fit_ssasymp),
    coef = map(model, ~ if (is.null(.x)) NULL else coef(.x)),
    pred = map2(model, data, ~ if (is.null(.x)) rep(NA_real_, nrow(.y)) else predict(.x, newdata = .y)),
    max_model_R2 = map_dbl(pred, ~ max(.x[-1], na.rm = TRUE)),
    Asym = map_dbl(coef, ~ if (is.null(.x)) NA_real_ else .x[["Asym"]]),
    R0   = map_dbl(coef, ~ if (is.null(.x)) NA_real_ else .x[["R0"]]),
    lrc  = map_dbl(coef, ~ if (is.null(.x)) NA_real_ else .x[["lrc"]]),
    basal = Asym,
    y_half = ((max_model_R2 - basal) / 2) + basal,
    x_half = log((R0 - Asym) / (y_half - Asym)) / exp(lrc)
  ) %>%
  dplyr::select(Source, basal, max_model_R2, y_half, x_half)

final_data_clean$Source <- factor(final_data_clean$Source, levels = levels_order)
half_decay$Source <- factor(half_decay$Source, levels = levels_order)

ggplot(data = final_data_clean, aes(x = mid, y = mean_R2_adj, colour = Source)) +
  geom_point(size = 3, shape = 21, fill = NA) +
  labs(x = "Distance (Kb)", y = expression(LD ~ (r^{2})), colour = "Population") +
  scale_x_continuous(limits = c(0, 50000), breaks = c(0, 1*10^4, 2*10^4, 3*10^4, 4*10^4, 5*10^4),
                     labels = c("0", "10", "20", "30", "40", "50")) +
  facet_wrap(. ~ Source, ncol = 5, nrow = 3) +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.8),
        
        axis.line = element_line(size = 0.8, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8)
  ) +
  geom_smooth(
    method = "nls",
    formula = y ~ SSasymp(x, Asym, R0, lrc),
    se = FALSE,
    aes(group = Source, color = Source)
  ) +
  geom_vline(
    data = half_decay,
    aes(xintercept = x_half),
    colour = "blue",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4)))

ggsave("LD_distance_LATB_50a5_meanR2_adj_V2.pdf", width = 5, height = 4)

ggplot(data = final_data_clean, aes(x = mid, y = mean_R2_adj, colour = Source)) +
  geom_point(size = 3, shape = 21, fill = NA) +
  scale_x_continuous(limits = c(0, 50000), breaks = c(0, 1*10^4, 2*10^4, 3*10^4, 4*10^4, 5*10^4),
                     labels = c("0", "10", "20", "30", "40", "50")) +
  labs(x = "Distance (Kb)", y = expression(LD ~ (r^{2})), colour = "Population") +
  
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.25, "cm"),
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_smooth(
    method = "nls",
    formula = y ~ SSasymp(x, Asym, R0, lrc),
    se = FALSE,
    aes(group = Source)
  ) +
  geom_vline(
    data = half_decay,
    aes(xintercept = x_half),
    colour = "blue",
    linetype = "dashed",
    inherit.aes = FALSE
  )

ggsave("LD_distance_SAU_50a5.pdf", width = 5, height = 4)
